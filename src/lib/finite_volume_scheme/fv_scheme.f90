module mod_finite_volume_schemes

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
  use mod_globals, only: debug_print
  use mod_units
  use mod_floating_point_utils, only: near_zero
  use mod_input, only: input_t
  use mod_reconstruction_factory, only: reconstruction_factory
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_flux_tensor, only: H => flux_tensor_t
  use mod_grid_factory, only: grid_factory
  use mod_bc_factory, only: bc_factory
  use mod_eos, only: eos
  use mod_source, only: source_t
  use mod_source_factory, only: source_factory
  use mod_evo_operator_factory, only: evo_operator_factory
  use mod_grid, only: grid_t
  use hdf5_interface, only: hdf5_file
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_nondimensionalization, only: set_scale_factors
  use mod_fluid, only: fluid_t

  implicit none
  private
  public :: finite_volume_scheme_t, make_fv_scheme

  type :: finite_volume_scheme_t
    !< Abstract representation of the finite volume scheme. This is essentially a puppeteer that
    !< manages the reconstruction, grid, and evolution operator to ultimately calculate the conserved
    !< state variables of each finite cell. The reconstruction, grid, and evolution implementations are passed
    !< on to decendants like the FVLEG scheme.

    character(len=32) :: title = ''
    integer(ik) :: iteration = 0
    real(rk) :: delta_t = 0.0_rk
    real(rk) :: time = 0.0_rk
    integer(ik) :: error_code = 0

    ! class(abstract_reconstruction_t), allocatable :: reconstruction_operator
    ! !< R_Omega reconstruction operator used to reconstruct the corners/midpoints based on the cell
    ! !< average (and gradient if high(er) order reconstruction used)

    ! class(abstract_evo_operator_t), allocatable :: evolution_operator
    ! !< Evolution operator to construct (rho, u, v, p) at each corner and midpoint

    class(grid_t), allocatable :: grid
    !< Grid class to hold geometric information (edge lengths, volumes, etc.)

    class(source_t), allocatable :: source_term
    class(fluid_t), allocatable :: fluid

  contains
    procedure, public :: initialize
    procedure, public :: reconstruct
    procedure, public :: apply_primitive_vars_bc
    procedure, public :: apply_reconstructed_state_bc
    procedure, public :: apply_gradient_bc
    procedure, public :: apply_source_terms
    procedure, public :: set_time
    final :: finalize
  end type

  interface make_fv_scheme
    module procedure :: constructor
  end interface

contains

  function constructor(input) result(fv)
    class(input_t), intent(in) :: input
    type(finite_volume_scheme_t), pointer :: fv

    allocate(fv)
    call fv%initialize(input)

  end function constructor

  subroutine initialize(self, input)
    !< Construct the finite volume local evolution Galerkin (fvleg) scheme
    class(finite_volume_scheme_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    ! class(boundary_condition_t), pointer :: bc => null()
    class(source_t), pointer :: source_term => null()
    class(grid_t), pointer :: grid => null()
    class(abstract_reconstruction_t), pointer :: r_omega => null()
    class(abstract_evo_operator_t), pointer :: E0 => null()
    type(hdf5_file) :: h5
    integer(ik) :: alloc_status
    integer(ik), dimension(:, :), allocatable :: ghost_layers
    character(32) :: str_buff
    alloc_status = 0

    call debug_print('Initializing finite_volume_scheme_t', __FILE__, __LINE__)

    if(input%restart_from_file) then
      call h5%initialize(filename=trim(input%restart_file), status='old', action='r')
      call h5%get('/time', self%time)
      call h5%get('/iteration', self%iteration)

      call h5%readattr('/time', 'units', str_buff)
      select case(trim(str_buff))
      case('ns')
        self%time = self%time * ns_to_s
      case('s')
        ! Do nothing, since seconds is what the code works in
      case default
        error stop "Unknown time units in .h5 file. Acceptable units are 'ns' or 's'."
      end select
      call h5%finalize()
    end if

    self%title = trim(input%title)

    grid => grid_factory(input)
    allocate(self%grid, source=grid, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%grid"
    deallocate(grid)

    ! Now that the grid is set, we can finish setting up the scale factors
    ! for non-dimensionalization. The grid sets the length scale automatically
    ! so that the smallest cell edge length is 1.
    call set_scale_factors(pressure_scale=input%reference_pressure, &
                           density_scale=input%reference_density)

    ! allocate(ghost_layers(4, 2))

    ! ! The ghost_layers array just tells the boundary condition which cell indices
    ! ! are tagged as boundaries. See boundary_condition_t%set_indices() for details
    ! ghost_layers(1, :) = [self%grid%ilo_bc_cell, self%grid%ilo_cell - 1] ! ilo
    ! ghost_layers(3, :) = [self%grid%jlo_bc_cell, self%grid%jlo_cell - 1] ! jlo

    ! ghost_layers(2, :) = [self%grid%ihi_cell + 1, self%grid%ihi_bc_cell] ! ihi
    ! ghost_layers(4, :) = [self%grid%jhi_cell + 1, self%grid%jhi_bc_cell] ! jhi

    ! ! Set boundary conditions
    ! bc => bc_factory(bc_type=input%plus_x_bc, location='+x', input=input, ghost_layers=ghost_layers)
    ! allocate(self%bc_plus_x, source=bc, stat=alloc_status)
    ! if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%bc_plus_x"
    ! deallocate(bc)

    ! bc => bc_factory(bc_type=input%plus_y_bc, location='+y', input=input, ghost_layers=ghost_layers)
    ! allocate(self%bc_plus_y, source=bc, stat=alloc_status)
    ! if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%bc_plus_y"
    ! deallocate(bc)

    ! bc => bc_factory(bc_type=input%minus_x_bc, location='-x', input=input, ghost_layers=ghost_layers)
    ! allocate(self%bc_minus_x, source=bc, stat=alloc_status)
    ! if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%bc_minus_x"
    ! deallocate(bc)

    ! bc => bc_factory(bc_type=input%minus_y_bc, location='-y', input=input, ghost_layers=ghost_layers)
    ! allocate(self%bc_minus_y, source=bc, stat=alloc_status)
    ! if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%bc_minus_y"
    ! deallocate(bc)

    if(input%enable_source_terms) then
      source_term => source_factory(input)
      allocate(self%source_term, source=source_term, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%source_term"
    end if

    write(*, '(a)') "Boundary Conditions"
    write(*, '(a)') "==================="
    write(*, '(3(a),i0,a)') "+x: ", trim(self%bc_plus_x%name), ' (priority = ', self%bc_plus_x%priority, ')'
    write(*, '(3(a),i0,a)') "-x: ", trim(self%bc_minus_x%name), ' (priority = ', self%bc_minus_x%priority, ')'
    write(*, '(3(a),i0,a)') "+y: ", trim(self%bc_plus_y%name), ' (priority = ', self%bc_plus_y%priority, ')'
    write(*, '(3(a),i0,a)') "-y: ", trim(self%bc_minus_y%name), ' (priority = ', self%bc_minus_y%priority, ')'
    write(*, *)

    r_omega => reconstruction_factory(input=input, grid_target=self%grid)
    allocate(self%reconstruction_operator, source=r_omega, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%reconstruction_operator"
    deallocate(r_omega)

    call self%reconstruction_operator%set_grid_pointer(self%grid)

    call debug_print('Making an E0 operator', __FILE__, __LINE__)
    E0 => evo_operator_factory(input=input, grid_target=self%grid, &
                               recon_operator_target=self%reconstruction_operator)

    allocate(self%evolution_operator, source=E0, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%evolution_operator"
    deallocate(E0)

    call self%evolution_operator%set_grid_pointer(grid_target=self%grid)
    call self%evolution_operator%set_reconstruction_operator_pointer(operator_target=self%reconstruction_operator)
  end subroutine initialize

  subroutine finalize(self)
    !< Implementation of the class cleanup
    type(finite_volume_scheme_t), intent(inout) :: self
    integer(ik) :: alloc_status

    alloc_status = 0

    call debug_print('Running finite_volume_scheme_t%finalize()', __FILE__, __LINE__)

    if(allocated(self%evolution_operator)) then
      deallocate(self%evolution_operator, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%evolution_operator"
    end if

    if(allocated(self%grid)) then
      deallocate(self%grid, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%grid"
    end if

    if(allocated(self%reconstruction_operator)) then
      deallocate(self%reconstruction_operator, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%reconstruction_operator"
    end if

    if(allocated(self%source_term)) then
      deallocate(self%source_term, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%source_term"
    end if

    if(allocated(self%bc_plus_x)) then
      deallocate(self%bc_plus_x, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%bc_plus_x"
    end if

    if(allocated(self%bc_plus_y)) then
      deallocate(self%bc_plus_y, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%bc_plus_y"
    end if

    if(allocated(self%bc_minus_x)) then
      deallocate(self%bc_minus_x, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%bc_minus_x"
    end if

    if(allocated(self%bc_minus_y)) then
      deallocate(self%bc_minus_y, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate finite_volume_scheme_t%bc_minus_y"
    end if

  end subroutine finalize

  subroutine apply_source_terms(self, conserved_vars, lbounds)
    class(finite_volume_scheme_t), intent(inout) :: self
    integer(ik), dimension(3), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(inout) :: conserved_vars

    if(allocated(self%source_term)) then

      ! if(self%source_term%ilo /= 0 .and. self%source_term%ihi /= 0) then
      !   if(max(self%source_term%ilo, self%source_term%ihi) > max(self%grid%ilo_cell, self%grid%ihi_cell)) then
      !     error stop "max(self%source_term%ilo, self%source_term%ihi) > max(self%grid%ilo_cell, self%grid%ihi_cell)"
      !   end if

      !   if(min(self%source_term%ilo, self%source_term%ihi) < min(self%grid%ilo_cell, self%grid%ihi_cell)) then
      !     error stop "min(self%source_term%ilo, self%source_term%ihi) > min(self%grid%ilo_cell, self%grid%ihi_cell)"
      !   end if
      ! end if

      ! if(self%source_term%jlo /= 0 .and. self%source_term%jhi /= 0) then
      !   if(max(self%source_term%jlo, self%source_term%jhi) > max(self%grid%jlo_cell, self%grid%jhi_cell)) then
      !     error stop "max(self%source_term%jlo, self%source_term%jhi) > max(self%grid%jlo_cell, self%grid%jhi_cell)"
      !   end if

      !   if(min(self%source_term%jlo, self%source_term%jhi) < min(self%grid%jlo_cell, self%grid%jhi_cell)) then
      !     error stop "min(self%source_term%jlo, self%source_term%jhi) > min(self%grid%jlo_cell, self%grid%jhi_cell)"
      !   end if
      ! end if

      call self%source_term%apply_source(conserved_vars=conserved_vars, lbounds=lbounds, time=self%time)
    end if

  end subroutine apply_source_terms

  subroutine reconstruct(self, primitive_var, lbounds, reconstructed_var, stage, name)
    !< Implementation of the FVLEG reconstruction.
    !< This reconstructs the entire grid at all the nodes/midpoints
    class(finite_volume_scheme_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    integer(ik), intent(in) :: stage !< stage in the time integration scheme, e.g. stage 2 in RK2
    character(len=*), intent(in) :: name !< what variable will be reconstructed
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in) :: primitive_var !< (i,j); cell primitive variable to reconstruct
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(out) :: reconstructed_var
    !< ((corner1:midpoint4), i, j); reconstructed variable, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint

    character(len=50) :: stage_name = ''

    write(stage_name, '(2(a, i0))') 'iter_', self%iteration, 'stage_', stage

    call debug_print('Running finite_volume_scheme_t%reconstruct()', __FILE__, __LINE__)
    call self%reconstruction_operator%reconstruct(primitive_var=primitive_var, &
                                                  reconstructed_var=reconstructed_var, &
                                                  lbounds=lbounds, &
                                                  name=name, &
                                                  stage_name=stage_name)
  end subroutine reconstruct

  subroutine set_time(self, time, delta_t, iteration)
    !< Set the time statistics
    class(finite_volume_scheme_t), intent(inout) :: self
    real(rk), intent(in) :: time
    real(rk), intent(in) :: delta_t
    integer(ik), intent(in) :: iteration

    self%time = time
    self%iteration = iteration
    self%delta_t = delta_t
    call self%bc_plus_x%set_time(time)
    call self%bc_minus_x%set_time(time)
    call self%bc_plus_y%set_time(time)
    call self%bc_minus_y%set_time(time)

    self%evolution_operator%time = time
    self%evolution_operator%time_step = delta_t
    self%evolution_operator%iteration = iteration

  end subroutine set_time

  subroutine apply_primitive_vars_bc(self, rho, u, v, p, lbounds)
    !< Apply the boundary conditions
    class(finite_volume_scheme_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: rho
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: u
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: v
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: p

    integer(ik) :: priority
    integer(ik) :: max_priority_bc !< highest goes first

    call debug_print('Running apply_primitive_var_bc', __FILE__, __LINE__)

    max_priority_bc = max(self%bc_plus_x%priority, self%bc_plus_y%priority, &
                          self%bc_minus_x%priority, self%bc_minus_y%priority)

    do priority = max_priority_bc, 0, -1

      if(self%bc_plus_x%priority == priority) then
        call self%bc_plus_x%apply_primitive_var_bc(rho=rho, u=u, v=v, p=p, lbounds=lbounds)
      end if

      if(self%bc_plus_y%priority == priority) then
        call self%bc_plus_y%apply_primitive_var_bc(rho=rho, u=u, v=v, p=p, lbounds=lbounds)
      end if

      if(self%bc_minus_x%priority == priority) then
        call self%bc_minus_x%apply_primitive_var_bc(rho=rho, u=u, v=v, p=p, lbounds=lbounds)
      end if

      if(self%bc_minus_y%priority == priority) then
        call self%bc_minus_y%apply_primitive_var_bc(rho=rho, u=u, v=v, p=p, lbounds=lbounds)
      end if

    end do

  end subroutine apply_primitive_vars_bc

  subroutine apply_gradient_bc(self)
    !< Apply the boundary conditions
    class(finite_volume_scheme_t), intent(inout) :: self

    integer(ik) :: priority
    integer(ik) :: max_priority_bc !< highest goes first
    integer(ik), dimension(2) :: lbounds

    call debug_print('Running finite_volume_scheme_t%apply_gradient_bc()', __FILE__, __LINE__)

    lbounds = lbound(self%reconstruction_operator%grad_x_rho)
    max_priority_bc = max(self%bc_plus_x%priority, self%bc_plus_y%priority, &
                          self%bc_minus_x%priority, self%bc_minus_y%priority)

    do priority = max_priority_bc, 0, -1

      if(self%bc_plus_x%priority == priority) then
        call self%bc_plus_x%apply_gradient_bc(grad_x=self%reconstruction_operator%grad_x_rho, &
                                              grad_y=self%reconstruction_operator%grad_y_rho, &
                                              lbounds=lbounds)
        call self%bc_plus_x%apply_gradient_bc(grad_x=self%reconstruction_operator%grad_x_p, &
                                              grad_y=self%reconstruction_operator%grad_y_p, &
                                              lbounds=lbounds)
      end if

      if(self%bc_plus_y%priority == priority) then
        call self%bc_plus_y%apply_gradient_bc(grad_x=self%reconstruction_operator%grad_x_rho, &
                                              grad_y=self%reconstruction_operator%grad_y_rho, &
                                              lbounds=lbounds)
        call self%bc_plus_y%apply_gradient_bc(grad_x=self%reconstruction_operator%grad_x_p, &
                                              grad_y=self%reconstruction_operator%grad_y_p, &
                                              lbounds=lbounds)
      end if

      if(self%bc_minus_x%priority == priority) then
        call self%bc_minus_x%apply_gradient_bc(grad_x=self%reconstruction_operator%grad_x_rho, &
                                               grad_y=self%reconstruction_operator%grad_y_rho, &
                                               lbounds=lbounds)
        call self%bc_minus_x%apply_gradient_bc(grad_x=self%reconstruction_operator%grad_x_p, &
                                               grad_y=self%reconstruction_operator%grad_y_p, &
                                               lbounds=lbounds)
      end if

      if(self%bc_minus_y%priority == priority) then
        call self%bc_minus_y%apply_gradient_bc(grad_x=self%reconstruction_operator%grad_x_rho, &
                                               grad_y=self%reconstruction_operator%grad_y_rho, &
                                               lbounds=lbounds)
        call self%bc_minus_y%apply_gradient_bc(grad_x=self%reconstruction_operator%grad_x_p, &
                                               grad_y=self%reconstruction_operator%grad_y_p, &
                                               lbounds=lbounds)
      end if

    end do

  end subroutine apply_gradient_bc

  subroutine apply_reconstructed_state_bc(self, recon_rho, recon_u, recon_v, recon_p, lbounds)
    !< Apply the boundary conditions
    class(finite_volume_scheme_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_rho
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_u
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_v
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_p
    integer(ik) :: priority
    integer(ik) :: max_priority_bc !< highest goes first

    call debug_print('Running apply_reconstructed_state_bc', __FILE__, __LINE__)

    max_priority_bc = max(self%bc_plus_x%priority, self%bc_plus_y%priority, &
                          self%bc_minus_x%priority, self%bc_minus_y%priority)

    do priority = max_priority_bc, 0, -1

      if(self%bc_plus_x%priority == priority) then
        call self%bc_plus_x%apply_reconstructed_state_bc(recon_rho=recon_rho, &
                                                         recon_u=recon_u, &
                                                         recon_v=recon_v, &
                                                         recon_p=recon_p, lbounds=lbounds)
      end if

      if(self%bc_plus_y%priority == priority) then
        call self%bc_plus_y%apply_reconstructed_state_bc(recon_rho=recon_rho, &
                                                         recon_u=recon_u, &
                                                         recon_v=recon_v, &
                                                         recon_p=recon_p, lbounds=lbounds)
      end if

      if(self%bc_minus_x%priority == priority) then
        call self%bc_minus_x%apply_reconstructed_state_bc(recon_rho=recon_rho, &
                                                          recon_u=recon_u, &
                                                          recon_v=recon_v, &
                                                          recon_p=recon_p, lbounds=lbounds)
      end if

      if(self%bc_minus_y%priority == priority) then
        call self%bc_minus_y%apply_reconstructed_state_bc(recon_rho=recon_rho, &
                                                          recon_u=recon_u, &
                                                          recon_v=recon_v, &
                                                          recon_p=recon_p, lbounds=lbounds)
      end if

    end do
  end subroutine apply_reconstructed_state_bc

end module mod_finite_volume_schemes
