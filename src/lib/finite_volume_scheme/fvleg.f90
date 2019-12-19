module mod_fvleg

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use mod_globals, only: debug_print
  use mod_floating_point_utils, only: near_zero
  use mod_surrogate, only: surrogate
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t
  use mod_input, only: input_t
  use mod_finite_volume_schemes, only: finite_volume_scheme_t
  use mod_reconstruction_factory, only: reconstruction_factory
  ! use mod_local_evo_operator, only: local_evo_operator_t
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_second_order_reconstruction, only: second_order_reconstruction_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_flux_tensor, only: H => flux_tensor_t
  use mod_grid_factory, only: grid_factory
  use mod_bc_factory, only: bc_factory
  use mod_evo_operator_factory, only: evo_operator_factory
  use mod_grid, only: grid_t
  use hdf5_interface, only: hdf5_file
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_time_integrator_factory, only: time_integrator_factory

  ! use mod_first_order_reconstruction, only: first_order_reconstruction_t

  implicit none
  private
  public :: fvleg_t, new_fvleg

  type, extends(finite_volume_scheme_t) :: fvleg_t
    !< Implementation of the finite volume local evolution Galerkin (FVLEG) scheme type

  contains
    procedure, public :: initialize
    procedure, public :: reconstruct => reconstruct_fvleg
    procedure, public :: evolve_domain
    procedure, public :: apply_conserved_vars_bc
    procedure, public :: apply_reconstructed_state_bc
    procedure, public :: apply_cell_gradient_bc
    procedure, public :: t => time_derivative
    procedure, private :: integrate_fluxes
    procedure, private :: initialize_from_hdf5
    procedure, private :: initialize_from_ini
    procedure, pass(lhs), public :: type_plus_type => add_fvleg
    procedure, pass(lhs), public :: type_minus_type => subtract_fvleg
    procedure, pass(lhs), public :: type_mul_real => fvleg_mul_real
    procedure, pass(rhs), public :: real_mul_type => real_mul_fvleg
    procedure, pass(lhs), public :: assign => assign_fvleg
    procedure, public :: force_finalization
    final :: finalize
  end type

  interface new_fvleg
    module procedure :: constructor
  end interface

contains

  function constructor(input) result(fvleg)
    class(input_t), intent(in) :: input
    type(fvleg_t), pointer :: fvleg

    allocate(fvleg)
    call fvleg%initialize(input)

  end function

  subroutine initialize(self, input)
    !< Construct the finite volume local evolution Galerkin (fvleg) scheme
    class(fvleg_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    class(boundary_condition_t), pointer :: bc => null()
    class(grid_t), pointer :: grid => null()
    class(strategy), pointer :: time_integrator => null()
    class(abstract_reconstruction_t), pointer :: r_omega => null()
    class(abstract_evo_operator_t), pointer :: E0 => null()

    integer(ik) :: alloc_status
    alloc_status = 0

    call debug_print('Initializing fvleg_t', __FILE__, __LINE__)

    self%title = input%title

    grid => grid_factory(input)
    allocate(self%grid, source=grid, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%grid"
    deallocate(grid)

    time_integrator => time_integrator_factory(input)
    allocate(self%time_integrator, source=time_integrator, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%time_integrator"
    deallocate(time_integrator)

    ! Set boundary conditions
    bc => bc_factory(bc_type=input%plus_x_bc, location='+x')
    allocate(self%bc_plus_x, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%bc_plus_x"
    deallocate(bc)

    bc => bc_factory(bc_type=input%plus_y_bc, location='+y')
    allocate(self%bc_plus_y, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%bc_plus_y"
    deallocate(bc)

    bc => bc_factory(bc_type=input%minus_x_bc, location='-x')
    allocate(self%bc_minus_x, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%bc_minus_x"
    deallocate(bc)

    bc => bc_factory(bc_type=input%minus_y_bc, location='-y')
    allocate(self%bc_minus_y, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%bc_minus_y"
    deallocate(bc)

    associate(imin=>self%grid%ilo_bc_cell, imax=>self%grid%ihi_bc_cell, &
              jmin=>self%grid%jlo_bc_cell, jmax=>self%grid%jhi_bc_cell)

      allocate(self%conserved_vars(4, imin:imax, jmin:jmax), stat=alloc_status)
      ! ((rho,u,v,p),i,j) Conserved variables for each cell
      if(alloc_status /= 0) then
        error stop "Unable to allocate fvleg_t%conserved_vars"
      end if

      allocate(self%reconstructed_state(4, 4, 2, imin:imax, jmin:jmax), stat=alloc_status)
      ! ((rho, u ,v, p), point, node/midpoint, i, j); this is a cell-based value, so imax=ni-1, etc
      if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%reconstructed_state"
    end associate

    r_omega => reconstruction_factory(input=input, grid_target=self%grid, &
                                      conserved_vars_target=self%conserved_vars, lbounds=lbound(self%conserved_vars))
    allocate(self%reconstruction_operator, source=r_omega, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%reconstruction_operator"
    deallocate(r_omega)

    call debug_print('Making an E0 operator', __FILE__, __LINE__)
    E0 => evo_operator_factory(input=input, grid_target=self%grid, &
                               recon_operator_target=self%reconstruction_operator, &
                               reconstructed_state_target=self%reconstructed_state, lbounds=lbound(self%reconstructed_state))

    allocate(self%evolution_operator, source=E0, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%evolution_operator"
    deallocate(E0)

    associate(imin_node=>self%grid%ilo_node, imax_node=>self%grid%ihi_node, &
              jmin_node=>self%grid%jlo_node, jmax_node=>self%grid%jhi_node, &
              imin_cell=>self%grid%ilo_node, imax_cell=>self%grid%ihi_node, &
              jmin_cell=>self%grid%jlo_node, jmax_cell=>self%grid%jhi_node)

      ! corners
      if(.not. allocated(self%evolved_corner_state)) then
        ! ((rho,u,v,p), i, j); Reconstructed U at each corner
        allocate(self%evolved_corner_state(4, imin_node:imax_node, jmin_node:jmax_node), stat=alloc_status)
        if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%evolved_corner_state"
      end if

      if(.not. allocated(self%corner_reference_state)) then
        ! ((rho, u ,v, p), i, j); Reference state (tilde) at each corner
        allocate(self%corner_reference_state(4, imin_node:imax_node, jmin_node:jmax_node), stat=alloc_status)
        if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%corner_reference_state"
      endif

      ! left/right midpoints
      if(.not. allocated(self%evolved_leftright_midpoints_state)) then
        ! ((rho,u,v,p), i, j); Reconstructed U at each midpoint on the left/right edges (edges 1 and 3)
        allocate(self%evolved_leftright_midpoints_state(4, imin_cell:imax_cell, jmin_node:jmax_node), stat=alloc_status)
        if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%evolved_leftright_midpoints_state"
      end if

      if(.not. allocated(self%leftright_midpoints_reference_state)) then
        ! ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the left/right edges (edges 1 and 3)
        allocate(self%leftright_midpoints_reference_state(4, imin_cell:imax_cell, jmin_node:jmax_node), stat=alloc_status)
        if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%leftright_midpoints_reference_state"
      end if

      ! down/up midpoints
      if(.not. allocated(self%evolved_downup_midpoints_state)) then
        ! ((rho,u,v,p), i, j); Reconstructed U at each midpoint on the down/up edges (edges 2 and 4)
        allocate(self%evolved_downup_midpoints_state(4, imin_node:imax_node, jmin_cell:jmax_cell), stat=alloc_status)
        if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%evolved_downup_midpoints_state"
      end if

      if(.not. allocated(self%downup_midpoints_reference_state)) then
        ! ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the down/up edges (edges 2 and 4)
        allocate(self%downup_midpoints_reference_state(4, imin_node:imax_node, jmin_cell:jmax_cell), stat=alloc_status)
        if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%downup_midpoints_reference_state"
      end if

    end associate

    if(input%read_init_cond_from_file) then
      call self%initialize_from_hdf5(input)
    else
      call self%initialize_from_ini(input)
    end if

    self%initiated = .true.
  end subroutine initialize

  subroutine initialize_from_hdf5(self, input)
    !< Initialize from an .hdf5 file. The conserved variables are already allocated appropriately from
    !< from the grid class, but this will just initialize them to the values found in the hdf5 file
    class(fvleg_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    type(hdf5_file) :: h5

    real(rk), dimension(:, :), allocatable :: density
    real(rk), dimension(:, :), allocatable :: x_velocity
    real(rk), dimension(:, :), allocatable :: y_velocity
    real(rk), dimension(:, :), allocatable :: pressure

    call debug_print('Initializing fvleg_t from hdf5', __FILE__, __LINE__)
    call h5%initialize(filename=input%initial_condition_file, status='old', action='r')
    call h5%get('/density', density)
    call h5%get('/x_velocity', x_velocity)
    call h5%get('/y_velocity', y_velocity)
    call h5%get('/pressure', pressure)
    call h5%finalize()

    if(any(near_zero(pressure))) then
      error stop "Some (or all) of the pressure array is ~0 in fvleg_t%initialize_from_hdf5"
    end if

    associate(imin=>self%grid%ilo_bc_cell, imax=>self%grid%ihi_bc_cell, &
              jmin=>self%grid%jlo_bc_cell, jmax=>self%grid%jhi_bc_cell)

      self%conserved_vars(1, imin:imax, jmin:jmax) = density    ! (1:imax,1:jmax)
      self%conserved_vars(2, imin:imax, jmin:jmax) = x_velocity ! (1:imax,1:jmax)
      self%conserved_vars(3, imin:imax, jmin:jmax) = y_velocity ! (1:imax,1:jmax)
      self%conserved_vars(4, imin:imax, jmin:jmax) = pressure   ! (1:imax,1:jmax)
    end associate

  end subroutine initialize_from_hdf5

  subroutine initialize_from_ini(self, input)
    !< Initialize from an .ini file. The conserved variables are already allocated appropriately from
    !< from the grid class, but this will just initialize them to the values found in the .ini file
    class(fvleg_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    call debug_print('Initializing fvleg_t from .ini', __FILE__, __LINE__)
    write(*, '(a,4(f0.3, 1x))') 'Initializing fvleg_t%conserved_vars to [rho,u,v,p]: ', &
      input%init_density, input%init_x_velocity, input%init_y_velocity, input%init_pressure
    self%conserved_vars(1, :, :) = input%init_density
    self%conserved_vars(2, :, :) = input%init_x_velocity
    self%conserved_vars(3, :, :) = input%init_y_velocity
    self%conserved_vars(4, :, :) = input%init_pressure

  end subroutine initialize_from_ini

  subroutine force_finalization(self)
    class(fvleg_t), intent(inout) :: self
    integer(ik) :: alloc_status

    call debug_print('Calling fvleg_t%force_finalization()', __FILE__, __LINE__)
    if(allocated(self%evolution_operator)) then
      deallocate(self%evolution_operator, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%evolution_operator"
    end if

    if(allocated(self%grid)) then
      deallocate(self%grid, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%grid"
    end if

    if(allocated(self%reconstruction_operator)) then
      deallocate(self%reconstruction_operator, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%reconstruction_operator"
    end if

    if(allocated(self%bc_plus_x)) then
      deallocate(self%bc_plus_x, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%bc_plus_x"
    end if

    if(allocated(self%bc_plus_y)) then
      deallocate(self%bc_plus_y, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%bc_plus_y"
    end if

    if(allocated(self%bc_minus_x)) then
      deallocate(self%bc_minus_x, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%bc_minus_x"
    end if

    if(allocated(self%bc_minus_y)) then
      deallocate(self%bc_minus_y, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%bc_minus_y"
    end if

    if(allocated(self%conserved_vars)) then
      deallocate(self%conserved_vars, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%conserved_vars"
    end if

    if(allocated(self%reconstructed_state)) then
      deallocate(self%reconstructed_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%reconstructed_state"
    end if

    if(allocated(self%evolved_corner_state)) then
      deallocate(self%evolved_corner_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%evolved_corner_state"
    end if

    if(allocated(self%corner_reference_state)) then
      deallocate(self%corner_reference_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%corner_reference_state"
    end if

    if(allocated(self%evolved_downup_midpoints_state)) then
      deallocate(self%evolved_downup_midpoints_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%evolved_downup_midpoints_state"
    end if

    if(allocated(self%evolved_leftright_midpoints_state)) then
      deallocate(self%evolved_leftright_midpoints_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%evolved_leftright_midpoints_state"
    end if

    if(allocated(self%downup_midpoints_reference_state)) then
      deallocate(self%downup_midpoints_reference_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%downup_midpoints_reference_state"
    end if

    if(allocated(self%leftright_midpoints_reference_state)) then
      deallocate(self%leftright_midpoints_reference_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%leftright_midpoints_reference_state"
    end if
  end subroutine

  subroutine finalize(self)
    !< Implementation of the class cleanup
    type(fvleg_t), intent(inout) :: self
    integer(ik) :: alloc_status

    alloc_status = 0

    call debug_print('Calling fvleg_t%finalize()', __FILE__, __LINE__)

    if(allocated(self%evolution_operator)) then
      deallocate(self%evolution_operator, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%evolution_operator"
    end if

    if(allocated(self%grid)) then
      deallocate(self%grid, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%grid"
    end if

    if(allocated(self%reconstruction_operator)) then
      deallocate(self%reconstruction_operator, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%reconstruction_operator"
    end if

    if(allocated(self%bc_plus_x)) then
      deallocate(self%bc_plus_x, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%bc_plus_x"
    end if

    if(allocated(self%bc_plus_y)) then
      deallocate(self%bc_plus_y, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%bc_plus_y"
    end if

    if(allocated(self%bc_minus_x)) then
      deallocate(self%bc_minus_x, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%bc_minus_x"
    end if

    if(allocated(self%bc_minus_y)) then
      deallocate(self%bc_minus_y, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%bc_minus_y"
    end if

    if(allocated(self%conserved_vars)) then
      deallocate(self%conserved_vars, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%conserved_vars"
    end if

    if(allocated(self%reconstructed_state)) then
      deallocate(self%reconstructed_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%reconstructed_state"
    end if

    if(allocated(self%evolved_corner_state)) then
      deallocate(self%evolved_corner_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%evolved_corner_state"
    end if

    if(allocated(self%corner_reference_state)) then
      deallocate(self%corner_reference_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%corner_reference_state"
    end if

    if(allocated(self%evolved_downup_midpoints_state)) then
      deallocate(self%evolved_downup_midpoints_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%evolved_downup_midpoints_state"
    end if

    if(allocated(self%evolved_leftright_midpoints_state)) then
      deallocate(self%evolved_leftright_midpoints_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%evolved_leftright_midpoints_state"
    end if

    if(allocated(self%downup_midpoints_reference_state)) then
      deallocate(self%downup_midpoints_reference_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%downup_midpoints_reference_state"
    end if

    if(allocated(self%leftright_midpoints_reference_state)) then
      deallocate(self%leftright_midpoints_reference_state, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%leftright_midpoints_reference_state"
    end if

  end subroutine finalize

  subroutine reconstruct_fvleg(self)
    !< Implementation of the FVLEG reconstruction.
    !< This reconstructs the entire grid at all the nodes/midpoints
    class(fvleg_t), intent(inout) :: self

    call debug_print('Reconstructing the domain', __FILE__, __LINE__)
    call self%reconstruction_operator%set_conserved_vars_pointer(conserved_vars=self%conserved_vars, &
                                                                 lbounds=lbound(self%conserved_vars))
    call self%reconstruction_operator%set_grid_pointer(self%grid)
    call self%reconstruction_operator%reconstruct_domain(reconstructed_domain=self%reconstructed_state, lbounds=lbound(self%reconstructed_state))
  end subroutine reconstruct_fvleg

  subroutine apply_conserved_vars_bc(self)
    !< Apply the boundary conditions
    class(fvleg_t), intent(inout) :: self

    call debug_print('Calling apply_conserved_var_bc', __FILE__, __LINE__)
    call self%bc_plus_x%apply_conserved_var_bc(conserved_vars=self%conserved_vars)
    call self%bc_plus_y%apply_conserved_var_bc(conserved_vars=self%conserved_vars)
    call self%bc_minus_x%apply_conserved_var_bc(conserved_vars=self%conserved_vars)
    call self%bc_minus_y%apply_conserved_var_bc(conserved_vars=self%conserved_vars)
  end subroutine apply_conserved_vars_bc

  subroutine apply_reconstructed_state_bc(self)
    !< Apply the boundary conditions
    class(fvleg_t), intent(inout) :: self

    call debug_print('Calling apply_reconstructed_state_bc', __FILE__, __LINE__)
    call self%bc_plus_x%apply_reconstructed_state_bc(reconstructed_state=self%reconstructed_state)
    call self%bc_plus_y%apply_reconstructed_state_bc(reconstructed_state=self%reconstructed_state)
    call self%bc_minus_x%apply_reconstructed_state_bc(reconstructed_state=self%reconstructed_state)
    call self%bc_minus_y%apply_reconstructed_state_bc(reconstructed_state=self%reconstructed_state)
  end subroutine apply_reconstructed_state_bc

  subroutine apply_cell_gradient_bc(self, r_omega)
    !< Apply the boundary conditions
    class(fvleg_t), intent(inout) :: self
    class(abstract_reconstruction_t), intent(in) :: r_omega
    select type(r_omega)
    class is(second_order_reconstruction_t)
      call debug_print('Calling apply_cell_gradient_bc', __FILE__, __LINE__)
      call self%bc_plus_x%apply_cell_gradient_bc(cell_gradient=self%reconstruction_operator%cell_gradient)
      call self%bc_plus_y%apply_cell_gradient_bc(cell_gradient=self%reconstruction_operator%cell_gradient)
      call self%bc_minus_x%apply_cell_gradient_bc(cell_gradient=self%reconstruction_operator%cell_gradient)
      call self%bc_minus_y%apply_cell_gradient_bc(cell_gradient=self%reconstruction_operator%cell_gradient)
    end select

  end subroutine apply_cell_gradient_bc

  function time_derivative(self) result(dU_dt)
    ! //TODO: make pure
    !< Implementation of dU/dt = 1/Omega_ij Sum(1,k) Integral(H . n dl)
    class(fvleg_t), intent(in) :: self
    class(integrand_t), allocatable :: dU_dt

    type(fvleg_t), allocatable :: local_dU_dt

    call debug_print('Finding dU/dt', __FILE__, __LINE__)

    allocate(local_dU_dt, source=self)

    ! if (allocated(local_dU_dt%conserved_vars)) deallocate(local_dU_dt%conserved_vars)
    ! allocate(local_dU_dt%conserved_vars, source=self%conserved_vars)

    ! if (allocated(local_dU_dt%evolved_corner_state)) deallocate(local_dU_dt%evolved_corner_state)
    ! allocate(local_dU_dt%evolved_corner_state, source=self%evolved_corner_state)

    ! if (allocated(local_dU_dt%corner_reference_state)) deallocate(local_dU_dt%corner_reference_state)
    ! allocate(local_dU_dt%corner_reference_state, source=self%corner_reference_state)

    ! if (allocated(local_dU_dt%evolved_downup_midpoints_state)) deallocate(local_dU_dt%evolved_downup_midpoints_state)
    ! allocate(local_dU_dt%evolved_downup_midpoints_state, source=self%evolved_downup_midpoints_state)

    ! if (allocated(local_dU_dt%evolved_leftright_midpoints_state)) deallocate(local_dU_dt%evolved_leftright_midpoints_state)
    ! allocate(local_dU_dt%evolved_leftright_midpoints_state, source=self%evolved_leftright_midpoints_state)

    ! ! if (allocated()) deallocate()
    ! allocate(local_dU_dt%downup_midpoints_reference_state, source=self%downup_midpoints_reference_state)

    ! ! if (allocated()) deallocate()
    ! allocate(local_dU_dt%leftright_midpoints_reference_state, source=self%leftright_midpoints_reference_state)

    ! ! if (allocated()) deallocate()
    ! allocate(local_dU_dt%reconstructed_state, source=self%reconstructed_state)

    ! if (allocated(local_dU_dt%title)) deallocate(local_dU_dt%title)
    ! allocate(character(len=len(self%title)) :: local_dU_dt%title)
    ! local_dU_dt%title = self%title

    ! ! if (allocated()) deallocate()
    ! allocate(local_dU_dt%reconstruction_operator, source=self%reconstruction_operator)

    ! ! if (allocated()) deallocate()
    ! ! allocate(local_dU_dt%evolution_operator, source=self%evolution_operator)
    ! local_dU_dt%evolution_operator=self%evolution_operator

    ! ! if (allocated(local_dU_dt%grid)) deallocate()
    ! allocate(local_dU_dt%grid, source=self%grid)

    ! if (allocated(local_dU_dt%bc_plus_x)) deallocate(local_dU_dt%bc_plus_x)
    ! allocate(local_dU_dt%bc_plus_x , source=self%bc_plus_x )
    ! ! local_dU_dt%bc_plus_x=self%bc_plus_x

    ! if (allocated(local_dU_dt%bc_plus_y)) deallocate(local_dU_dt%bc_plus_y)
    ! allocate(local_dU_dt%bc_plus_y , source=self%bc_plus_y )

    ! if (allocated(local_dU_dt%bc_minus_x)) deallocate(local_dU_dt%bc_minus_x)
    ! allocate(local_dU_dt%bc_minus_x, source=self%bc_minus_x)

    ! if (allocated(local_dU_dt%bc_minus_y)) deallocate(local_dU_dt%bc_minus_y)
    ! allocate(local_dU_dt%bc_minus_y, source=self%bc_minus_y)

    call local_dU_dt%set_temp(calling_function='time_derivative (local_dU_dt)', line=__LINE__)

    call local_dU_dt%calculate_reference_state()
    call local_dU_dt%apply_conserved_vars_bc()
    call local_dU_dt%reconstruct()                    ! reconstruct only the real domain // TODO: make sure init applies ghost conditions at the start?
    call local_dU_dt%apply_cell_gradient_bc(self%reconstruction_operator)
    call local_dU_dt%apply_reconstructed_state_bc()   ! ghost cells now pick up reconstructed and conserved states
    call local_dU_dt%evolve_domain()                  ! now evolve the flow (with b.c. and reconstructed state info)
    ! call local_dU_dt%reconstruction_operator%nullify_pointer_members()

    ! Set the RHS of Eq. 3 in the main text
    local_dU_dt%conserved_vars = local_dU_dt%integrate_fluxes()

    call move_alloc(local_dU_dt, dU_dt)

  end function time_derivative

  subroutine evolve_domain(self)
    class(fvleg_t), intent(inout) :: self

    call debug_print('Evolving the domain', __FILE__, __LINE__)

    ! Set the pointers for the evolution operator (makes the code w/in it easier)
    call self%evolution_operator%set_grid_pointer(grid_target=self%grid)
    call self%evolution_operator%set_reconstruction_operator_pointer(operator_target=self%reconstruction_operator)
    call self%evolution_operator%set_reconstructed_state_pointer(reconstructed_state_target=self%reconstructed_state, &
                                                                 lbounds=lbound(self%reconstructed_state))

    ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of left/right edge vectors
    call debug_print('Evolving left/right midpoints', __FILE__, __LINE__)
    call self%evolution_operator%evolve_leftright_midpoints(reference_state=self%leftright_midpoints_reference_state, &
                                                            evolved_state=self%evolved_leftright_midpoints_state)

    ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of down/up edge vectors
    call debug_print('Evolving down/up midpoints', __FILE__, __LINE__)
    call self%evolution_operator%evolve_downup_midpoints(reference_state=self%downup_midpoints_reference_state, &
                                                         evolved_state=self%evolved_downup_midpoints_state)

    ! Evolve, i.e. E0(R_omega), at all corner nodes
    call debug_print('Evolving corner nodes', __FILE__, __LINE__)
    call self%evolution_operator%evolve_corners(reference_state=self%corner_reference_state, &
                                                evolved_state=self%evolved_corner_state)

    ! call self%evolution_operator%nullify_pointer_members()
  end subroutine

  function integrate_fluxes(self) result(rhs)
    !< Evaluate the fluxes along the edges. This is equation 13 in the paper
    use mod_flux_tensor, only: operator(.dot.)
    class(fvleg_t), intent(in) :: self
    real(rk), dimension(4) :: edge_flux

    integer(ik) :: ilo, ihi, jlo, jhi
    real(rk), dimension(:, :, :), allocatable :: rhs
    integer(ik) :: i, j

    ilo = self%grid%ilo_cell
    jlo = self%grid%jlo_cell
    ihi = self%grid%ihi_cell
    jhi = self%grid%jhi_cell

    call debug_print('Integrating fluxes', __FILE__, __LINE__)
    allocate(rhs, source=self%conserved_vars)
    rhs = 0.0_rk

    ! do j = jlo, jhi
    !   do i = ilo, ihi
    do concurrent(j=jlo:jhi)
      do concurrent(i=jlo:jhi)

        edge_flux = 0.0_rk
        ! print *, 'i, j:', i, j
        ! Edge 1 (bottom)
        ! print *, 'bottom'
        associate(E0_R_omega_k1=>self%evolved_corner_state(:, i, j), &
                  E0_R_omega_kc=>self%evolved_leftright_midpoints_state(:, i, j), &
                  E0_R_omega_k2=>self%evolved_corner_state(:, i + 1, j), &
                  n_hat=>self%grid%cell_edge_norm_vectors(:, 1, i, j), &
                  delta_l=>self%grid%cell_edge_lengths(1, i, j))

          ! Eq. 13, for edge 1
          edge_flux = edge_flux + &
                      ( &
                      ((H(E0_R_omega_k1) + &
                        H(E0_R_omega_kc) * 4.0_rk + &
                        H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

        end associate

        ! Edge 2 (right)
        ! print *, 'right'
        associate(E0_R_omega_k1=>self%evolved_corner_state(:, i + 1, j), &
                  E0_R_omega_kc=>self%evolved_downup_midpoints_state(:, i + 1, j), &
                  E0_R_omega_k2=>self%evolved_corner_state(:, i + 1, j + 1), &
                  n_hat=>self%grid%cell_edge_norm_vectors(:, 2, i, j), &
                  delta_l=>self%grid%cell_edge_lengths(2, i, j))

          ! Eq. 13, for edge 2
          edge_flux = edge_flux + &
                      ( &
                      ((H(E0_R_omega_k1) + &
                        H(E0_R_omega_kc) * 4.0_rk + &
                        H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

        end associate

        ! Edge 3 (top)
        ! print *, 'top'
        associate(E0_R_omega_k1=>self%evolved_corner_state(:, i + 1, j + 1), &
                  E0_R_omega_kc=>self%evolved_leftright_midpoints_state(:, i, j + 1), &
                  E0_R_omega_k2=>self%evolved_corner_state(:, i, j + 1), &
                  n_hat=>self%grid%cell_edge_norm_vectors(:, 3, i, j), &
                  delta_l=>self%grid%cell_edge_lengths(3, i, j))

          ! Eq. 13, for edge 3
          edge_flux = edge_flux + &
                      ( &
                      ((H(E0_R_omega_k1) + &
                        H(E0_R_omega_kc) * 4.0_rk + &
                        H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

        end associate

        ! Edge 4 (left)
        ! print *, 'left'
        associate(E0_R_omega_k1=>self%evolved_corner_state(:, i, j + 1), &
                  E0_R_omega_kc=>self%evolved_downup_midpoints_state(:, i, j), &
                  E0_R_omega_k2=>self%evolved_corner_state(:, i, j), &
                  n_hat=>self%grid%cell_edge_norm_vectors(:, 4, i, j), &
                  delta_l=>self%grid%cell_edge_lengths(4, i, j))

          ! Eq. 13, for edge 4
          edge_flux = edge_flux + &
                      ( &
                      ((H(E0_R_omega_k1) + &
                        H(E0_R_omega_kc) * 4.0_rk + &
                        H(E0_R_omega_k2)) .dot.n_hat) * (delta_l / 6.0_rk))

        end associate

        rhs(:, i, j) = (-1.0_rk / self%grid%cell_volume(i, j)) * edge_flux

      end do ! i
    end do ! j

  end function integrate_fluxes

  function subtract_fvleg(lhs, rhs) result(difference)
    !< Implementation of the (-) operator for the fvleg type
    class(fvleg_t), intent(in) :: lhs
    class(integrand_t), intent(in) :: rhs

    class(integrand_t), allocatable :: difference
    type(fvleg_t), allocatable :: local_difference

    call debug_print('Calling fvleg_t%subtract_fvleg()', __FILE__, __LINE__)

    ! call lhs%guard_temp(calling_function='subtract_fvleg', line=__LINE__)

    select type(rhs)
    class is(fvleg_t)
      allocate(local_difference, source=lhs)
      local_difference%conserved_vars = lhs%conserved_vars - rhs%conserved_vars
    class default
      error stop 'fvleg_t%subtract_fvleg: unsupported rhs class'
    end select

    call local_difference%set_temp(calling_function='subtract_fvleg', line=__LINE__)
    call move_alloc(local_difference, difference)
    ! call lhs%clean_temp(calling_function='subtract_fvleg', line=__LINE__)
  end function subtract_fvleg

  function add_fvleg(lhs, rhs) result(sum)
    !< Implementation of the (+) operator for the fvleg type
    class(fvleg_t), intent(in) :: lhs
    class(integrand_t), intent(in) :: rhs

    class(integrand_t), allocatable :: sum
    type(fvleg_t), allocatable :: local_sum

    call debug_print('Calling fvleg_t%add_fvleg()', __FILE__, __LINE__)
    ! call lhs%guard_temp(calling_function='add_fvleg (lhs)', line=__LINE__)

    select type(rhs)
    class is(fvleg_t)
      allocate(local_sum, source=lhs)
      local_sum%conserved_vars = lhs%conserved_vars + rhs%conserved_vars
    class default
      error stop 'fvleg_t%add_fvleg: unsupported rhs class'
    end select

    call local_sum%set_temp(calling_function='add_fvleg (local_sum)', line=__LINE__)
    call move_alloc(local_sum, sum)
    ! call lhs%clean_temp(calling_function='add_fvleg (lhs)', line=__LINE__)
  end function add_fvleg

  function fvleg_mul_real(lhs, rhs) result(product)
    !< Implementation of the fvleg_t * real operation
    class(fvleg_t), intent(in) :: lhs
    real(rk), intent(in) :: rhs
    class(integrand_t), allocatable :: product

    type(fvleg_t), allocatable :: local_product
    integer(ik) :: alloc_status

    call debug_print('Calling fvleg_t%fvleg_mul_real()', __FILE__, __LINE__)
    ! call lhs%guard_temp(calling_function='fvleg_mul_real (lhs)', line=__LINE__)
    ! print*, 'allocated(local_product): ', allocated(local_product)
    ! print*, 'allocated(lhs): ', lhs%initiated
    ! print*, 'allocated(lhs%conserved_vars):', allocated(lhs%conserved_vars)

    allocate(local_product, source=lhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate local_product in fvleg_t%fvleg_mul_real"

    local_product%conserved_vars = lhs%conserved_vars * rhs

    call local_product%set_temp(calling_function='fvleg_mul_real (local_product)', line=__LINE__)
    call move_alloc(local_product, product)
    ! call lhs%clean_temp(calling_function='fvleg_mul_real (lhs)', line=__LINE__)
  end function fvleg_mul_real

  function real_mul_fvleg(lhs, rhs) result(product)
    !< Implementation of the real * fvleg_t operation
    class(fvleg_t), intent(in) :: rhs
    real(rk), intent(in) :: lhs
    class(integrand_t), allocatable :: product

    type(fvleg_t), allocatable :: local_product
    integer(ik) :: alloc_status

    call debug_print('Calling fvleg_t%real_mul_fvleg()', __FILE__, __LINE__)
    ! call rhs%guard_temp(calling_function='real_mul_fvleg (rhs)', line=__LINE__)

    ! print *, 'rhs%initiated', rhs%initiated
    allocate(local_product, source=rhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate local_product in fvleg_t%fvleg_mul_real"

    local_product%conserved_vars = rhs%conserved_vars * lhs

    call local_product%set_temp(calling_function='real_mul_fvleg (local_product)', line=__LINE__)
    call move_alloc(local_product, product)
    ! call rhs%clean_temp(calling_function='real_mul_fvleg (rhs)', line=__LINE__)
  end function real_mul_fvleg

  subroutine assign_fvleg(lhs, rhs)
    !< Implementation of the (=) operator for the fvleg type. e.g. lhs = rhs
    class(fvleg_t), intent(inout) :: lhs
    class(integrand_t), intent(in) :: rhs
    integer(ik) :: alloc_status

    alloc_status = 0
    call debug_print('Calling assign_fvleg_t', __FILE__, __LINE__)
    call rhs%guard_temp(calling_function='assign_fvleg (rhs)', line=__LINE__)

    select type(rhs)
    class is(fvleg_t)
      ! allocate(lhs, source=rhs)
      ! lhs = rhs%duplicate()
      ! ! allocate(regular_2d_grid_t :: lhs%grid)
      ! allocate(lhs%grid, source=rhs%grid, stat=alloc_status)
      ! if(alloc_status /= 0) error stop "Unable to allocate lhs%grid from rhs%grid in fvleg_t%assign"
      lhs%grid = rhs%grid

      ! if (allocated(lhs%bc_plus_x)) deallocate(lhs%bc_plus_x)
      ! allocate(lhs%bc_plus_x, source=rhs%bc_plus_x)

      ! if (allocated(lhs%bc_plus_y)) deallocate(lhs%bc_plus_y)
      ! allocate(lhs%bc_plus_y, source=rhs%bc_plus_y)

      ! if (allocated(lhs%bc_minus_x)) deallocate(lhs%bc_minus_x)
      ! allocate(lhs%bc_minus_x, source=rhs%bc_minus_x)

      ! if (allocated(lhs%bc_minus_y)) deallocate(lhs%bc_minus_y)
      ! allocate(lhs%bc_minus_y, source=rhs%bc_minus_y)

      lhs%bc_plus_x = rhs%bc_plus_x
      lhs%bc_plus_y = rhs%bc_plus_y
      lhs%bc_minus_x = rhs%bc_minus_x
      lhs%bc_minus_y = rhs%bc_minus_y

      ! if(.not. allocated(lhs%evolution_operator)) then
      !   allocate(lhs%evolution_operator, source=rhs%evolution_operator, stat=alloc_status)
      !   if(alloc_status /= 0) then
      !     error stop "Unable to allocate lhs%evolution_operator from rhs%evolution_operator in fvleg_t%assign"
      !   end if
      ! end if
      ! call lhs%evolution_operator%set_grid_pointer(lhs%grid)
      lhs%evolution_operator = rhs%evolution_operator

      ! if(.not. allocated(lhs%time_integrator)) then
      !   allocate(lhs%time_integrator, source=rhs%time_integrator, stat=alloc_status)
      !   if(alloc_status /= 0) then
      !     error stop "Unable to allocate lhs%time_integrator from rhs%time_integrator in fvleg_t%assign"
      !   endif
      ! end if
      lhs%time_integrator = rhs%time_integrator

      ! if(.not. allocated(lhs%reconstruction_operator)) then
      !   allocate(lhs%reconstruction_operator, source=rhs%reconstruction_operator, stat=alloc_status)
      !   if(alloc_status /= 0) then
      !     error stop "Unable to allocate lhs%reconstruction_operator "// &
      !       "from rhs%reconstruction_operator in fvleg_t%assign"
      !   endif
      ! end if
      lhs%reconstruction_operator = rhs%reconstruction_operator

      ! if(.not. allocated(lhs%conserved_vars)) then
      !   allocate(lhs%conserved_vars, mold=rhs%conserved_vars, stat=alloc_status)
      !   if(alloc_status /= 0) then
      !     error stop "Unable to allocate lhs%conserved_vars from "// &
      !       "rhs%conserved_vars in fvleg_t%assign"
      !   endif
      ! end if
      lhs%conserved_vars = rhs%conserved_vars

      ! if(.not. allocated(lhs%evolved_corner_state)) then
      !   allocate(lhs%evolved_corner_state, &
      !            mold=rhs%evolved_corner_state, stat=alloc_status)
      !   if(alloc_status /= 0) then
      !     error stop "Unable to allocate lhs%evolved_corner_state "// &
      !       "from rhs%evolved_corner_state in fvleg_t%assign"
      !   endif
      ! endif
      lhs%evolved_corner_state = rhs%evolved_corner_state

      ! if(.not. allocated(lhs%corner_reference_state)) then
      !   allocate(lhs%corner_reference_state, &
      !            mold=rhs%corner_reference_state, stat=alloc_status)
      !   if(alloc_status /= 0) then
      !     error stop "Unable to allocate lhs%corner_reference_state "// &
      !       "from rhs%corner_reference_state in fvleg_t%assign"
      !   endif
      ! endif
      lhs%corner_reference_state = rhs%corner_reference_state

      ! if(.not. allocated(lhs%evolved_downup_midpoints_state)) then
      !   allocate(lhs%evolved_downup_midpoints_state, &
      !            mold=rhs%evolved_downup_midpoints_state, stat=alloc_status)
      !   if(alloc_status /= 0) then
      !     error stop "Unable to allocate lhs%evolved_downup_midpoints_state "// &
      !       "from rhs%evolved_downup_midpoints_state in fvleg_t%assign"
      !   end if
      ! endif
      lhs%evolved_downup_midpoints_state = rhs%evolved_downup_midpoints_state

      ! if(.not. allocated(lhs%evolved_leftright_midpoints_state)) then
      !   allocate(lhs%evolved_leftright_midpoints_state, &
      !            mold=rhs%evolved_leftright_midpoints_state, stat=alloc_status)
      !   if(alloc_status /= 0) then
      !     error stop "Unable to allocate lhs%evolved_leftright_midpoints_state "// &
      !       "from rhs%evolved_leftright_midpoints_state in fvleg_t%assign"
      !   end if
      ! endif
      lhs%evolved_leftright_midpoints_state = rhs%evolved_leftright_midpoints_state

      ! if(.not. allocated(lhs%downup_midpoints_reference_state)) then
      !   allocate(lhs%downup_midpoints_reference_state, &
      !            mold=rhs%downup_midpoints_reference_state, stat=alloc_status)
      !   if(alloc_status /= 0) then
      !     error stop "Unable to allocate lhs%downup_midpoints_reference_state "// &
      !       "from rhs%downup_midpoints_reference_state in fvleg_t%assign"
      !   end if
      ! endif
      lhs%downup_midpoints_reference_state = rhs%downup_midpoints_reference_state

      ! if(.not. allocated(lhs%leftright_midpoints_reference_state)) then
      !   allocate(lhs%leftright_midpoints_reference_state, &
      !            mold=rhs%leftright_midpoints_reference_state, stat=alloc_status)
      !   if(alloc_status /= 0) then
      !     error stop "Unable to allocate lhs%leftright_midpoints_reference_state "// &
      !       "from rhs%leftright_midpoints_reference_state in fvleg_t%assign"
      !   end if
      ! endif
      lhs%leftright_midpoints_reference_state = rhs%leftright_midpoints_reference_state

      ! if(.not. allocated(lhs%reconstructed_state)) then
      !   allocate(lhs%reconstructed_state, mold=rhs%reconstructed_state, stat=alloc_status)
      !   if(alloc_status /= 0) then
      !     error stop "Unable to allocate lhs%reconstructed_state from "// &
      !       "rhs%reconstructed_state in fvleg_t%assign"
      !   end if
      ! endif
      lhs%reconstructed_state = rhs%reconstructed_state
      lhs%initiated = .true.

    class default
      error stop 'fvleg_t%assign_fvleg: unsupported class'
    end select
    call rhs%clean_temp(calling_function='assign_fvleg (rhs)', line=__LINE__)
  end subroutine assign_fvleg

end module mod_fvleg
