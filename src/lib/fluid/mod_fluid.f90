module mod_fluid
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use, intrinsic :: ieee_arithmetic

#ifdef USE_OPENMP
  use omp_lib
#endif /* USE_OPENMP */

  use mod_globals, only: debug_print, print_evolved_cell_data, print_recon_data
  use mod_nondimensionalization, only: rho_0, v_0, p_0, e_0
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_floating_point_utils, only: near_zero
  use mod_flux_array, only: get_fluxes, flux_array_t
  use mod_surrogate, only: surrogate
  use mod_units
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t
  use mod_surrogate, only: surrogate
  use mod_grid, only: grid_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use hdf5_interface, only: hdf5_file
  use mod_finite_volume_schemes, only: finite_volume_scheme_t
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_flux_tensor, only: operator(.dot.), H => flux_tensor_t
  use mod_time_integrator_factory, only: time_integrator_factory

  implicit none

  private
  public :: fluid_t, new_fluid

  type, extends(integrand_t) :: fluid_t
    real(rk), dimension(:, :), allocatable :: rho    !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: rho_u  !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: rho_v  !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: rho_E  !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: u      !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: v      !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: p      !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: cs     !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: mach   !< (i, j); Conserved quantities
    logical :: prim_vars_updated = .false.
  contains
    procedure, public :: initialize
    procedure, private :: initialize_from_ini
    procedure, private :: initialize_from_hdf5
    procedure, public :: t => time_derivative
    procedure, public :: residual_smoother
    procedure, public :: sanity_check
    procedure, public :: calculate_derived_quantities
    procedure, public :: apply_boundary_conditions
    procedure, nopass, private :: flux_edges
    procedure, pass(lhs), public :: type_plus_type => add_fluid
    procedure, pass(lhs), public :: type_minus_type => subtract_fluid
    procedure, pass(lhs), public :: type_mul_real => fluid_mul_real
    procedure, pass(rhs), public :: real_mul_type => real_mul_fluid
    procedure, pass(lhs), public :: assign => assign_fluid
    procedure, public :: force_finalization
    procedure, private, nopass :: add_fields
    final :: finalize
  end type fluid_t

contains

  function new_fluid(input, finite_volume_scheme) result(fluid)
    !< Fluid constructor
    class(input_t), intent(in) :: input
    class(finite_volume_scheme_t), intent(in) :: finite_volume_scheme
    type(fluid_t), pointer :: fluid

    allocate(fluid)
    call fluid%initialize(input, finite_volume_scheme)
  end function new_fluid

  subroutine initialize(self, input, finite_volume_scheme)
    class(fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(finite_volume_scheme_t), intent(in) :: finite_volume_scheme
    class(strategy), pointer :: time_integrator => null()

    integer(ik) :: alloc_status, i, j, ilo, ihi, jlo, jhi

    alloc_status = 0
    call debug_print('Initializing fluid_t', __FILE__, __LINE__)

    associate(imin=>finite_volume_scheme%grid%ilo_bc_cell, &
              imax=>finite_volume_scheme%grid%ihi_bc_cell, &
              jmin=>finite_volume_scheme%grid%jlo_bc_cell, &
              jmax=>finite_volume_scheme%grid%jhi_bc_cell)

      allocate(self%rho(imin:imax, jmin:jmax))
      allocate(self%u(imin:imax, jmin:jmax))
      allocate(self%v(imin:imax, jmin:jmax))
      allocate(self%p(imin:imax, jmin:jmax))
      allocate(self%rho_u(imin:imax, jmin:jmax))
      allocate(self%rho_v(imin:imax, jmin:jmax))
      allocate(self%rho_E(imin:imax, jmin:jmax))
      allocate(self%cs(imin:imax, jmin:jmax))
      allocate(self%mach(imin:imax, jmin:jmax))
    end associate

    time_integrator => time_integrator_factory(input)
    allocate(self%time_integrator, source=time_integrator, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%time_integrator"
    deallocate(time_integrator)

    if(input%read_init_cond_from_file .or. input%restart_from_file) then
      call self%initialize_from_hdf5(input, finite_volume_scheme)
    else
      call self%initialize_from_ini(input)
    end if

    call eos%sound_speed(p=self%p, rho=self%rho, cs=self%cs)

    ilo = lbound(self%rho, dim=1)
    ihi = ubound(self%rho, dim=1)
    jlo = lbound(self%rho, dim=2)
    jhi = ubound(self%rho, dim=2)

    !$omp parallel default(none), &
    !$omp private(i, j), &
    !$omp firstprivate(ilo, ihi, jlo, jhi), &
    !$omp shared(self)
    !$omp simd
    do j = jlo, jhi
      do i = ilo, ihi
        self%mach(i, j) = sqrt(self%u(i, j)**2 + self%v(i, j)**2) / self%cs(i, j)
      end do
    end do
    !$omp end simd
    !$omp end parallel

    self%prim_vars_updated = .true.

  end subroutine initialize

  subroutine initialize_from_hdf5(self, input, finite_volume_scheme)
    !< Initialize from an .hdf5 file. The conserved variables are already allocated appropriately from
    !< from the grid class, but this will just initialize them to the values found in the hdf5 file
    class(fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(finite_volume_scheme_t), intent(in) :: finite_volume_scheme
    type(hdf5_file) :: h5
    logical :: file_exists
    character(:), allocatable :: filename
    character(32) :: str_buff = ''

    real(rk), dimension(:, :), allocatable :: density
    real(rk), dimension(:, :), allocatable :: x_velocity
    real(rk), dimension(:, :), allocatable :: y_velocity
    real(rk), dimension(:, :), allocatable :: pressure

    call debug_print('Initializing fluid_t from hdf5', __FILE__, __LINE__)

    if(input%restart_from_file) then
      filename = trim(input%restart_file)
    else
      filename = trim(input%initial_condition_file)
    end if

    file_exists = .false.
    inquire(file=filename, exist=file_exists)

    if(.not. file_exists) then
      write(*, '(a)') 'Error in finite_volume_scheme_t%initialize_from_hdf5(); file not found: "'//filename//'"'
      error stop 'Error in finite_volume_scheme_t%initialize_from_hdf5(); file not found, exiting...'
    end if

    call h5%initialize(filename=filename, status='old', action='r')

    call h5%get('/density', density)
    if(input%restart_from_file) then
      call h5%readattr('/density', 'units', str_buff)
      select case(trim(str_buff))
      case('g/cc')
      case default
        error stop "Unknown density units in .h5 file. Acceptable units are 'g/cc'."
      end select
    end if

    call h5%get('/x_velocity', x_velocity)
    call h5%get('/y_velocity', y_velocity)
    if(input%restart_from_file) then
      call h5%readattr('/x_velocity', 'units', str_buff)
      select case(trim(str_buff))
      case('km/s')
        x_velocity = x_velocity * km_per_sec_to_cm_per_sec
        y_velocity = y_velocity * km_per_sec_to_cm_per_sec
      case('cm/s')
      case default
        error stop "Unknown velocity units in .h5 file. Acceptable units are 'km/s' and 'cm/s'."
      end select
    end if

    call h5%get('/pressure', pressure)
    if(input%restart_from_file) then
      call h5%readattr('/pressure', 'units', str_buff)
      select case(trim(str_buff))
      case('barye')
      case('Mbar')
        pressure = pressure * mega_bar_to_barye
      case default
        error stop "Unknown density units in .h5 file. Acceptable units are 'g/cc'."
      end select
    end if

    call h5%finalize()

    ! Non-dimensionalize
    density = density / rho_0
    x_velocity = x_velocity / v_0
    y_velocity = y_velocity / v_0
    pressure = pressure / p_0

    if(any(near_zero(pressure))) then
      error stop "Some (or all) of the pressure array is ~0 in fluid_t%initialize_from_hdf5"
    end if

    if(any(near_zero(density))) then
      error stop "Some (or all) of the density array is ~0 in fluid_t%initialize_from_hdf5"
    end if

    ! associate(imin=>finite_volume_scheme%grid%ilo_bc_cell, imax=>finite_volume_scheme%grid%ihi_bc_cell, &
    !           jmin=>finite_volume_scheme%grid%jlo_bc_cell, jmax=>finite_volume_scheme%grid%jhi_bc_cell)

    self%rho = density
    self%u = x_velocity
    self%v = y_velocity
    self%p = pressure
    ! end associate
    call eos%primitive_to_conserved(rho=self%rho, u=self%u, v=self%v, p=self%p, &
                                    rho_u=self%rho_u, rho_v=self%rho_v, rho_E=self%rho_E)

    write(*, '(a)') 'Initial fluid stats'
    write(*, '(a)') '==================================================='
    write(*, '(a, f0.4)') 'EOS Gamma:                     ', eos%get_gamma()
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max density    [non-dim]: ', minval(density), maxval(density)
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max x-velocity [non-dim]: ', minval(x_velocity), maxval(x_velocity)
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max y-velocity [non-dim]: ', minval(x_velocity), maxval(y_velocity)
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max pressure   [non-dim]: ', minval(pressure), maxval(pressure)
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max density    [dim]:     ', minval(density) * rho_0, maxval(density) * rho_0
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max x-velocity [dim]:     ', minval(x_velocity) * v_0, maxval(x_velocity) * v_0
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max y-velocity [dim]:     ', minval(y_velocity) * v_0, maxval(y_velocity) * v_0
    write(*, '(a, 2(es10.3, 1x))') 'Min/Max pressure   [dim]:     ', minval(pressure) * p_0, maxval(pressure) * p_0
    write(*, '(a)') '==================================================='
    write(*, *)

  end subroutine initialize_from_hdf5

  subroutine initialize_from_ini(self, input)
    !< Initialize from an .ini file. The conserved variables are already allocated appropriately from
    !< from the grid class, but this will just initialize them to the values found in the .ini file
    class(fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    real(rk) :: density, x_velocity, y_velocity, pressure

    ! Non-dimensionalize
    density = input%init_density / rho_0
    x_velocity = input%init_x_velocity / v_0
    y_velocity = input%init_y_velocity / v_0
    pressure = input%init_pressure / p_0

    call debug_print('Initializing fluid_t from .ini', __FILE__, __LINE__)
    write(*, '(a,4(f0.3, 1x))') 'Initializing with [rho,u,v,p]: ', &
      input%init_density, input%init_x_velocity, input%init_y_velocity, input%init_pressure

    self%rho = density
    self%u = x_velocity
    self%v = y_velocity
    self%p = pressure
    call eos%primitive_to_conserved(rho=self%rho, u=self%u, v=self%v, p=self%p, &
                                    rho_u=self%rho_u, rho_v=self%rho_v, rho_E=self%rho_E)

    if(near_zero(input%init_pressure)) then
      error stop "Some (or all) of the pressure array is ~0 in fluid_t%initialize_from_hdf5"
    end if

    if(near_zero(input%init_density)) then
      error stop "Some (or all) of the density array is ~0 in fluid_t%initialize_from_hdf5"
    end if

  end subroutine initialize_from_ini

  subroutine force_finalization(self)
    class(fluid_t), intent(inout) :: self

    call debug_print('Running fluid_t%force_finalization()', __FILE__, __LINE__)
    if(allocated(self%rho)) deallocate(self%rho)
    if(allocated(self%u)) deallocate(self%u)
    if(allocated(self%v)) deallocate(self%v)
    if(allocated(self%p)) deallocate(self%p)
    if(allocated(self%rho_u)) deallocate(self%rho_u)
    if(allocated(self%rho_v)) deallocate(self%rho_v)
    if(allocated(self%rho_E)) deallocate(self%rho_E)
    if(allocated(self%cs)) deallocate(self%cs)
    if(allocated(self%mach)) deallocate(self%mach)
    if(allocated(self%time_integrator)) deallocate(self%time_integrator)
  end subroutine force_finalization

  subroutine finalize(self)
    type(fluid_t), intent(inout) :: self

    call debug_print('Running fluid_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%rho)) deallocate(self%rho)
    if(allocated(self%u)) deallocate(self%u)
    if(allocated(self%v)) deallocate(self%v)
    if(allocated(self%p)) deallocate(self%p)
    if(allocated(self%rho_u)) deallocate(self%rho_u)
    if(allocated(self%rho_v)) deallocate(self%rho_v)
    if(allocated(self%rho_E)) deallocate(self%rho_E)
    if(allocated(self%cs)) deallocate(self%cs)
    if(allocated(self%mach)) deallocate(self%mach)
    if(allocated(self%time_integrator)) deallocate(self%time_integrator)
  end subroutine finalize

  function time_derivative(self, fv) result(d_dt)
    !< Implementation of dU/dt

    class(fluid_t), intent(in) :: self
    class(finite_volume_scheme_t), intent(inout) :: fv

    ! Locals
    class(integrand_t), allocatable :: d_dt !< dU/dt (integrand_t to satisfy parent interface)
    type(fluid_t), allocatable :: local_d_dt !< dU/dt
    integer(ik) :: alloc_status, error_code, i, j

    real(rk), dimension(:, :), allocatable :: evolved_corner_rho !< (i,j); Reconstructed rho at the corners
    real(rk), dimension(:, :), allocatable :: evolved_corner_u   !< (i,j); Reconstructed u at the corners
    real(rk), dimension(:, :), allocatable :: evolved_corner_v   !< (i,j); Reconstructed v at the corners
    real(rk), dimension(:, :), allocatable :: evolved_corner_p   !< (i,j); Reconstructed p at the corners
    real(rk), dimension(:, :), allocatable :: evolved_lr_mid_rho !< (i,j); Reconstructed rho at the left/right midpoints
    real(rk), dimension(:, :), allocatable :: evolved_lr_mid_u   !< (i,j); Reconstructed u at the left/right midpoints
    real(rk), dimension(:, :), allocatable :: evolved_lr_mid_v   !< (i,j); Reconstructed v at the left/right midpoints
    real(rk), dimension(:, :), allocatable :: evolved_lr_mid_p   !< (i,j); Reconstructed p at the left/right midpoints
    real(rk), dimension(:, :), allocatable :: evolved_du_mid_rho !< (i,j); Reconstructed rho at the down/up midpoints
    real(rk), dimension(:, :), allocatable :: evolved_du_mid_u   !< (i,j); Reconstructed u at the down/up midpoints
    real(rk), dimension(:, :), allocatable :: evolved_du_mid_v   !< (i,j); Reconstructed v at the down/up midpoints
    real(rk), dimension(:, :), allocatable :: evolved_du_mid_p   !< (i,j); Reconstructed p at the down/up midpoints

    real(rk), dimension(:, :, :), allocatable :: rho
    real(rk), dimension(:, :, :), allocatable :: u
    real(rk), dimension(:, :, :), allocatable :: v
    real(rk), dimension(:, :, :), allocatable :: p
    !< ((rho, u, v, p), i, j); Primitive variables at each cell center

    real(rk), dimension(:, :, :), allocatable :: evolved_corner_state
    !< ((rho, u, v, p), i, j); Reconstructed primitive variables at each corner

    ! Indexing the midpoints is a pain, so they're split by the up/down edges and left/right edges
    real(rk), dimension(:, :, :), allocatable :: evolved_downup_midpoints_state
    !< ((rho, u, v, p), i, j); Reconstructed primitive variables at each midpoint defined by vectors that go up/down
    !< (edges 2 (right edge) and 4 (left edge))

    real(rk), dimension(:, :, :), allocatable :: evolved_leftright_midpoints_state
    !< ((rho, u, v, p), i, j); Reconstructed primitive variables at each midpoint defined by vectors that go left/right
    !< (edges 1 (bottom edge) and 3 (top edge))

    real(rk), dimension(:, :, :), allocatable, target :: rho_recon_state
    !< ((corner1:midpoint4), i, j); reconstructed density, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint
    real(rk), dimension(:, :, :), allocatable, target :: u_recon_state
    !< ((corner1:midpoint4), i, j); reconstructed x-velocity, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint
    real(rk), dimension(:, :, :), allocatable, target :: v_recon_state
    !< ((corner1:midpoint4), i, j); reconstructed y-velocity, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint
    real(rk), dimension(:, :, :), allocatable, target :: p_recon_state
    !< ((corner1:midpoint4), i, j); reconstructed pressure, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint

    integer(ik), dimension(2) :: bounds

    call debug_print('Running fluid_t%time_derivative()', __FILE__, __LINE__)

    associate(imin=>fv%grid%ilo_bc_cell, imax=>fv%grid%ihi_bc_cell, &
              jmin=>fv%grid%jlo_bc_cell, jmax=>fv%grid%jhi_bc_cell)

      allocate(rho_recon_state(1:8, imin:imax, jmin:jmax))
      allocate(u_recon_state(1:8, imin:imax, jmin:jmax))
      allocate(v_recon_state(1:8, imin:imax, jmin:jmax))
      allocate(p_recon_state(1:8, imin:imax, jmin:jmax))
    end associate

    associate(imin_node=>fv%grid%ilo_node, imax_node=>fv%grid%ihi_node, &
              jmin_node=>fv%grid%jlo_node, jmax_node=>fv%grid%jhi_node, &
              imin_cell=>fv%grid%ilo_cell, imax_cell=>fv%grid%ihi_cell, &
              jmin_cell=>fv%grid%jlo_cell, jmax_cell=>fv%grid%jhi_cell)

      ! corners
      allocate(evolved_corner_state(4, imin_node:imax_node, jmin_node:jmax_node), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate evolved_corner_state in fluid_t%time_derivative()"

      ! left/right midpoints
      allocate(evolved_leftright_midpoints_state(4, imin_cell:imax_cell, jmin_node:jmax_node), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate evolved_leftright_midpoints_state in fluid_t%time_derivative()"

      ! down/up midpoints
      allocate(evolved_downup_midpoints_state(4, imin_node:imax_node, jmin_cell:jmax_cell), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate evolved_downup_midpoints_state"

      allocate(evolved_corner_rho(imin_node:imax_node, jmin_node:jmax_node))
      allocate(evolved_corner_u(imin_node:imax_node, jmin_node:jmax_node))
      allocate(evolved_corner_v(imin_node:imax_node, jmin_node:jmax_node))
      allocate(evolved_corner_p(imin_node:imax_node, jmin_node:jmax_node))

      allocate(evolved_lr_mid_rho(imin_cell:imax_cell, jmin_node:jmax_node))
      allocate(evolved_lr_mid_u(imin_cell:imax_cell, jmin_node:jmax_node))
      allocate(evolved_lr_mid_v(imin_cell:imax_cell, jmin_node:jmax_node))
      allocate(evolved_lr_mid_p(imin_cell:imax_cell, jmin_node:jmax_node))

      allocate(evolved_du_mid_rho(imin_node:imax_node, jmin_cell:jmax_cell))
      allocate(evolved_du_mid_u(imin_node:imax_node, jmin_cell:jmax_cell))
      allocate(evolved_du_mid_v(imin_node:imax_node, jmin_cell:jmax_cell))
      allocate(evolved_du_mid_p(imin_node:imax_node, jmin_cell:jmax_cell))

    end associate

    allocate(local_d_dt, source=self)

    if(.not. local_d_dt%prim_vars_updated) then
      call local_d_dt%calculate_derived_quantities()
    end if

    bounds = lbound(local_d_dt%rho)

    call local_d_dt%apply_boundary_conditions(fv)

    ! Now we can reconstruct the entire domain
    call fv%reconstruct(primitive_var=local_d_dt%rho, lbounds=bounds, &
                        reconstructed_var=rho_recon_state)

    call fv%reconstruct(primitive_var=local_d_dt%u, lbounds=bounds, &
                        reconstructed_var=u_recon_state)
    call fv%reconstruct(primitive_var=local_d_dt%v, lbounds=bounds, &
                        reconstructed_var=v_recon_state)
    call fv%reconstruct(primitive_var=local_d_dt%p, lbounds=bounds, &
                        reconstructed_var=p_recon_state)

    ! Apply the reconstructed state to the ghost layers
    call fv%apply_reconstructed_state_bc(recon_rho=rho_recon_state, recon_u=u_recon_state, &
                                         recon_v=v_recon_state, recon_p=p_recon_state, lbounds=bounds)

    call fv%evolution_operator%set_reconstructed_state_pointers(rho=rho_recon_state, &
                                                                u=u_recon_state, &
                                                                v=v_recon_state, &
                                                                p=p_recon_state, &
                                                                lbounds=bounds)
    ! error stop
    ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of down/up edge vectors
    call debug_print('Evolving down/up midpoints', __FILE__, __LINE__)
    bounds = lbound(evolved_du_mid_rho)
    call fv%evolution_operator%evolve(evolved_rho=evolved_du_mid_rho, evolved_u=evolved_du_mid_u, &
                                      evolved_v=evolved_du_mid_v, evolved_p=evolved_du_mid_p, &
                                      location='down/up midpoint', &
                                      lbounds=bounds, &
                                      error_code=error_code)

    ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of left/right edge vectors
    call debug_print('Evolving left/right midpoints', __FILE__, __LINE__)
    bounds = lbound(evolved_lr_mid_rho)
    call fv%evolution_operator%evolve(evolved_rho=evolved_lr_mid_rho, evolved_u=evolved_lr_mid_u, &
                                      evolved_v=evolved_lr_mid_v, evolved_p=evolved_lr_mid_p, &
                                      location='left/right midpoint', &
                                      lbounds=bounds, &
                                      error_code=error_code)

    ! Evolve, i.e. E0(R_omega), at all corner nodes
    call debug_print('Evolving corner nodes', __FILE__, __LINE__)
    bounds = lbound(evolved_corner_rho)
    call fv%evolution_operator%evolve(evolved_rho=evolved_corner_rho, evolved_u=evolved_corner_u, &
                                      evolved_v=evolved_corner_v, evolved_p=evolved_corner_p, &
                                      location='corner', &
                                      lbounds=bounds, &
                                      error_code=error_code)

    if(error_code /= 0) then
      fv%error_code = error_code
    end if

    nullify(fv%evolution_operator%reconstructed_rho)
    nullify(fv%evolution_operator%reconstructed_u)
    nullify(fv%evolution_operator%reconstructed_v)
    nullify(fv%evolution_operator%reconstructed_p)

    call self%flux_edges(grid=fv%grid, &
                         evolved_corner_rho=evolved_corner_rho, &
                         evolved_corner_u=evolved_corner_u, &
                         evolved_corner_v=evolved_corner_v, &
                         evolved_corner_p=evolved_corner_p, &
                         evolved_lr_mid_rho=evolved_lr_mid_rho, &
                         evolved_lr_mid_u=evolved_lr_mid_u, &
                         evolved_lr_mid_v=evolved_lr_mid_v, &
                         evolved_lr_mid_p=evolved_lr_mid_p, &
                         evolved_du_mid_rho=evolved_du_mid_rho, &
                         evolved_du_mid_u=evolved_du_mid_u, &
                         evolved_du_mid_v=evolved_du_mid_v, &
                         evolved_du_mid_p=evolved_du_mid_p, &
                         d_rho_dt=local_d_dt%rho, &
                         d_rhou_dt=local_d_dt%rho_u, &
                         d_rhov_dt=local_d_dt%rho_v, &
                         d_rhoE_dt=local_d_dt%rho_E)

    call move_alloc(local_d_dt, d_dt)

    call d_dt%set_temp(calling_function='fluid_t%time_derivative (d_dt)', line=__LINE__)

    ! Now deallocate everything
    deallocate(evolved_corner_rho)
    deallocate(evolved_corner_u)
    deallocate(evolved_corner_v)
    deallocate(evolved_corner_p)
    deallocate(evolved_lr_mid_rho)
    deallocate(evolved_lr_mid_u)
    deallocate(evolved_lr_mid_v)
    deallocate(evolved_lr_mid_p)
    deallocate(evolved_du_mid_rho)
    deallocate(evolved_du_mid_u)
    deallocate(evolved_du_mid_v)
    deallocate(evolved_du_mid_p)

  end function time_derivative

  subroutine apply_boundary_conditions(self, fv)
    !< Apply the primitive variable boundary conditions. This is separated out from the
    !< time derivative procedure, b/c it is helpful to view the bc in the contour files.
    !< This is also called by the integrand_t%integrate procedure after all of the
    !< time integration.

    class(fluid_t), intent(inout) :: self
    class(finite_volume_scheme_t), intent(inout) :: fv
    integer(ik), dimension(2) :: bounds

    bounds = lbound(self%rho)
    call fv%apply_primitive_vars_bc(rho=self%rho, &
                                    u=self%u, &
                                    v=self%v, &
                                    p=self%p, lbounds=bounds)

  end subroutine

  subroutine calculate_derived_quantities(self)
    !< Find derived quantities like sound speed, mach number, primitive variables

    class(fluid_t), intent(inout) :: self
    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi

    ilo = lbound(self%rho, dim=1)
    ihi = ubound(self%rho, dim=1)
    jlo = lbound(self%rho, dim=2)
    jhi = ubound(self%rho, dim=2)

    call eos%conserved_to_primitive(rho=self%rho, rho_u=self%rho_u, &
                                    rho_v=self%rho_v, rho_E=self%rho_E, &
                                    u=self%u, v=self%v, p=self%p)

    call eos%sound_speed(p=self%p, rho=self%rho, cs=self%cs)

    !$omp parallel default(none) &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp shared(self)
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        self%mach(i, j) = sqrt(self%u(i, j)**2 + self%v(i, j)**2) / self%cs(i, j)
      end do
    end do
    !$omp end do
    !$omp end parallel

    self%prim_vars_updated = .true.

  end subroutine

  subroutine flux_edges(grid, &
                        evolved_corner_rho, evolved_corner_u, evolved_corner_v, evolved_corner_p, &
                        evolved_lr_mid_rho, evolved_lr_mid_u, evolved_lr_mid_v, evolved_lr_mid_p, &
                        evolved_du_mid_rho, evolved_du_mid_u, evolved_du_mid_v, evolved_du_mid_p, &
                        d_rho_dt, d_rhou_dt, d_rhov_dt, d_rhoE_dt)
    !< Evaluate the fluxes along the edges. This is equation 13 in the paper
    class(grid_t), intent(in) :: grid

    real(rk), dimension(grid%ilo_node:, grid%jlo_node:), contiguous, intent(in) :: evolved_corner_rho !< (i,j); Reconstructed rho at the corners
    real(rk), dimension(grid%ilo_node:, grid%jlo_node:), contiguous, intent(in) :: evolved_corner_u   !< (i,j); Reconstructed u at the corners
    real(rk), dimension(grid%ilo_node:, grid%jlo_node:), contiguous, intent(in) :: evolved_corner_v   !< (i,j); Reconstructed v at the corners
    real(rk), dimension(grid%ilo_node:, grid%jlo_node:), contiguous, intent(in) :: evolved_corner_p   !< (i,j); Reconstructed p at the corners
    real(rk), dimension(grid%ilo_cell:, grid%jlo_node:), contiguous, intent(in) :: evolved_lr_mid_rho !< (i,j); Reconstructed rho at the left/right midpoints
    real(rk), dimension(grid%ilo_cell:, grid%jlo_node:), contiguous, intent(in) :: evolved_lr_mid_u   !< (i,j); Reconstructed u at the left/right midpoints
    real(rk), dimension(grid%ilo_cell:, grid%jlo_node:), contiguous, intent(in) :: evolved_lr_mid_v   !< (i,j); Reconstructed v at the left/right midpoints
    real(rk), dimension(grid%ilo_cell:, grid%jlo_node:), contiguous, intent(in) :: evolved_lr_mid_p   !< (i,j); Reconstructed p at the left/right midpoints
    real(rk), dimension(grid%ilo_node:, grid%jlo_cell:), contiguous, intent(in) :: evolved_du_mid_rho !< (i,j); Reconstructed rho at the down/up midpoints
    real(rk), dimension(grid%ilo_node:, grid%jlo_cell:), contiguous, intent(in) :: evolved_du_mid_u   !< (i,j); Reconstructed u at the down/up midpoints
    real(rk), dimension(grid%ilo_node:, grid%jlo_cell:), contiguous, intent(in) :: evolved_du_mid_v   !< (i,j); Reconstructed v at the down/up midpoints
    real(rk), dimension(grid%ilo_node:, grid%jlo_cell:), contiguous, intent(in) :: evolved_du_mid_p   !< (i,j); Reconstructed p at the down/up midpoints
    real(rk), dimension(grid%ilo_bc_cell:, grid%jlo_bc_cell:), contiguous, intent(inout) :: d_rho_dt
    real(rk), dimension(grid%ilo_bc_cell:, grid%jlo_bc_cell:), contiguous, intent(inout) :: d_rhou_dt
    real(rk), dimension(grid%ilo_bc_cell:, grid%jlo_bc_cell:), contiguous, intent(inout) :: d_rhov_dt
    real(rk), dimension(grid%ilo_bc_cell:, grid%jlo_bc_cell:), contiguous, intent(inout) :: d_rhoE_dt

    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: i, j, edge, xy
    real(rk), dimension(2, 4) :: n_hat
    real(rk), dimension(4) :: delta_l
    integer(ik), dimension(2) :: bounds

    type(flux_array_t) :: corner_fluxes
    type(flux_array_t) :: downup_mid_fluxes
    type(flux_array_t) :: leftright_mid_fluxes

    real(rk) :: bottom_rho_flux
    real(rk) :: bottom_rhou_flux
    real(rk) :: bottom_rhov_flux
    real(rk) :: bottom_rhoE_flux

    real(rk) :: right_rho_flux
    real(rk) :: right_rhou_flux
    real(rk) :: right_rhov_flux
    real(rk) :: right_rhoE_flux

    real(rk) :: top_rho_flux
    real(rk) :: top_rhou_flux
    real(rk) :: top_rhov_flux
    real(rk) :: top_rhoE_flux

    real(rk) :: left_rho_flux
    real(rk) :: left_rhou_flux
    real(rk) :: left_rhov_flux
    real(rk) :: left_rhoE_flux

    call debug_print('Running fluid_t%flux_edges()', __FILE__, __LINE__)

    ! Get the flux arrays for each corner node or midpoint
    ilo = grid%ilo_node; ihi = grid%ihi_node
    jlo = grid%jlo_node; jhi = grid%jhi_node
    bounds = [ilo, jlo]
    corner_fluxes = get_fluxes(rho=evolved_corner_rho, u=evolved_corner_u, v=evolved_corner_v, &
                               p=evolved_corner_p, lbounds=bounds)

    ilo = grid%ilo_cell; ihi = grid%ihi_cell
    jlo = grid%jlo_node; jhi = grid%jhi_node
    bounds = [ilo, jlo]
    leftright_mid_fluxes = get_fluxes(rho=evolved_lr_mid_rho, u=evolved_lr_mid_u, v=evolved_lr_mid_v, &
                                      p=evolved_lr_mid_p, lbounds=bounds)

    ilo = grid%ilo_node; ihi = grid%ihi_node
    jlo = grid%jlo_cell; jhi = grid%jhi_cell
    bounds = [ilo, jlo]
    downup_mid_fluxes = get_fluxes(rho=evolved_du_mid_rho, u=evolved_du_mid_u, v=evolved_du_mid_v, &
                                   p=evolved_du_mid_p, lbounds=bounds)

    ilo = grid%ilo_cell; ihi = grid%ihi_cell
    jlo = grid%jlo_cell; jhi = grid%jhi_cell
    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j, delta_l, n_hat) &
    !$omp private(bottom_rho_flux, bottom_rhou_flux, bottom_rhov_flux, bottom_rhoE_flux) &
    !$omp private(right_rho_flux, right_rhou_flux, right_rhov_flux, right_rhoE_flux) &
    !$omp private(left_rho_flux, left_rhou_flux, left_rhov_flux, left_rhoE_flux) &
    !$omp private(top_rho_flux, top_rhou_flux, top_rhov_flux, top_rhoE_flux) &
    !$omp shared(grid, corner_fluxes, leftright_mid_fluxes, downup_mid_fluxes) &
    !$omp shared(d_rho_dt, d_rhou_dt, d_rhov_dt, d_rhoE_dt)
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi

        delta_l = grid%cell_edge_lengths(:, i, j)

        do edge = 1, 4
          do xy = 1, 2
            n_hat(xy, edge) = grid%cell_edge_norm_vectors(xy, edge, i, j)
          end do
        end do

        associate(F_c1=>corner_fluxes%F(:, i, j), &
                  G_c1=>corner_fluxes%G(:, i, j), &
                  F_c2=>corner_fluxes%F(:, i + 1, j), &
                  G_c2=>corner_fluxes%G(:, i + 1, j), &
                  F_c3=>corner_fluxes%F(:, i + 1, j + 1), &
                  G_c3=>corner_fluxes%G(:, i + 1, j + 1), &
                  F_c4=>corner_fluxes%F(:, i, j + 1), &
                  G_c4=>corner_fluxes%G(:, i, j + 1), &
                  F_m1=>leftright_mid_fluxes%F(:, i, j), &
                  G_m1=>leftright_mid_fluxes%G(:, i, j), &
                  F_m2=>downup_mid_fluxes%F(:, i + 1, j), &
                  G_m2=>downup_mid_fluxes%G(:, i + 1, j), &
                  F_m3=>leftright_mid_fluxes%F(:, i, j + 1), &
                  G_m3=>leftright_mid_fluxes%G(:, i, j + 1), &
                  G_m4=>downup_mid_fluxes%G(:, i, j), &
                  F_m4=>downup_mid_fluxes%F(:, i, j), &
                  n_hat_1=>grid%cell_edge_norm_vectors(:, 1, i, j), &
                  n_hat_2=>grid%cell_edge_norm_vectors(:, 2, i, j), &
                  n_hat_3=>grid%cell_edge_norm_vectors(:, 3, i, j), &
                  n_hat_4=>grid%cell_edge_norm_vectors(:, 4, i, j))

          bottom_rho_flux = ((F_c1(1) + 4.0_rk * F_m1(1) + F_c2(1)) * n_hat_1(1) + &
                             (G_c1(1) + 4.0_rk * G_m1(1) + G_c2(1)) * n_hat_1(2)) * (delta_l(1) / 6.0_rk)
          bottom_rhou_flux = ((F_c1(2) + 4.0_rk * F_m1(2) + F_c2(2)) * n_hat_1(1) + &
                              (G_c1(2) + 4.0_rk * G_m1(2) + G_c2(2)) * n_hat_1(2)) * (delta_l(1) / 6.0_rk)
          bottom_rhov_flux = ((F_c1(3) + 4.0_rk * F_m1(3) + F_c2(3)) * n_hat_1(1) + &
                              (G_c1(3) + 4.0_rk * G_m1(3) + G_c2(3)) * n_hat_1(2)) * (delta_l(1) / 6.0_rk)
          bottom_rhoE_flux = ((F_c1(4) + 4.0_rk * F_m1(4) + F_c2(4)) * n_hat_1(1) + &
                              (G_c1(4) + 4.0_rk * G_m1(4) + G_c2(4)) * n_hat_1(2)) * (delta_l(1) / 6.0_rk)

          right_rho_flux = ((F_c2(1) + 4.0_rk * F_m2(1) + F_c3(1)) * n_hat_2(1) + &
                            (G_c2(1) + 4.0_rk * G_m2(1) + G_c3(1)) * n_hat_2(2)) * (delta_l(2) / 6.0_rk)
          right_rhou_flux = ((F_c2(2) + 4.0_rk * F_m2(2) + F_c3(2)) * n_hat_2(1) + &
                             (G_c2(2) + 4.0_rk * G_m2(2) + G_c3(2)) * n_hat_2(2)) * (delta_l(2) / 6.0_rk)
          right_rhov_flux = ((F_c2(3) + 4.0_rk * F_m2(3) + F_c3(3)) * n_hat_2(1) + &
                             (G_c2(3) + 4.0_rk * G_m2(3) + G_c3(3)) * n_hat_2(2)) * (delta_l(2) / 6.0_rk)
          right_rhoE_flux = ((F_c2(4) + 4.0_rk * F_m2(4) + F_c3(4)) * n_hat_2(1) + &
                             (G_c2(4) + 4.0_rk * G_m2(4) + G_c3(4)) * n_hat_2(2)) * (delta_l(2) / 6.0_rk)

          top_rho_flux = ((F_c3(1) + 4.0_rk * F_m3(1) + F_c4(1)) * n_hat_3(1) + &
                          (G_c3(1) + 4.0_rk * G_m3(1) + G_c4(1)) * n_hat_3(2)) * (delta_l(3) / 6.0_rk)
          top_rhou_flux = ((F_c3(2) + 4.0_rk * F_m3(2) + F_c4(2)) * n_hat_3(1) + &
                           (G_c3(2) + 4.0_rk * G_m3(2) + G_c4(2)) * n_hat_3(2)) * (delta_l(3) / 6.0_rk)
          top_rhov_flux = ((F_c3(3) + 4.0_rk * F_m3(3) + F_c4(3)) * n_hat_3(1) + &
                           (G_c3(3) + 4.0_rk * G_m3(3) + G_c4(3)) * n_hat_3(2)) * (delta_l(3) / 6.0_rk)
          top_rhoE_flux = ((F_c3(4) + 4.0_rk * F_m3(4) + F_c4(4)) * n_hat_3(1) + &
                           (G_c3(4) + 4.0_rk * G_m3(4) + G_c4(4)) * n_hat_3(2)) * (delta_l(3) / 6.0_rk)

          left_rho_flux = ((F_c4(1) + 4.0_rk * F_m4(1) + F_c1(1)) * n_hat_4(1) + &
                           (G_c4(1) + 4.0_rk * G_m4(1) + G_c1(1)) * n_hat_4(2)) * (delta_l(4) / 6.0_rk)
          left_rhou_flux = ((F_c4(2) + 4.0_rk * F_m4(2) + F_c1(2)) * n_hat_4(1) + &
                            (G_c4(2) + 4.0_rk * G_m4(2) + G_c1(2)) * n_hat_4(2)) * (delta_l(4) / 6.0_rk)
          left_rhov_flux = ((F_c4(3) + 4.0_rk * F_m4(3) + F_c1(3)) * n_hat_4(1) + &
                            (G_c4(3) + 4.0_rk * G_m4(3) + G_c1(3)) * n_hat_4(2)) * (delta_l(4) / 6.0_rk)
          left_rhoE_flux = ((F_c4(4) + 4.0_rk * F_m4(4) + F_c1(4)) * n_hat_4(1) + &
                            (G_c4(4) + 4.0_rk * G_m4(4) + G_c1(4)) * n_hat_4(2)) * (delta_l(4) / 6.0_rk)

          d_rho_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * &
                           (bottom_rho_flux + right_rho_flux + top_rho_flux + left_rho_flux)

          d_rhou_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * &
                            (bottom_rhou_flux + right_rhou_flux + top_rhou_flux + left_rhou_flux)

          d_rhov_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * &
                            (bottom_rhov_flux + right_rhov_flux + top_rhov_flux + left_rhov_flux)

          d_rhoE_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * &
                            (bottom_rhoE_flux + right_rhoE_flux + top_rhoE_flux + left_rhoE_flux)
        end associate

      end do ! i
    end do ! j
    !$omp end do
    !$omp end parallel

    ilo = grid%ilo_bc_cell; ihi = grid%ihi_bc_cell
    jlo = grid%jlo_bc_cell; jhi = grid%jhi_bc_cell
    d_rho_dt(ilo, :) = 0.0_rk
    d_rho_dt(ihi, :) = 0.0_rk
    d_rho_dt(:, jlo) = 0.0_rk
    d_rho_dt(:, jhi) = 0.0_rk

    d_rhou_dt(ilo, :) = 0.0_rk
    d_rhou_dt(ihi, :) = 0.0_rk
    d_rhou_dt(:, jlo) = 0.0_rk
    d_rhou_dt(:, jhi) = 0.0_rk

    d_rhov_dt(ilo, :) = 0.0_rk
    d_rhov_dt(ihi, :) = 0.0_rk
    d_rhov_dt(:, jlo) = 0.0_rk
    d_rhov_dt(:, jhi) = 0.0_rk

    d_rhoE_dt(ilo, :) = 0.0_rk
    d_rhoE_dt(ihi, :) = 0.0_rk
    d_rhoE_dt(:, jlo) = 0.0_rk
    d_rhoE_dt(:, jhi) = 0.0_rk

  end subroutine flux_edges

  subroutine residual_smoother(self)
    class(fluid_t), intent(inout) :: self

    ! integer(ik) :: ilo = 0
    integer(ik) :: ihi = 0
    integer(ik) :: jlo = 0
    integer(ik) :: jhi = 0
    ! integer(ik) :: i, j
    ! real(rk) :: eps
    ! real(rk) :: u, v
    ! real(rk), dimension(:, :), allocatable :: sound_speed

    ! call self%get_sound_speed(sound_speed)

    ! ilo = lbound(self%conserved_vars, dim=2)
    ! ihi = ubound(self%conserved_vars, dim=2)
    ! jlo = lbound(self%conserved_vars, dim=3)
    ! jhi = ubound(self%conserved_vars, dim=3)

    ! !$omp parallel default(none) private(i,j,u,v,ilo,ihi,jlo,jhi)
    ! !$omp do
    ! do j = jlo, jhi
    !   do i = ilo, ihi
    !     u = self%conserved_vars(2, i, j) / self%conserved_vars(1, i, j)
    !     v = self%conserved_vars(3, i, j) / self%conserved_vars(1, i, j)

    !     if(abs(u / sound_speed(i, j)) < 1e-8_rk) then
    !       self%conserved_vars(2, i, j) = 0.0_rk
    !     end if

    !     if(abs(v / sound_speed(i, j)) < 1e-8_rk) then
    !       self%conserved_vars(3, i, j) = 0.0_rk
    !     end if
    !   end do
    ! end do
    ! !$omp end do
    ! !$omp end parallel

    ! ! print*, minval((self%conserved_vars(2, :,:) / self%conserved_vars(1, :,:))/sound_speed)
    ! ! print*, maxval((self%conserved_vars(2, :,:) / self%conserved_vars(1, :,:))/sound_speed)
    ! ! print*, minval((self%conserved_vars(3, :,:) / self%conserved_vars(1, :,:))/sound_speed)
    ! ! print*, maxval((self%conserved_vars(3, :,:) / self%conserved_vars(1, :,:))/sound_speed)

    ! deallocate(sound_speed)

  end subroutine residual_smoother

  function subtract_fluid(lhs, rhs) result(difference)
    !< Implementation of the (-) operator for the fluid type
    class(fluid_t), intent(in) :: lhs
    class(integrand_t), intent(in) :: rhs

    class(integrand_t), allocatable :: difference
    type(fluid_t), allocatable :: local_difference

    call debug_print('Running fluid_t%subtract_fluid()', __FILE__, __LINE__)

    select type(rhs)
    class is(fluid_t)
      allocate(local_difference, source=lhs)
      call subtract_fields(a=lhs%rho, b=rhs%rho, c=local_difference%rho) ! c=a+b
      call subtract_fields(a=lhs%rho_u, b=rhs%rho_u, c=local_difference%rho_u) ! c=a+b
      call subtract_fields(a=lhs%rho_v, b=rhs%rho_v, c=local_difference%rho_v) ! c=a+b
      call subtract_fields(a=lhs%rho_E, b=rhs%rho_E, c=local_difference%rho_E) ! c=a+b
      local_difference%prim_vars_updated = .false.
    class default
      error stop 'fluid_t%subtract_fluid: unsupported rhs class'
    end select

    call move_alloc(local_difference, difference)
    call difference%set_temp(calling_function='subtract_fluid (difference)', line=__LINE__)
  end function subtract_fluid

  function add_fluid(lhs, rhs) result(sum)
    !< Implementation of the (+) operator for the fluid type
    class(fluid_t), intent(in) :: lhs
    class(integrand_t), intent(in) :: rhs

    class(integrand_t), allocatable :: sum
    type(fluid_t), allocatable :: local_sum

    call debug_print('Running fluid_t%add_fluid()', __FILE__, __LINE__)

    select type(rhs)
    class is(fluid_t)
      allocate(local_sum, source=lhs)
      call add_fields(a=lhs%rho, b=rhs%rho, c=local_sum%rho) ! c=a+b
      call add_fields(a=lhs%rho_u, b=rhs%rho_u, c=local_sum%rho_u) ! c=a+b
      call add_fields(a=lhs%rho_v, b=rhs%rho_v, c=local_sum%rho_v) ! c=a+b
      call add_fields(a=lhs%rho_E, b=rhs%rho_E, c=local_sum%rho_E) ! c=a+b
      local_sum%prim_vars_updated = .false.
    class default
      error stop 'fluid_t%add_fluid: unsupported rhs class'
    end select

    call local_sum%set_temp(calling_function='fluid_t%add_fluid(local_sum)', line=__LINE__)
    call move_alloc(local_sum, sum)
    call sum%set_temp(calling_function='fluid_t%add_fluid(sum)', line=__LINE__)
  end function add_fluid

  subroutine add_fields(a, b, c)
    !< Dumb routine for a vectorized version of c = a + b
    real(rk), dimension(:, :), contiguous, intent(in) :: a
    real(rk), dimension(:, :), contiguous, intent(in) :: b
    real(rk), dimension(:, :), contiguous, intent(inout) :: c
    integer(ik) :: i, j, k
    integer(ik) :: ilo = 0
    integer(ik) :: ihi = 0
    integer(ik) :: jlo = 0
    integer(ik) :: jhi = 0

    ilo = lbound(a, dim=1)
    ihi = ubound(a, dim=1)
    jlo = lbound(a, dim=2)
    jhi = ubound(a, dim=2)

    !$omp parallel default(none) &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp shared(a,b,c)
    !$omp do simd
    do j = jlo, jhi
      do i = ilo, ihi
        c(i, j) = a(i, j) + b(i, j)
      end do
    end do
    !$omp end do simd
    !$omp end parallel

  end subroutine add_fields

  subroutine subtract_fields(a, b, c)
    !< Dumb routine for a vectorized version of c = a - b
    real(rk), dimension(:, :), contiguous, intent(in) :: a
    real(rk), dimension(:, :), contiguous, intent(in) :: b
    real(rk), dimension(:, :), contiguous, intent(inout) :: c
    integer(ik) :: i, j, k
    integer(ik) :: ilo = 0
    integer(ik) :: ihi = 0
    integer(ik) :: jlo = 0
    integer(ik) :: jhi = 0

    ilo = lbound(a, dim=1)
    ihi = ubound(a, dim=1)
    jlo = lbound(a, dim=2)
    jhi = ubound(a, dim=2)

    !$omp parallel default(none) &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp shared(a,b,c)
    !$omp do simd
    do j = jlo, jhi
      do i = ilo, ihi
        c(i, j) = a(i, j) - b(i, j)
      end do
    end do
    !$omp end do simd
    !$omp end parallel

  end subroutine subtract_fields

  subroutine mult_fields(a, b, c)
    !< Dumb routine for a vectorized version of c = a * b
    real(rk), dimension(:, :), contiguous, intent(in) :: a
    real(rk), dimension(:, :), contiguous, intent(in) :: b
    real(rk), dimension(:, :), contiguous, intent(inout) :: c
    integer(ik) :: i, j, k
    integer(ik) :: ilo = 0
    integer(ik) :: ihi = 0
    integer(ik) :: jlo = 0
    integer(ik) :: jhi = 0

    ilo = lbound(a, dim=1)
    ihi = ubound(a, dim=1)
    jlo = lbound(a, dim=2)
    jhi = ubound(a, dim=2)

    !$omp parallel default(none) &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp shared(a,b,c)
    !$omp do simd
    do j = jlo, jhi
      do i = ilo, ihi
        c(i, j) = a(i, j) * b(i, j)
      end do
    end do
    !$omp end do simd
    !$omp end parallel

  end subroutine mult_fields

  subroutine mult_field_by_real(a, b, c)
    !< Dumb routine for a vectorized version of c = a * b
    real(rk), dimension(:, :), contiguous, intent(in) :: a
    real(rk), intent(in) :: b
    real(rk), dimension(:, :), contiguous, intent(inout) :: c
    integer(ik) :: i, j, k
    integer(ik) :: ilo = 0
    integer(ik) :: ihi = 0
    integer(ik) :: jlo = 0
    integer(ik) :: jhi = 0

    ilo = lbound(a, dim=1)
    ihi = ubound(a, dim=1)
    jlo = lbound(a, dim=2)
    jhi = ubound(a, dim=2)

    !$omp parallel default(none) &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp shared(a,b,c)
    !$omp do simd
    do j = jlo, jhi
      do i = ilo, ihi
        c(i, j) = a(i, j) * b
      end do
    end do
    !$omp end do simd
    !$omp end parallel
  end subroutine mult_field_by_real

  function fluid_mul_real(lhs, rhs) result(product)
    !< Implementation of the fluid * real operation
    class(fluid_t), intent(in) :: lhs
    real(rk), intent(in) :: rhs
    class(integrand_t), allocatable :: product

    class(fluid_t), allocatable :: local_product
    integer(ik) :: alloc_status

    call debug_print('Running fluid_t%fluid_mul_real()', __FILE__, __LINE__)

    allocate(local_product, source=lhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate local_product in fluid_t%fluid_mul_real"

    call mult_field_by_real(a=lhs%rho, b=rhs, c=local_product%rho)
    call mult_field_by_real(a=lhs%rho_u, b=rhs, c=local_product%rho_u)
    call mult_field_by_real(a=lhs%rho_v, b=rhs, c=local_product%rho_v)
    call mult_field_by_real(a=lhs%rho_E, b=rhs, c=local_product%rho_E)
    local_product%prim_vars_updated = .false.

    call move_alloc(local_product, product)
    call product%set_temp(calling_function='fluid_mul_real (product)', line=__LINE__)
  end function fluid_mul_real

  function real_mul_fluid(lhs, rhs) result(product)
    !< Implementation of the real * fluid operation
    class(fluid_t), intent(in) :: rhs
    real(rk), intent(in) :: lhs
    class(integrand_t), allocatable :: product

    type(fluid_t), allocatable :: local_product
    integer(ik) :: alloc_status

    call debug_print('Running fluid_t%real_mul_fluid()', __FILE__, __LINE__)

    allocate(local_product, source=rhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate local_product in fluid_t%real_mul_fluid"

    local_product%time_integrator = rhs%time_integrator
    call mult_field_by_real(a=rhs%rho, b=lhs, c=local_product%rho)
    call mult_field_by_real(a=rhs%rho_u, b=lhs, c=local_product%rho_u)
    call mult_field_by_real(a=rhs%rho_v, b=lhs, c=local_product%rho_v)
    call mult_field_by_real(a=rhs%rho_E, b=lhs, c=local_product%rho_E)
    local_product%prim_vars_updated = .false.

    call move_alloc(local_product, product)
    call product%set_temp(calling_function='real_mul_fluid (product)', line=__LINE__)
  end function real_mul_fluid

  subroutine assign_fluid(lhs, rhs)
    !< Implementation of the (=) operator for the fluid type. e.g. lhs = rhs
    class(fluid_t), intent(inout) :: lhs
    class(integrand_t), intent(in) :: rhs
    integer(ik) :: alloc_status
    integer(ik) :: i, j, k
    integer(ik) :: ilo = 0
    integer(ik) :: ihi = 0
    integer(ik) :: jlo = 0
    integer(ik) :: jhi = 0

    alloc_status = 0
    call debug_print('Running fluid_t%assign_fluid', __FILE__, __LINE__)

    call rhs%guard_temp(calling_function='assign_fluid (rhs)', line=__LINE__)
    select type(rhs)
    class is(fluid_t)
      lhs%time_integrator = rhs%time_integrator
      lhs%rho = rhs%rho
      lhs%rho_u = rhs%rho_u
      lhs%rho_v = rhs%rho_v
      lhs%rho_E = rhs%rho_E
    class default
      error stop 'Error in fluid_t%assign_fluid: unsupported class'
    end select

    call rhs%clean_temp(calling_function='assign_fluid (rhs)', line=__LINE__)
  end subroutine assign_fluid

  subroutine sanity_check(self, error_code)
    !< Run checks on the conserved variables. Density and pressure need to be > 0. No NaNs or Inifinite numbers either.
    class(fluid_t), intent(in) :: self
    integer(ik) :: i, j, ilo, jlo, ihi, jhi
    logical :: invalid_numbers, negative_numbers
    integer(ik), intent(out) :: error_code

    error_code = 0

    negative_numbers = .false.
    invalid_numbers = .false.

    ilo = lbound(self%rho, dim=1)
    ihi = ubound(self%rho, dim=1)
    jlo = lbound(self%rho, dim=2)
    jhi = ubound(self%rho, dim=2)

    !$omp parallel default(none) &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp shared(self)
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(self%rho(i, j) < 0.0_rk) then
          error stop "Error: Negative density found in fluid_t%sanity_check()"
        end if
      end do
    end do
    !$omp end do nowait

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(self%p(i, j) < 0.0_rk) then
          error stop "Error: Negative pressure found in fluid_t%sanity_check()"
        end if
      end do
    end do
    !$omp end do nowait

    ! NaN checks

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%rho(i, j))) then
          error stop "Error: NaN density found in fluid_t%sanity_check()"
        end if
      end do
    end do
    !$omp end do nowait
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%u(i, j))) then
          error stop "Error: NaN x-velocity found in fluid_t%sanity_check()"
        end if
      end do
    end do
    !$omp end do nowait
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%v(i, j))) then
          error stop "Error: NaN y-velocity found in fluid_t%sanity_check()"
        end if
      end do
    end do
    !$omp end do nowait
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%p(i, j))) then
          write(std_err, '(a, i0, ", ", i0, a)') "NaN pressure found at (", i, j, ")"
          error stop "Error: NaN pressure found in fluid_t%sanity_check()"
        end if
      end do
    end do
    !$omp end do nowait

    !$omp end parallel

  end subroutine sanity_check

end module mod_fluid
