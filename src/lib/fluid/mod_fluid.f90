module mod_fluid
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use, intrinsic :: ieee_arithmetic

#ifdef USE_OPENMP
  use omp_lib
#endif /* USE_OPENMP */

  use mod_globals, only: debug_print, print_evolved_cell_data, print_recon_data
  use mod_nondimensionalization, only: scale_factors_set, rho_0, v_0, p_0, e_0, t_0
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_floating_point_utils, only: near_zero, nearly_equal, neumaier_sum, neumaier_sum_2, neumaier_sum_3, neumaier_sum_4
  use mod_functional, only: operator(.sort.)
  use mod_flux_array, only: get_fluxes, flux_array_t
  use mod_units
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_grid, only: grid_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use hdf5_interface, only: hdf5_file
  use mod_flux_tensor, only: operator(.dot.), H => flux_tensor_t
  use mod_riemann_solver, only: riemann_solver_t

  implicit none

  private
  public :: fluid_t, new_fluid

  logical, parameter :: filter_small_mach = .false.

  type :: fluid_t
    private ! make all private by default

    real(rk), dimension(:, :), allocatable, public :: rho    !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable, public :: rho_u  !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable, public :: rho_v  !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable, public :: rho_E  !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable, public :: u      !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable, public :: v      !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable, public :: p      !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable, public :: cs     !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable, public :: mach_u   !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable, public :: mach_v   !< (i, j); Conserved quantities

    class(riemann_solver_t), allocatable :: solver !< solver scheme used to flux quantities at cell interfaces

    ! Time variables
    character(len=10) :: time_integration_scheme = 'ssp_rk2'
    real(rk) :: time = 0.0_rk !< current simulation time
    real(rk) :: dt = 0.0_rk   !< time step
    integer(ik) :: iteration = 0 !< current iteration number

    logical, public :: prim_vars_updated = .false.
    logical :: smooth_residuals = .true.

    ! Residual history
    character(len=32) :: residual_hist_file = 'residual_hist.csv'
    logical :: residual_hist_header_written = .false.
  contains
    ! Private methods
    private
    procedure :: initialize_from_ini
    procedure :: initialize_from_hdf5
    procedure :: residual_smoother
    procedure :: write_residual_history
    procedure :: calculate_derived_quantities
    procedure :: sanity_check
    procedure :: ssp_rk2
    procedure :: ssp_rk3
    procedure, nopass :: add_fields
    procedure, nopass :: mult_fields
    procedure, nopass :: subtract_fields

    ! Operators
    procedure, pass(lhs), public :: add_fluid
    procedure, pass(lhs), public :: subtract_fluid
    procedure, pass(lhs), public :: fluid_mul_real
    procedure, pass(rhs), public :: real_mul_fluid
    procedure, pass(lhs), public :: assign_fluid

    ! Public methods
    procedure, public :: initialize
    procedure, public :: integrate
    procedure, public :: t => time_derivative
    procedure, public :: force_finalization

    ! Finalizer
    final :: finalize

    ! Map operators to corresponding procedures
    generic :: operator(+) => add_fluid
    generic :: operator(-) => subtract_fluid
    generic :: operator(*) => real_mul_fluid, fluid_mul_real
    generic :: assignment(=) => assign_fluid
  end type fluid_t

contains

  function new_fluid(input, grid) result(fluid)
    !< Fluid constructor
    class(input_t), intent(in) :: input
    class(grid_t), intent(in) :: grid
    type(fluid_t), pointer :: fluid

    allocate(fluid)
    call fluid%initialize(input, grid)
  end function new_fluid

  subroutine initialize(self, input, grid)
    class(fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in) :: grid

    integer(ik) :: alloc_status, i, j, ilo, ihi, jlo, jhi, io

    alloc_status = 0
    call debug_print('Initializing fluid_t', __FILE__, __LINE__)

    if(.not. scale_factors_set) then
      error stop "Error in fluid_t%initialize(), global non-dimensional "// &
        "scale factors haven't been set yet. These need to be set before fluid initialization"
    end if

    associate(imin=>grid%ilo_bc_cell, &
              imax=>grid%ihi_bc_cell, &
              jmin=>grid%jlo_bc_cell, &
              jmax=>grid%jhi_bc_cell)

      allocate(self%rho(imin:imax, jmin:jmax))
      allocate(self%u(imin:imax, jmin:jmax))
      allocate(self%v(imin:imax, jmin:jmax))
      allocate(self%p(imin:imax, jmin:jmax))
      allocate(self%rho_u(imin:imax, jmin:jmax))
      allocate(self%rho_v(imin:imax, jmin:jmax))
      allocate(self%rho_E(imin:imax, jmin:jmax))
      allocate(self%cs(imin:imax, jmin:jmax))
      allocate(self%mach_u(imin:imax, jmin:jmax))
      allocate(self%mach_v(imin:imax, jmin:jmax))
    end associate

    self%smooth_residuals = input%smooth_residuals

    self%time_integration_scheme = trim(input%time_integration_strategy)
    ! time_integrator => time_integrator_factory(input)
    ! allocate(self%time_integrator, source=time_integrator, stat=alloc_status)
    ! if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%time_integrator"
    ! deallocate(time_integrator)

    open(newunit=io, file=trim(self%residual_hist_file), status='replace')
    write(io, '(a)') 'iteration,time,rho,rho_u,rho_v,rho_E'
    close(io)

    if(input%read_init_cond_from_file .or. input%restart_from_file) then
      call self%initialize_from_hdf5(input)
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
        self%mach_u(i, j) = self%u(i, j) / self%cs(i, j)
        self%mach_v(i, j) = self%v(i, j) / self%cs(i, j)
      end do
    end do
    !$omp end simd
    !$omp end parallel

    self%prim_vars_updated = .true.

  end subroutine initialize

  subroutine initialize_from_hdf5(self, input)
    !< Initialize from an .hdf5 file. The conserved variables are already allocated appropriately from
    !< from the grid class, but this will just initialize them to the values found in the hdf5 file
    class(fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
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

    self%rho = density
    self%u = x_velocity
    self%v = y_velocity
    self%p = pressure

    call eos%primitive_to_conserved(rho=self%rho, u=self%u, v=self%v, p=self%p, &
                                    rho_u=self%rho_u, rho_v=self%rho_v, rho_E=self%rho_E)

    write(*, '(a)') 'Initial fluid stats'
    write(*, '(a)') '==================================================='
    write(*, '(a, f0.4)') 'EOS Gamma:                     ', eos%get_gamma()
    write(*, '(a, 2(es16.6, 1x))') 'Min/Max density    [non-dim]: ', minval(density), maxval(density)
    write(*, '(a, 2(es16.6, 1x))') 'Min/Max x-velocity [non-dim]: ', minval(x_velocity), maxval(x_velocity)
    write(*, '(a, 2(es16.6, 1x))') 'Min/Max y-velocity [non-dim]: ', minval(x_velocity), maxval(y_velocity)
    write(*, '(a, 2(es16.6, 1x))') 'Min/Max pressure   [non-dim]: ', minval(pressure), maxval(pressure)
    write(*, '(a, 2(es16.6, 1x))') 'Min/Max density    [dim]:     ', minval(density) * rho_0, maxval(density) * rho_0
    write(*, '(a, 2(es16.6, 1x))') 'Min/Max x-velocity [dim]:     ', minval(x_velocity) * v_0, maxval(x_velocity) * v_0
    write(*, '(a, 2(es16.6, 1x))') 'Min/Max y-velocity [dim]:     ', minval(y_velocity) * v_0, maxval(y_velocity) * v_0
    write(*, '(a, 2(es16.6, 1x))') 'Min/Max pressure   [dim]:     ', minval(pressure) * p_0, maxval(pressure) * p_0
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
    if(allocated(self%mach_u)) deallocate(self%mach_u)
    if(allocated(self%mach_v)) deallocate(self%mach_v)
    ! if(allocated(self%time_integrator)) deallocate(self%time_integrator)
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
    if(allocated(self%mach_u)) deallocate(self%mach_u)
    if(allocated(self%mach_v)) deallocate(self%mach_v)
    ! if(allocated(self%time_integrator)) deallocate(self%time_integrator)
  end subroutine finalize

  subroutine write_residual_history(self)
    !< This writes out the change in residual to a file for convergence history monitoring. This
    !< should not be called on a single instance of fluid_t. It should be called something like the following:
    !<
    !< dU_dt = U%t(fv, stage=1)          ! 1st stage
    !< dU1_dt = U_1%t(fv, stage=2)       ! 2nd stage
    !< R = dU1_dt - dU_dt                ! Difference in the stages, eg residuals
    !< call R%write_residual_history(fv) ! Now write out the difference
    !<

    class(fluid_t), intent(inout) :: self

    real(rk) :: rho_diff   !< difference in the rho residual
    real(rk) :: rho_u_diff !< difference in the rhou residual
    real(rk) :: rho_v_diff !< difference in the rhov residual
    real(rk) :: rho_E_diff !< difference in the rhoE residual
    integer(ik) :: io

    rho_diff = maxval(abs(self%rho))
    rho_u_diff = maxval(abs(self%rho_u))
    rho_v_diff = maxval(abs(self%rho_v))
    rho_E_diff = maxval(abs(self%rho_E))

    open(newunit=io, file=trim(self%residual_hist_file), status='old', position="append")
    write(io, '(i0, ",", 5(es16.6, ","))') self%iteration, self%time * t_0, rho_diff, rho_u_diff, rho_v_diff, rho_E_diff
    close(io)

  end subroutine

  subroutine integrate(self, dt, grid)
    !< Integrate in time
    class(fluid_t), intent(inout) :: self
    real(rk), intent(in) :: dt !< time step
    class(grid_t), intent(in) :: grid !< grid class - the solver needs grid topology
    type(fluid_t), allocatable :: d_dt !< dU/dt

    self%time = self%time + dt

    select case(trim(self%time_integration_scheme))
    case('ssp_rk2')
      call self%ssp_rk2(grid)
    case('ssp_rk3')
      call self%ssp_rk3(grid)
    case default
      error stop "Error: Unknown time integration scheme in fluid_t%integrate()"
    end select
  end subroutine integrate

  function time_derivative(self, grid, stage) result(d_dt)
    !< Implementation of the time derivative

    ! Inputs/Output
    class(fluid_t), intent(inout) :: self
    class(grid_t), intent(in) :: grid     !< grid class - the solver needs grid topology
    type(fluid_t), allocatable :: d_dt    !< dU/dt
    integer(ik), intent(in) :: stage      !< stage in the time integration scheme

    ! Locals
    integer(ik), dimension(2) :: lbounds

    lbounds = lbound(self%rho)

    call self%solver%solve(time=self%time, &
                           grid=grid, lbounds=lbounds, &
                           rho=self%rho, u=self%u, v=self%v, p=self%p, &
                           rho_u=self%rho_u, rho_v=self%rho_v, rho_E=self%rho_E)
  end function time_derivative

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
        self%mach_u(i, j) = self%u(i, j) / self%cs(i, j)
        self%mach_v(i, j) = self%v(i, j) / self%cs(i, j)
      end do
    end do
    !$omp end do
    !$omp end parallel

    self%prim_vars_updated = .true.

  end subroutine calculate_derived_quantities

  subroutine flux_edges(grid, &
                        evolved_corner_rho, evolved_corner_u, evolved_corner_v, evolved_corner_p, &
                        evolved_lr_mid_rho, evolved_lr_mid_u, evolved_lr_mid_v, evolved_lr_mid_p, &
                        evolved_du_mid_rho, evolved_du_mid_u, evolved_du_mid_v, evolved_du_mid_p, &
                        d_rho_dt, d_rhou_dt, d_rhov_dt, d_rhoE_dt)
    !< Evaluate the fluxes along the edges. This is equation 13 in the paper
    ! FIXME: Move this into a separate class/module - this will allow easier unit testing
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
    integer(ik) :: i, j, k, edge, xy
    real(rk), dimension(2, 4) :: n_hat
    real(rk), dimension(4) :: delta_l
    integer(ik), dimension(2) :: bounds

    type(flux_array_t) :: corner_fluxes
    type(flux_array_t) :: downup_mid_fluxes
    type(flux_array_t) :: leftright_mid_fluxes

    real(rk) :: f_sum, g_sum
    real(rk) :: diff

    real(rk), dimension(4) :: bottom_flux
    real(rk), dimension(4) :: right_flux
    real(rk), dimension(4) :: top_flux
    real(rk), dimension(4) :: left_flux

    real(rk) :: rho_flux, ave_rho_flux
    real(rk) :: rhou_flux, ave_rhou_flux
    real(rk) :: rhov_flux, ave_rhov_flux
    real(rk) :: rhoE_flux, ave_rhoE_flux
    real(rk) :: threshold

    real(rk), parameter :: FLUX_EPS = 1e-13_rk !epsilon(1.0_rk)
    real(rk), parameter :: REL_THRESHOLD = 1e-5_rk

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
    !$omp private(i, j, k, delta_l, n_hat) &
    !$omp private(bottom_flux, right_flux, left_flux, top_flux) &
    !$omp private(f_sum, g_sum, diff) &
    !$omp shared(grid, corner_fluxes, leftright_mid_fluxes, downup_mid_fluxes) &
    !$omp shared(d_rho_dt, d_rhou_dt, d_rhov_dt, d_rhoE_dt) &
    !$omp private(rho_flux, ave_rho_flux, rhou_flux, ave_rhou_flux, rhov_flux, ave_rhov_flux, rhoE_flux, ave_rhoE_flux, threshold)
    !$omp do simd
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

          !  Bottom
          ! bottom_flux (rho, rhou, rhov, rhoE)
          do k = 1, 4
            f_sum = sum([F_c1(k), 4.0_rk * F_m1(k), F_c2(k)]) * n_hat_1(1)
            g_sum = sum([G_c1(k), 4.0_rk * G_m1(k), G_c2(k)]) * n_hat_1(2)
            ! if(abs(f_sum + g_sum) < epsilon(1.0_rk)) then !FIXME: is this the right way to do it?
            !   bottom_flux(k) = 0.0_rk
            ! else
            bottom_flux(k) = (f_sum + g_sum) * (delta_l(1) / 6.0_rk)
            ! end if
          end do

          !  Right
          ! right_flux (rho, rhou, rhov, rhoE)
          do k = 1, 4
            f_sum = sum([F_c2(k), 4.0_rk * F_m2(k), F_c3(k)]) * n_hat_2(1)
            g_sum = sum([G_c2(k), 4.0_rk * G_m2(k), G_c3(k)]) * n_hat_2(2)
            ! if(abs(f_sum + g_sum) < epsilon(1.0_rk)) then
            !   right_flux(k) = 0.0_rk
            ! else
            right_flux(k) = (f_sum + g_sum) * (delta_l(2) / 6.0_rk)
            ! end if
          end do

          ! Top
          ! top_flux (rho, rhou, rhov, rhoE)
          do k = 1, 4
            f_sum = sum([F_c3(k), 4.0_rk * F_m3(k), F_c4(k)]) * n_hat_3(1)
            g_sum = sum([G_c3(k), 4.0_rk * G_m3(k), G_c4(k)]) * n_hat_3(2)
            ! if(abs(f_sum + g_sum) < epsilon(1.0_rk)) then
            !   top_flux(k) = 0.0_rk
            ! else
            top_flux(k) = (f_sum + g_sum) * (delta_l(3) / 6.0_rk)
            ! end if
          end do

          ! Left
          ! left_flux (rho, rhou, rhov, rhoE)
          do k = 1, 4
            f_sum = sum([F_c4(k), 4.0_rk * F_m4(k), F_c1(k)]) * n_hat_4(1)
            g_sum = sum([G_c4(k), 4.0_rk * G_m4(k), G_c1(k)]) * n_hat_4(2)
            ! if(abs(f_sum + g_sum) < epsilon(1.0_rk)) then
            ! left_flux(k) = 0.0_rk
            ! else
            left_flux(k) = (f_sum + g_sum) * (delta_l(4) / 6.0_rk)
            ! end if
          end do

          ave_rho_flux = sum([abs(left_flux(1)), abs(right_flux(1)), abs(top_flux(1)), abs(bottom_flux(1))]) / 4.0_rk
          ave_rhou_flux = sum([abs(left_flux(2)), abs(right_flux(2)), abs(top_flux(2)), abs(bottom_flux(2))]) / 4.0_rk
          ave_rhov_flux = sum([abs(left_flux(3)), abs(right_flux(3)), abs(top_flux(3)), abs(bottom_flux(3))]) / 4.0_rk
          ave_rhoE_flux = sum([abs(left_flux(4)), abs(right_flux(4)), abs(top_flux(4)), abs(bottom_flux(4))]) / 4.0_rk

          ! Now sum them all up together
          rho_flux = neumaier_sum_4([left_flux(1), right_flux(1), top_flux(1), bottom_flux(1)])
          rhou_flux = neumaier_sum_4([left_flux(2), right_flux(2), top_flux(2), bottom_flux(2)])
          rhov_flux = neumaier_sum_4([left_flux(3), right_flux(3), top_flux(3), bottom_flux(3)])
          rhoE_flux = neumaier_sum_4([left_flux(4), right_flux(4), top_flux(4), bottom_flux(4)])

          threshold = abs(ave_rho_flux) * REL_THRESHOLD
          if(abs(rho_flux) < threshold .or. abs(rho_flux) < epsilon(1.0_rk)) then
            rho_flux = 0.0_rk
          end if

          threshold = abs(ave_rhou_flux) * REL_THRESHOLD
          if(abs(rhou_flux) < threshold .or. abs(rhou_flux) < epsilon(1.0_rk)) then
            rhou_flux = 0.0_rk
          end if

          threshold = abs(ave_rhov_flux) * REL_THRESHOLD
          if(abs(rhov_flux) < threshold .or. abs(rhov_flux) < epsilon(1.0_rk)) then
            rhov_flux = 0.0_rk
          end if

          threshold = abs(ave_rhoE_flux) * REL_THRESHOLD
          if(abs(rhoE_flux) < threshold .or. abs(rhoE_flux) < epsilon(1.0_rk)) then
            rhoE_flux = 0.0_rk
          end if

          d_rho_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * rho_flux
          d_rhou_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * rhou_flux
          d_rhov_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * rhov_flux
          d_rhoE_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * rhoE_flux

          ! if(abs(rhov_flux) > 0.0_rk) then
          !   print *, 'i,j', i, j
          !   write(*, '(a, es16.6)') 'rhov_flux:                        : ', rhov_flux
          !   ! write(*, '(a, es16.6)') 'FLUX_EPS                          : ', FLUX_EPS
          !   ! write(*, '(a, es16.6)') 'ave_rhov_flux                     : ', ave_rhov_flux
          !   ! write(*, '(a, es16.6)') 'REL_THRESHOLD                     : ', REL_THRESHOLD
          !   ! write(*, '(a, es16.6)') 'abs(ave_rhov_flux) * REL_THRESHOLD: ', abs(ave_rhov_flux) * REL_THRESHOLD
          !   ! write(*, '(a, 4(es16.6, 1x))') ' -> ', left_flux(3), right_flux(3), top_flux(3), bottom_flux(3)

          !   print*, 'bottom flux'
          !   print*, 'n_hat_1: ', n_hat_1
          !   write(*, '(a, 4(es16.6))') "F_c1:          ", F_c1
          !   write(*, '(a, 4(es16.6))') "4.0_rk * F_m1: ", 4.0_rk * F_m1
          !   write(*, '(a, 4(es16.6))') "F_c2:          ", F_c2
          !   write(*, '(a, 4(es16.6))') "G_c1:          ", G_c1
          !   write(*, '(a, 4(es16.6))') "4.0_rk * G_m1: ", 4.0_rk * G_m1
          !   write(*, '(a, 4(es16.6))') "G_c2:          ", G_c2
          !   print*, '(delta_l(1) / 6.0_rk) ', (delta_l(1) / 6.0_rk)
          !   do k = 1, 4
          !     f_sum = sum([F_c1(k), 4.0_rk * F_m1(k), F_c2(k)]) * n_hat_1(1)
          !     g_sum = sum([G_c1(k), 4.0_rk * G_m1(k), G_c2(k)]) * n_hat_1(2)
          !     print*, 'f_sum', k, f_sum
          !     print*, 'g_sum', k, g_sum
          !   end do

          !   error stop
          ! end if

        end associate

      end do ! i
    end do ! j
    !$omp end do simd
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

    ! print *, 'fluxed vars differences'
    ! write(*, '(a, 6(es16.6))') "d/dt rho   : ", maxval(d_rho_dt(204, jlo + 1:jhi - 1)) - minval(d_rho_dt(204, jlo + 1:jhi - 1))
    ! write(*, '(a, 6(es16.6))') "d/dt rho u : ", maxval(d_rhou_dt(204, jlo + 1:jhi - 1)) - minval(d_rhou_dt(204, jlo + 1:jhi - 1))
    ! write(*, '(a, 6(es16.6))') "d/dt rho v : ", maxval(d_rhov_dt(204, jlo + 1:jhi - 1)) - minval(d_rhov_dt(204, jlo + 1:jhi - 1))
    ! write(*, '(a, 6(es16.6))') "d/dt rho E : ", maxval(d_rhoE_dt(204, jlo + 1:jhi - 1)) - minval(d_rhoE_dt(204, jlo + 1:jhi - 1))
    ! print *

  end subroutine flux_edges

  subroutine residual_smoother(self)
    class(fluid_t), intent(inout) :: self

    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: i, j
    real(rk), parameter :: EPS = 5e-14_rk

    if(self%smooth_residuals) then
      ilo = lbound(self%rho, dim=1) + 1
      ihi = ubound(self%rho, dim=1) - 1
      jlo = lbound(self%rho, dim=2) + 1
      jhi = ubound(self%rho, dim=2) - 1

      !$omp parallel default(none) &
      !$omp shared(self) &
      !$omp firstprivate(ilo,ihi,jlo,jhi) &
      !$omp private(i,j)
      !$omp do
      do j = jlo, jhi
        do i = ilo, ihi
          if(abs(self%rho(i, j) - self%rho(i - 1, j)) < EPS) self%rho(i, j) = self%rho(i - 1, j)
          if(abs(self%rho_u(i, j) - self%rho_u(i - 1, j)) < EPS) self%rho_u(i, j) = self%rho_u(i - 1, j)
          if(abs(self%rho_v(i, j) - self%rho_v(i - 1, j)) < EPS) self%rho_v(i, j) = self%rho_v(i - 1, j)
          if(abs(self%rho_E(i, j) - self%rho_E(i - 1, j)) < EPS) self%rho_E(i, j) = self%rho_E(i - 1, j)
        end do
      end do
      !$omp end do
      !$omp do
      do j = jlo, jhi
        do i = ilo, ihi
          if(abs(self%rho(i, j) - self%rho(i, j - 1)) < EPS) self%rho(i, j) = self%rho(i, j - 1)
          if(abs(self%rho_u(i, j) - self%rho_u(i, j - 1)) < EPS) self%rho_u(i, j) = self%rho_u(i, j - 1)
          if(abs(self%rho_v(i, j) - self%rho_v(i, j - 1)) < EPS) self%rho_v(i, j) = self%rho_v(i, j - 1)
          if(abs(self%rho_E(i, j) - self%rho_E(i, j - 1)) < EPS) self%rho_E(i, j) = self%rho_E(i, j - 1)
        end do
      end do
      !$omp end do
      !$omp end parallel
    end if

  end subroutine residual_smoother

  function subtract_fluid(lhs, rhs) result(difference)
    !< Implementation of the (-) operator for the fluid type
    class(fluid_t), intent(in) :: lhs
    class(fluid_t), intent(in) :: rhs
    type(fluid_t), allocatable :: difference

    call debug_print('Running fluid_t%subtract_fluid()', __FILE__, __LINE__)

    allocate(difference, source=lhs)
    call subtract_fields(a=lhs%rho, b=rhs%rho, c=difference%rho) ! c=a+b
    call subtract_fields(a=lhs%rho_u, b=rhs%rho_u, c=difference%rho_u) ! c=a+b
    call subtract_fields(a=lhs%rho_v, b=rhs%rho_v, c=difference%rho_v) ! c=a+b
    call subtract_fields(a=lhs%rho_E, b=rhs%rho_E, c=difference%rho_E) ! c=a+b
    difference%prim_vars_updated = .false.

    ! call difference%set_temp(calling_function='subtract_fluid (difference)', line=__LINE__)
  end function subtract_fluid

  function add_fluid(lhs, rhs) result(sum)
    !< Implementation of the (+) operator for the fluid type
    class(fluid_t), intent(in) :: lhs
    class(fluid_t), intent(in) :: rhs
    type(fluid_t), allocatable :: sum

    call debug_print('Running fluid_t%add_fluid()', __FILE__, __LINE__)

    allocate(sum, source=lhs)
    call add_fields(a=lhs%rho, b=rhs%rho, c=sum%rho) ! c=a+b
    call add_fields(a=lhs%rho_u, b=rhs%rho_u, c=sum%rho_u) ! c=a+b
    call add_fields(a=lhs%rho_v, b=rhs%rho_v, c=sum%rho_v) ! c=a+b
    call add_fields(a=lhs%rho_E, b=rhs%rho_E, c=sum%rho_E) ! c=a+b
    sum%prim_vars_updated = .false.

    ! call sum%set_temp(calling_function='fluid_t%add_fluid(sum)', line=__LINE__)
    ! call move_alloc(sum, sum)
    ! call sum%set_temp(calling_function='fluid_t%add_fluid(sum)', line=__LINE__)
  end function add_fluid

  subroutine add_fields(a, b, c)
    !< Dumb routine for a vectorized version of c = a + b
    real(rk), dimension(:, :), contiguous, intent(in) :: a
    real(rk), dimension(:, :), contiguous, intent(in) :: b
    real(rk), dimension(:, :), contiguous, intent(inout) :: c
    real(rk) :: diff, threshold
    integer(ik) :: i, j
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
    !$omp private(i, j, diff, threshold) &
    !$omp shared(a,b,c)
    !$omp do simd
    do j = jlo, jhi
      do i = ilo, ihi
        c(i, j) = a(i, j) + b(i, j)
        diff = abs(c(i, j) - a(i, j))
        threshold = abs(a(i, j) + b(i, j)) * 1e-9_rk
        if(diff < threshold) then
          ! if (diff > 0.0_rk) write(*,'(8(es16.6))') diff, threshold, c(i,j), a(i,j), b(i,j)
          c(i, j) = a(i, j)
        end if
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
    real(rk) :: diff, threshold
    integer(ik) :: i, j
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
    !$omp private(i, j, diff, threshold) &
    !$omp shared(a,b,c)
    !$omp do simd
    do j = jlo, jhi
      do i = ilo, ihi
        c(i, j) = a(i, j) - b(i, j)
        diff = abs(c(i, j) - a(i, j))
        threshold = abs(a(i, j) + b(i, j)) * 1e-9_rk
        if(diff < threshold) then
          c(i, j) = a(i, j)
        end if
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
    integer(ik) :: i, j
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
    integer(ik) :: i, j
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
    type(fluid_t), allocatable :: product

    integer(ik) :: alloc_status

    call debug_print('Running fluid_t%fluid_mul_real()', __FILE__, __LINE__)

    allocate(product, source=lhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate product in fluid_t%fluid_mul_real"

    call mult_field_by_real(a=lhs%rho, b=rhs, c=product%rho)
    call mult_field_by_real(a=lhs%rho_u, b=rhs, c=product%rho_u)
    call mult_field_by_real(a=lhs%rho_v, b=rhs, c=product%rho_v)
    call mult_field_by_real(a=lhs%rho_E, b=rhs, c=product%rho_E)
    product%prim_vars_updated = .false.

    ! call move_alloc(product, product)
    ! call product%set_temp(calling_function='fluid_mul_real (product)', line=__LINE__)
  end function fluid_mul_real

  function real_mul_fluid(lhs, rhs) result(product)
    !< Implementation of the real * fluid operation
    class(fluid_t), intent(in) :: rhs
    real(rk), intent(in) :: lhs
    type(fluid_t), allocatable :: product
    integer(ik) :: alloc_status

    call debug_print('Running fluid_t%real_mul_fluid()', __FILE__, __LINE__)

    allocate(product, source=rhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate product in fluid_t%real_mul_fluid"

    ! product%time_integrator = rhs%time_integrator
    call mult_field_by_real(a=rhs%rho, b=lhs, c=product%rho)
    call mult_field_by_real(a=rhs%rho_u, b=lhs, c=product%rho_u)
    call mult_field_by_real(a=rhs%rho_v, b=lhs, c=product%rho_v)
    call mult_field_by_real(a=rhs%rho_E, b=lhs, c=product%rho_E)
    product%prim_vars_updated = .false.

    ! call move_alloc(product, product)
    ! call product%set_temp(calling_function='real_mul_fluid (product)', line=__LINE__)
  end function real_mul_fluid

  subroutine assign_fluid(lhs, rhs)
    !< Implementation of the (=) operator for the fluid type. e.g. lhs = rhs
    class(fluid_t), intent(inout) :: lhs
    type(fluid_t), intent(in) :: rhs

    call debug_print('Running fluid_t%assign_fluid', __FILE__, __LINE__)

    ! call rhs%guard_temp(calling_function='assign_fluid (rhs)', line=__LINE__)

    ! lhs%time_integrator = rhs%time_integrator
    lhs%residual_hist_header_written = rhs%residual_hist_header_written
    lhs%rho = rhs%rho
    lhs%rho_u = rhs%rho_u
    lhs%rho_v = rhs%rho_v
    lhs%rho_E = rhs%rho_E

    ! call rhs%clean_temp(calling_function='assign_fluid (rhs)', line=__LINE__)
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
          write(std_err, '(a, i0, ", ", i0, a)') "Negative density found at (", i, j, ")"
          error stop "Negative density!"
        end if
      end do
    end do
    !$omp end do nowait

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(self%p(i, j) < 0.0_rk) then
          write(std_err, '(a, i0, ", ", i0, a)') "Negative pressure found at (", i, j, ")"
          error stop "Negative pressure!"
        end if
      end do
    end do
    !$omp end do nowait

    ! NaN checks

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%rho(i, j))) then
          write(std_err, '(a, i0, ", ", i0, a)') "NaN density found at (", i, j, ")"
          error stop "Error: NaN density found in fluid_t%sanity_check()"
        end if
      end do
    end do
    !$omp end do nowait
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%u(i, j))) then
          write(std_err, '(a, i0, ", ", i0, a)') "NaN x-velocity found at (", i, j, ")"
          error stop "Error: NaN x-velocity found in fluid_t%sanity_check()"
        end if
      end do
    end do
    !$omp end do nowait
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%v(i, j))) then
          write(std_err, '(a, i0, ", ", i0, a)') "NaN y-velocity found at (", i, j, ")"
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

  subroutine ssp_rk3(U, grid)
    !< Strong-stability preserving Runge-Kutta 3rd order
    class(fluid_t), intent(inout) :: U
    class(grid_t), intent(in) :: grid

    type(fluid_t), allocatable :: U_1 !< first stage
    type(fluid_t), allocatable :: U_2 !< second stage
    type(fluid_t), allocatable :: R !< hist

    allocate(U_1, source=U)
    allocate(U_2, source=U)
    allocate(R, source=U)

    associate(dt=>U%dt)
      ! 1st stage
      U_1 = U + dt * U%t(grid, stage=1)
      call U_1%residual_smoother()

      ! 2nd stage
      U_2 = (3.0_rk / 4.0_rk) * U &
            + (1.0_rk / 4.0_rk) * U_1 &
            + (1.0_rk / 4.0_rk) * dt * U_1%t(grid, stage=2)
      call U_2%residual_smoother()

      ! Final stage
      U = (1.0_rk / 3.0_rk) * U &
          + (2.0_rk / 3.0_rk) * U_2 &
          + (2.0_rk / 3.0_rk) * dt * U_2%t(grid, stage=3)
      call U%residual_smoother()

      ! Convergence history
      R = U - U_1
      call R%write_residual_history()
    end associate

    deallocate(R)
    deallocate(U_1)
    deallocate(U_2)

  end subroutine ssp_rk3

  subroutine ssp_rk2(U, grid)
    !< Strong-stability preserving Runge-Kutta 2nd order
    class(fluid_t), intent(inout) :: U
    class(grid_t), intent(in) :: grid

    type(fluid_t), allocatable :: U_1 !< first stage
    type(fluid_t), allocatable :: R !< hist

    allocate(U_1, source=U)
    allocate(R, source=U)

    ! 1st stage
    associate(dt=>U%dt)
      call debug_print('Running ssp_rk2_t 1st stage', __FILE__, __LINE__)
      U_1 = U + U%t(grid, stage=1) * dt
      call U_1%residual_smoother()

      ! Final stage
      call debug_print('Running ssp_rk2_t 2nd stage', __FILE__, __LINE__)
      U = 0.5_rk * U + 0.5_rk * U_1 + &
          (0.5_rk * dt) * U_1%t(grid, stage=2)
      call U%residual_smoother()

      ! Convergence history
      R = U - U_1
      call R%write_residual_history()
    end associate

    deallocate(R)
    deallocate(U_1)

  end subroutine ssp_rk2

end module mod_fluid
