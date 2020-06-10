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

  logical, parameter :: filter_small_mach = .false.

  type, extends(integrand_t) :: fluid_t
    real(rk), dimension(:, :), allocatable :: rho    !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: rho_u  !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: rho_v  !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: rho_E  !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: u      !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: v      !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: p      !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: cs     !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: mach_u   !< (i, j); Conserved quantities
    real(rk), dimension(:, :), allocatable :: mach_v   !< (i, j); Conserved quantities
    logical :: prim_vars_updated = .false.
    logical :: smooth_residuals = .true.
    character(len=32) :: residual_hist_file = 'residual_hist.csv'
    logical :: residual_hist_header_written = .false.
  contains
    procedure, public :: initialize
    procedure, private :: initialize_from_ini
    procedure, private :: initialize_from_hdf5
    procedure, public :: t => time_derivative
    procedure, public :: residual_smoother
    procedure, public :: write_residual_history
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

    integer(ik) :: alloc_status, i, j, ilo, ihi, jlo, jhi, io

    alloc_status = 0
    call debug_print('Initializing fluid_t', __FILE__, __LINE__)

    if(.not. scale_factors_set) then
      error stop "Error in fluid_t%initialize(), global non-dimensional "// &
        "scale factors haven't been set yet. These need to be set before fluid initialization"
    end if

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
      allocate(self%mach_u(imin:imax, jmin:jmax))
      allocate(self%mach_v(imin:imax, jmin:jmax))
    end associate

    self%smooth_residuals = input%smooth_residuals

    time_integrator => time_integrator_factory(input)
    allocate(self%time_integrator, source=time_integrator, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%time_integrator"
    deallocate(time_integrator)

    open(newunit=io, file=trim(self%residual_hist_file), status='replace')
    write(io, '(a)') 'iteration,time,rho,rho_u,rho_v,rho_E'
    close(io)

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
        self%mach_u(i, j) = self%u(i, j) / self%cs(i, j)
        self%mach_v(i, j) = self%v(i, j) / self%cs(i, j)
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
    ! error stop
    ! end associate
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
    if(allocated(self%mach_u)) deallocate(self%mach_u)
    if(allocated(self%mach_v)) deallocate(self%mach_v)
    if(allocated(self%time_integrator)) deallocate(self%time_integrator)
  end subroutine finalize

  subroutine write_residual_history(self, fv)
    !< This writes out the change in residual to a file for convergence history monitoring. This
    !< should not be called on a single instance of fluid_t. It should be called something like the following:
    !<
    !< dU_dt = U%t(fv, stage=1)          ! 1st stage
    !< dU1_dt = U_1%t(fv, stage=2)       ! 2nd stage
    !< R = dU1_dt - dU_dt                ! Difference in the stages, eg residuals
    !< call R%write_residual_history(fv) ! Now write out the difference
    !<

    class(fluid_t), intent(inout) :: self
    class(finite_volume_scheme_t), intent(in) :: fv

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
    write(io, '(i0, ",", 5(es16.6, ","))') fv%iteration, fv%time * t_0, rho_diff, rho_u_diff, rho_v_diff, rho_E_diff
    close(io)

  end subroutine

  function time_derivative(self, fv, stage) result(d_dt)
    !< Implementation of dU/dt

    class(fluid_t), intent(in) :: self
    class(finite_volume_scheme_t), intent(inout) :: fv
    integer(ik), intent(in) :: stage !< which stage in the time integration scheme are we in, e.g. RK2 stage 1

    ! Locals
    class(integrand_t), allocatable :: d_dt !< dU/dt (integrand_t to satisfy parent interface)
    type(fluid_t), allocatable :: local_d_dt !< dU/dt
    integer(ik) :: error_code

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
    call debug_print('Reconstructing density', __FILE__, __LINE__)
    call fv%reconstruct(primitive_var=local_d_dt%rho, lbounds=bounds, &
                        reconstructed_var=rho_recon_state, name='rho', stage=stage)

    call debug_print('Reconstructing x-velocity', __FILE__, __LINE__)
    call fv%reconstruct(primitive_var=local_d_dt%u, lbounds=bounds, &
                        reconstructed_var=u_recon_state, name='u', stage=stage)

    call debug_print('Reconstructing y-velocity', __FILE__, __LINE__)
    call fv%reconstruct(primitive_var=local_d_dt%v, lbounds=bounds, &
                        reconstructed_var=v_recon_state, name='v', stage=stage)

    call debug_print('Reconstructing pressure', __FILE__, __LINE__)
    call fv%reconstruct(primitive_var=local_d_dt%p, lbounds=bounds, &
                        reconstructed_var=p_recon_state, name='p', stage=stage)

    ! Apply the reconstructed state to the ghost layers
    call fv%apply_reconstructed_state_bc(recon_rho=rho_recon_state, recon_u=u_recon_state, &
                                         recon_v=v_recon_state, recon_p=p_recon_state, lbounds=bounds)

    call fv%evolution_operator%set_reconstructed_state_pointers(rho=rho_recon_state, &
                                                                u=u_recon_state, &
                                                                v=v_recon_state, &
                                                                p=p_recon_state, &
                                                                lbounds=bounds)

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

    deallocate(rho_recon_state)
    deallocate(u_recon_state)
    deallocate(v_recon_state)
    deallocate(p_recon_state)

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

  end subroutine apply_boundary_conditions

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

    call debug_print('Running fluid_t%assign_fluid', __FILE__, __LINE__)

    call rhs%guard_temp(calling_function='assign_fluid (rhs)', line=__LINE__)
    select type(rhs)
    class is(fluid_t)
      lhs%time_integrator = rhs%time_integrator
      lhs%residual_hist_header_written = rhs%residual_hist_header_written
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

end module mod_fluid
