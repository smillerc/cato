module mod_base_fluid
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
  public :: base_fluid_t

  logical, parameter :: filter_small_mach = .false.

  type, extends(integrand_t) :: base_fluid_t
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
    ! Common, inherited procedures
    procedure, public :: initialize
    procedure, public :: residual_smoother
    procedure, public :: write_residual_history
    procedure, public :: sanity_check
    procedure, public :: calculate_derived_quantities

    ! Operators
    procedure, pass(lhs), public :: type_plus_type => add_fluid
    procedure, pass(lhs), public :: type_minus_type => subtract_fluid
    procedure, pass(lhs), public :: type_mul_real => fluid_mul_real
    procedure, pass(rhs), public :: real_mul_type => real_mul_fluid
    procedure, pass(lhs), public :: assign => assign_fluid

    ! Private procedures
    procedure, private :: initialize_from_ini
    procedure, private :: initialize_from_hdf5
    procedure, private, nopass :: add_fields
    procedure, private, nopass :: subtract_fields
    procedure, private, nopass :: mult_fields

    ! Deferred procedures
    ! procedure(basic), deferred :: force_finalization
  end type base_fluid_t

  abstract interface
    subroutine basic(self)
      import :: base_fluid_t
      class(base_fluid_t), intent(inout) :: self
    end subroutine
  end interface

contains

  ! function new_fluid(input, finite_volume_scheme) result(fluid)
  !   !< Fluid constructor
  !   class(input_t), intent(in) :: input
  !   class(finite_volume_scheme_t), intent(in) :: finite_volume_scheme
  !   type(base_fluid_t), pointer :: fluid

  !   allocate(fluid)
  !   call fluid%initialize(input, finite_volume_scheme)
  ! end function new_fluid

  subroutine initialize(self, input, finite_volume_scheme)
    class(base_fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(finite_volume_scheme_t), intent(in) :: finite_volume_scheme
    class(strategy), pointer :: time_integrator => null()

    integer(ik) :: alloc_status, i, j, ilo, ihi, jlo, jhi, io

    alloc_status = 0
    call debug_print('Initializing base_fluid_t', __FILE__, __LINE__)

    if(.not. scale_factors_set) then
      error stop "Error in base_fluid_t%initialize(), global non-dimensional "// &
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
    class(base_fluid_t), intent(inout) :: self
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

    call debug_print('Initializing base_fluid_t from hdf5', __FILE__, __LINE__)

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
      error stop "Some (or all) of the pressure array is ~0 in base_fluid_t%initialize_from_hdf5"
    end if

    if(any(near_zero(density))) then
      error stop "Some (or all) of the density array is ~0 in base_fluid_t%initialize_from_hdf5"
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
    class(base_fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    real(rk) :: density, x_velocity, y_velocity, pressure

    ! Non-dimensionalize
    density = input%init_density / rho_0
    x_velocity = input%init_x_velocity / v_0
    y_velocity = input%init_y_velocity / v_0
    pressure = input%init_pressure / p_0

    call debug_print('Initializing base_fluid_t from .ini', __FILE__, __LINE__)
    write(*, '(a,4(f0.3, 1x))') 'Initializing with [rho,u,v,p]: ', &
      input%init_density, input%init_x_velocity, input%init_y_velocity, input%init_pressure

    self%rho = density
    self%u = x_velocity
    self%v = y_velocity
    self%p = pressure
    call eos%primitive_to_conserved(rho=self%rho, u=self%u, v=self%v, p=self%p, &
                                    rho_u=self%rho_u, rho_v=self%rho_v, rho_E=self%rho_E)

    if(near_zero(input%init_pressure)) then
      error stop "Some (or all) of the pressure array is ~0 in base_fluid_t%initialize_from_hdf5"
    end if

    if(near_zero(input%init_density)) then
      error stop "Some (or all) of the density array is ~0 in base_fluid_t%initialize_from_hdf5"
    end if

  end subroutine initialize_from_ini

  subroutine write_residual_history(self, fv)
    !< This writes out the change in residual to a file for convergence history monitoring. This
    !< should not be called on a single instance of base_fluid_t. It should be called something like the following:
    !<
    !< dU_dt = U%t(fv, stage=1)          ! 1st stage
    !< dU1_dt = U_1%t(fv, stage=2)       ! 2nd stage
    !< R = dU1_dt - dU_dt                ! Difference in the stages, eg residuals
    !< call R%write_residual_history(fv) ! Now write out the difference
    !<

    class(base_fluid_t), intent(inout) :: self
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

  end subroutine write_residual_history

  subroutine calculate_derived_quantities(self)
    !< Find derived quantities like sound speed, mach number, primitive variables

    class(base_fluid_t), intent(inout) :: self
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

  subroutine residual_smoother(self)
    class(base_fluid_t), intent(inout) :: self

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
    class(base_fluid_t), intent(in) :: lhs
    class(integrand_t), intent(in) :: rhs

    class(integrand_t), allocatable :: difference
    type(base_fluid_t), allocatable :: local_difference

    call debug_print('Running base_fluid_t%subtract_fluid()', __FILE__, __LINE__)

    select type(rhs)
    class is(base_fluid_t)
      allocate(local_difference, source=lhs)
      call subtract_fields(a=lhs%rho, b=rhs%rho, c=local_difference%rho) ! c=a+b
      call subtract_fields(a=lhs%rho_u, b=rhs%rho_u, c=local_difference%rho_u) ! c=a+b
      call subtract_fields(a=lhs%rho_v, b=rhs%rho_v, c=local_difference%rho_v) ! c=a+b
      call subtract_fields(a=lhs%rho_E, b=rhs%rho_E, c=local_difference%rho_E) ! c=a+b
      local_difference%prim_vars_updated = .false.
    class default
      error stop 'base_fluid_t%subtract_fluid: unsupported rhs class'
    end select

    call move_alloc(local_difference, difference)
    call difference%set_temp(calling_function='subtract_fluid (difference)', line=__LINE__)
  end function subtract_fluid

  function add_fluid(lhs, rhs) result(sum)
    !< Implementation of the (+) operator for the fluid type
    class(base_fluid_t), intent(in) :: lhs
    class(integrand_t), intent(in) :: rhs

    class(integrand_t), allocatable :: sum  ! FIXME: GFortran doesn't like polymorphic allocatable function results - big memory leak here... :(
    type(base_fluid_t), allocatable :: local_sum

    call debug_print('Running base_fluid_t%add_fluid()', __FILE__, __LINE__)

    select type(rhs)
    class is(base_fluid_t)
      allocate(local_sum, source=lhs)
      call add_fields(a=lhs%rho, b=rhs%rho, c=local_sum%rho) ! c=a+b
      call add_fields(a=lhs%rho_u, b=rhs%rho_u, c=local_sum%rho_u) ! c=a+b
      call add_fields(a=lhs%rho_v, b=rhs%rho_v, c=local_sum%rho_v) ! c=a+b
      call add_fields(a=lhs%rho_E, b=rhs%rho_E, c=local_sum%rho_E) ! c=a+b
      local_sum%prim_vars_updated = .false.
    class default
      error stop 'base_fluid_t%add_fluid: unsupported rhs class'
    end select

    call local_sum%set_temp(calling_function='base_fluid_t%add_fluid(local_sum)', line=__LINE__)
    call move_alloc(local_sum, sum)
    call sum%set_temp(calling_function='base_fluid_t%add_fluid(sum)', line=__LINE__)
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
    class(base_fluid_t), intent(in) :: lhs
    real(rk), intent(in) :: rhs
    class(integrand_t), allocatable :: product

    class(base_fluid_t), allocatable :: local_product
    integer(ik) :: alloc_status

    call debug_print('Running base_fluid_t%fluid_mul_real()', __FILE__, __LINE__)

    allocate(local_product, source=lhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate local_product in base_fluid_t%fluid_mul_real"

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
    real(rk), intent(in) :: lhs
    class(base_fluid_t), intent(in) :: rhs
    type(integrand_t), allocatable :: product

    type(base_fluid_t), allocatable :: local_product
    integer(ik) :: alloc_status

    call debug_print('Running base_fluid_t%real_mul_fluid()', __FILE__, __LINE__)

    allocate(local_product, source=rhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate local_product in base_fluid_t%real_mul_fluid"

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
    class(base_fluid_t), intent(inout) :: lhs
    class(integrand_t), intent(in) :: rhs

    call debug_print('Running base_fluid_t%assign_fluid', __FILE__, __LINE__)

    call rhs%guard_temp(calling_function='assign_fluid (rhs)', line=__LINE__)
    select type(rhs)
    class is(base_fluid_t)
      lhs%time_integrator = rhs%time_integrator
      lhs%residual_hist_header_written = rhs%residual_hist_header_written
      lhs%rho = rhs%rho
      lhs%rho_u = rhs%rho_u
      lhs%rho_v = rhs%rho_v
      lhs%rho_E = rhs%rho_E
    class default
      error stop 'Error in base_fluid_t%assign_fluid: unsupported class'
    end select

    call rhs%clean_temp(calling_function='assign_fluid (rhs)', line=__LINE__)
  end subroutine assign_fluid

  subroutine sanity_check(self, error_code)
    !< Run checks on the conserved variables. Density and pressure need to be > 0. No NaNs or Inifinite numbers either.
    class(base_fluid_t), intent(in) :: self
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
          error stop "Error: NaN density found in base_fluid_t%sanity_check()"
        end if
      end do
    end do
    !$omp end do nowait
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%u(i, j))) then
          write(std_err, '(a, i0, ", ", i0, a)') "NaN x-velocity found at (", i, j, ")"
          error stop "Error: NaN x-velocity found in base_fluid_t%sanity_check()"
        end if
      end do
    end do
    !$omp end do nowait
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%v(i, j))) then
          write(std_err, '(a, i0, ", ", i0, a)') "NaN y-velocity found at (", i, j, ")"
          error stop "Error: NaN y-velocity found in base_fluid_t%sanity_check()"
        end if
      end do
    end do
    !$omp end do nowait
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%p(i, j))) then
          write(std_err, '(a, i0, ", ", i0, a)') "NaN pressure found at (", i, j, ")"
          error stop "Error: NaN pressure found in base_fluid_t%sanity_check()"
        end if
      end do
    end do
    !$omp end do nowait

    !$omp end parallel

  end subroutine sanity_check

end module mod_base_fluid
