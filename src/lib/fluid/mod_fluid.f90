! MIT License
! Copyright (c) 2019 Sam Miller
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module mod_fluid
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use, intrinsic :: ieee_arithmetic

#ifdef USE_OPENMP
  use omp_lib
#endif /* USE_OPENMP */

  use mod_globals, only: debug_print, print_evolved_cell_data, print_recon_data, n_ghost_layers
  use mod_nondimensionalization, only: scale_factors_set, rho_0, v_0, p_0, e_0, t_0
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_floating_point_utils, only: near_zero, nearly_equal, neumaier_sum, neumaier_sum_2, neumaier_sum_3, neumaier_sum_4
  use mod_functional, only: operator(.sort.)
  use mod_units
  use mod_distinguisher, only: distinguish
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_grid, only: grid_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use hdf5_interface, only: hdf5_file
  use mod_flux_tensor, only: operator(.dot.), H => flux_tensor_t
  use mod_flux_solver, only: flux_solver_t
  use mod_ausm_plus_solver, only: ausm_plus_solver_t
  use mod_fvleg_solver, only: fvleg_solver_t
  use mod_m_ausmpw_plus_solver, only: m_ausmpw_plus_solver_t
  use mod_ausmpw_plus_solver, only: ausmpw_plus_solver_t
  use mod_slau_solver, only: slau_solver_t

  implicit none

  private
  public :: fluid_t, new_fluid

  logical, parameter :: filter_small_mach = .false.

  type :: fluid_t
    !< Fluid solver physics package

    private ! make all private by default

    real(rk), dimension(:, :), allocatable, public :: rho    !< (i, j); Conserved quantities
    !dir$ attributes align:__ALIGNBYTES__ :: rho
    real(rk), dimension(:, :), allocatable, public :: rho_u  !< (i, j); Conserved quantities
    !dir$ attributes align:__ALIGNBYTES__ :: rho_u
    real(rk), dimension(:, :), allocatable, public :: rho_v  !< (i, j); Conserved quantities
    !dir$ attributes align:__ALIGNBYTES__ :: rho_v
    real(rk), dimension(:, :), allocatable, public :: rho_E  !< (i, j); Conserved quantities
    !dir$ attributes align:__ALIGNBYTES__ :: rho_E
    real(rk), dimension(:, :), allocatable, public :: u      !< (i, j); Conserved quantities
    !dir$ attributes align:__ALIGNBYTES__ :: u
    real(rk), dimension(:, :), allocatable, public :: v      !< (i, j); Conserved quantities
    !dir$ attributes align:__ALIGNBYTES__ :: v
    real(rk), dimension(:, :), allocatable, public :: p      !< (i, j); Conserved quantities
    !dir$ attributes align:__ALIGNBYTES__ :: p
    real(rk), dimension(:, :), allocatable, public :: cs     !< (i, j); Conserved quantities
    !dir$ attributes align:__ALIGNBYTES__ :: cs
    real(rk), dimension(:, :), allocatable, public :: mach_u   !< (i, j); Conserved quantities
    !dir$ attributes align:__ALIGNBYTES__ :: mach_u
    real(rk), dimension(:, :), allocatable, public :: mach_v   !< (i, j); Conserved quantities
    !dir$ attributes align:__ALIGNBYTES__ :: mach_v

    integer(ik), dimension(:, :), allocatable, public :: continuous_sensor !< (i, j); Is the cell continuous, or (non)linear discontinuous?
    !dir$ attributes align:__ALIGNBYTES__ :: continuous_sensor

    class(flux_solver_t), allocatable :: solver !< solver scheme used to flux quantities at cell interfaces

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
    procedure :: calculate_derived_quantities
    procedure :: sanity_check
    procedure :: ssp_rk2
    procedure :: ssp_rk3
    procedure :: get_continuity_sensor
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
    procedure, public :: set_time
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
    class(flux_solver_t), pointer :: solver => null()

    integer(ik) :: alloc_status, i, j, ilo, ihi, jlo, jhi, io

    alloc_status = 0
    call debug_print('Calling fluid_t%initialize()', __FILE__, __LINE__)

    if(.not. scale_factors_set) then
      error stop "Error in fluid_t%initialize(), global non-dimensional "// &
        "scale factors haven't been set yet. These need to be set before fluid initialization"
    end if

    associate(imin => grid%ilo_bc_cell, &
              imax => grid%ihi_bc_cell, &
              jmin => grid%jlo_bc_cell, &
              jmax => grid%jhi_bc_cell)

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

    select case(trim(input%flux_solver))
    case('FVLEG')
      allocate(fvleg_solver_t :: solver)
    case('AUSM+-u', 'AUSM+-up', 'AUSM+-up_all_speed')
      allocate(ausm_plus_solver_t :: solver)
    case('M-AUSMPW+')
      allocate(m_ausmpw_plus_solver_t :: solver)
    case('AUSMPW+')
      allocate(ausmpw_plus_solver_t :: solver)
    case('SLAU', 'SLAU2', 'SD-SLAU', 'SD-SLAU2')
      allocate(slau_solver_t :: solver)
    case default
      write(std_err, '(a)') "Invalid flux solver in fluid_t%initializte(). It must be one of the following: "// &
        "['FVLEG', 'AUSM+-u','AUSM+-a','AUSM+-up','AUSM+-up_all_speed', 'AUSMPW+', 'M-AUSMPW+', 'SLAU', 'SLAU2', 'SD-SLAU', 'SD-SLAU2'], "// &
        "the input was: '"//trim(input%flux_solver)//"'"

      error stop "Invalid flux solver in fluid_t%initializte(). It must be one of the following: "// &
    "['FVLEG', 'AUSM+-u','AUSM+-a','AUSM+-up','AUSM+-up_all_speed', 'AUSMPW+', 'M-AUSMPW+', 'SLAU', 'SLAU2', 'SD-SLAU', 'SD-SLAU2']"
    end select

    call solver%initialize(input)
    allocate(self%solver, source=solver)
    deallocate(solver)

    open(newunit=io, file=trim(self%residual_hist_file), status='replace')
    write(io, '(a)') 'iteration,time,rho,rho_u,rho_v,rho_E'
    close(io)

    if(input%read_init_cond_from_file .or. input%restart_from_file) then
      call self%initialize_from_hdf5(input)
    else
      call self%initialize_from_ini(input)
    end if

    call self%get_continuity_sensor()
    call eos%sound_speed(p=self%p, rho=self%rho, cs=self%cs)

    ilo = lbound(self%rho, dim=1)
    ihi = ubound(self%rho, dim=1)
    jlo = lbound(self%rho, dim=2)
    jhi = ubound(self%rho, dim=2)

    !$omp parallel default(none), &
    !$omp private(i, j), &
    !$omp firstprivate(ilo, ihi, jlo, jhi), &
    !$omp shared(self)
    !$omp do
    do j = jlo, jhi
      !$omp simd
      do i = ilo, ihi
        self%mach_u(i, j) = self%u(i, j) / self%cs(i, j)
        self%mach_v(i, j) = self%v(i, j) / self%cs(i, j)
      end do
    end do
    !$omp end do
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
      write(*, '(a)') 'Error in fluid_t%initialize_from_hdf5(); file not found: "'//filename//'"'
      error stop 'Error in fluid_t%initialize_from_hdf5(); file not found, exiting...'
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
    if(allocated(self%solver)) deallocate(self%solver)
  end subroutine finalize

  subroutine set_time(self, time, dt, iteration)
    !< Set the time statistics
    class(fluid_t), intent(inout) :: self
    real(rk), intent(in) :: time          !< simulation time
    real(rk), intent(in) :: dt            !< time-step
    integer(ik), intent(in) :: iteration  !< iteration

    call debug_print('Running fluid_t%set_time()', __FILE__, __LINE__)
    self%time = time
    self%iteration = iteration
    self%dt = dt
  end subroutine set_time

  subroutine integrate(self, dt, grid)
    !< Integrate in time
    class(fluid_t), intent(inout) :: self
    real(rk), intent(in) :: dt !< time step
    class(grid_t), intent(in) :: grid !< grid class - the solver needs grid topology

    call debug_print('Running fluid_t%integerate()', __FILE__, __LINE__)
    self%time = self%time + dt
    self%dt = dt

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

    call debug_print('Running fluid_t%time_derivative()', __FILE__, __LINE__)

    lbounds = lbound(self%rho)
    allocate(d_dt, source=self)
    call self%solver%solve(dt=self%dt, &
                           grid=grid, lbounds=lbounds, &
                           rho=self%rho, u=self%u, v=self%v, p=self%p, &
                           d_rho_dt=d_dt%rho, &
                           d_rho_u_dt=d_dt%rho_u, &
                           d_rho_v_dt=d_dt%rho_v, &
                           d_rho_E_dt=d_dt%rho_E)

  end function time_derivative

  subroutine calculate_derived_quantities(self)
    !< Find derived quantities like sound speed, mach number, primitive variables

    class(fluid_t), intent(inout) :: self
    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi

    call debug_print('Running fluid_t%calculate_derived_quantities()', __FILE__, __LINE__)
    ilo = lbound(self%rho, dim=1)
    ihi = ubound(self%rho, dim=1)
    jlo = lbound(self%rho, dim=2)
    jhi = ubound(self%rho, dim=2)

    call eos%conserved_to_primitive(rho=self%rho, rho_u=self%rho_u, &
                                    rho_v=self%rho_v, rho_E=self%rho_E, &
                                    u=self%u, v=self%v, p=self%p)

    call self%get_continuity_sensor()
    call eos%sound_speed(p=self%p, rho=self%rho, cs=self%cs)

    !$omp parallel default(none) &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp shared(self)
    !$omp do
    do j = jlo, jhi
      !$omp simd
      do i = ilo, ihi
        self%mach_u(i, j) = self%u(i, j) / self%cs(i, j)
        self%mach_v(i, j) = self%v(i, j) / self%cs(i, j)
      end do
    end do
    !$omp end do
    !$omp end parallel

    self%prim_vars_updated = .true.

  end subroutine calculate_derived_quantities

  subroutine get_continuity_sensor(self)
    !< Run a distinguishing step that determines which regions in the domain are smooth (continuous),
    !< or discontinuous (linear or non-linear)
    class(fluid_t), intent(inout) :: self

    call distinguish(rho=self%rho, u=self%u, v=self%v, p=self%p, continuity_sensor=self%continuous_sensor)
  end subroutine get_continuity_sensor

  subroutine residual_smoother(self)
    class(fluid_t), intent(inout) :: self

    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: i, j
    real(rk), parameter :: EPS = 5e-14_rk

    call debug_print('Running fluid_t%residual_smoother()', __FILE__, __LINE__)
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
    !$omp do
    do j = jlo, jhi
#ifdef __SIMD_ALIGN_OMP__
      !$omp simd aligned(a, b, c:__ALIGNBYTES__)
#else
      !$omp simd
#endif
      do i = ilo, ihi
        c(i, j) = a(i, j) + b(i, j)
        diff = abs(c(i, j) - a(i, j))
        threshold = abs(a(i, j) + b(i, j)) * 1e-9_rk
        if(diff < threshold) c(i, j) = a(i, j)
      end do
    end do
    !$omp end do
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
    !$omp do
    do j = jlo, jhi
#ifdef __SIMD_ALIGN_OMP__
      !$omp simd aligned(a, b, c:__ALIGNBYTES__)
#else
      !$omp simd
#endif
      do i = ilo, ihi
        c(i, j) = a(i, j) - b(i, j)
        diff = abs(c(i, j) - a(i, j))
        threshold = abs(a(i, j) + b(i, j)) * 1e-9_rk
        if(diff < threshold) then
          c(i, j) = a(i, j)
        end if
      end do
    end do
    !$omp end do
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
    !$omp do
    do j = jlo, jhi
#ifdef __SIMD_ALIGN_OMP__
      !$omp simd aligned(a, b, c:__ALIGNBYTES__)
#else
      !$omp simd
#endif
      do i = ilo, ihi
        c(i, j) = a(i, j) * b(i, j)
      end do
    end do
    !$omp end do
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
    !$omp do
    do j = jlo, jhi
#ifdef __SIMD_ALIGN_OMP__
      !$omp simd aligned(a, c:__ALIGNBYTES__)
#else
      !$omp simd
#endif
      do i = ilo, ihi
        c(i, j) = a(i, j) * b
      end do
    end do
    !$omp end do
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
    integer(ik) :: error_code

    call debug_print('Running fluid_t%assign_fluid()', __FILE__, __LINE__)

    ! call rhs%guard_temp(calling_function='assign_fluid (rhs)', line=__LINE__)

    ! lhs%time_integrator = rhs%time_integrator
    lhs%residual_hist_header_written = rhs%residual_hist_header_written
    lhs%rho = rhs%rho
    lhs%rho_u = rhs%rho_u
    lhs%rho_v = rhs%rho_v
    lhs%rho_E = rhs%rho_E
    if(allocated(lhs%solver)) deallocate(lhs%solver)
    allocate(lhs%solver, source=rhs%solver)
    lhs%time_integration_scheme = rhs%time_integration_scheme
    lhs%time = rhs%time
    lhs%dt = rhs%dt
    lhs%iteration = rhs%iteration
    lhs%prim_vars_updated = rhs%prim_vars_updated
    lhs%smooth_residuals = rhs%smooth_residuals
    lhs%residual_hist_file = rhs%residual_hist_file
    lhs%residual_hist_header_written = rhs%residual_hist_header_written

    ! call lhs%sanity_check(error_code)
    call lhs%calculate_derived_quantities()
  end subroutine assign_fluid

  subroutine sanity_check(self, error_code)
    !< Run checks on the conserved variables. Density and pressure need to be > 0. No NaNs or Inifinite numbers either.
    class(fluid_t), intent(in) :: self
    integer(ik) :: i, j, ilo, jlo, ihi, jhi
    logical :: invalid_numbers, negative_numbers
    integer(ik), intent(out) :: error_code

    call debug_print('Running fluid_t%sanity_check()', __FILE__, __LINE__)

    error_code = 0

    negative_numbers = .false.
    invalid_numbers = .false.

    ilo = lbound(self%rho, dim=1) + n_ghost_layers
    ihi = ubound(self%rho, dim=1) - n_ghost_layers
    jlo = lbound(self%rho, dim=2) + n_ghost_layers
    jhi = ubound(self%rho, dim=2) - n_ghost_layers

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
    integer(ik) :: error_code
    real(rk) :: dt

    call debug_print('Running fluid_t%ssp_rk3()', __FILE__, __LINE__)

    dt = U%dt
    allocate(U_1, source=U)
    allocate(U_2, source=U)
    allocate(R, source=U)

    ! 1st stage
    U_1 = U + U%t(grid, stage=1) * dt
    call U_1%residual_smoother()

    ! 2nd stage
    U_2 = U * (3.0_rk / 4.0_rk) &
          + U_1 * (1.0_rk / 4.0_rk) &
          + U_1%t(grid, stage=2) * ((1.0_rk / 4.0_rk) * dt)
    call U_2%residual_smoother()

    ! Final stage
    U = U * (1.0_rk / 3.0_rk) &
        + U_2 * (2.0_rk / 3.0_rk) &
        + U_2%t(grid, stage=3) * ((2.0_rk / 3.0_rk) * dt)
    call U%residual_smoother()
    call U%sanity_check(error_code)
    ! Convergence history
    call write_residual_history(first_stage=U_1, last_stage=U)

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
    integer(ik) :: error_code

    call debug_print('Running fluid_t%rk2()', __FILE__, __LINE__)

    allocate(U_1, source=U)
    allocate(R, source=U)

    ! 1st stage
    associate(dt => U%dt)
      call debug_print(new_line('a')//'Running fluid_t%ssp_rk2_t() 1st stage'//new_line('a'), __FILE__, __LINE__)
      U_1 = U + U%t(grid, stage=1) * dt
      call U_1%residual_smoother()

      ! Final stage
      call debug_print(new_line('a')//'Running fluid_t%ssp_rk2_t() 2nd stage'//new_line('a'), __FILE__, __LINE__)
      U = U * 0.5_rk + U_1 * 0.5_rk + &
          U_1%t(grid, stage=2) * (0.5_rk * dt)
      call U%residual_smoother()
      call U%sanity_check(error_code)

      ! ! Convergence history
      call write_residual_history(first_stage=U_1, last_stage=U)
    end associate

    deallocate(R)
    deallocate(U_1)

  end subroutine ssp_rk2

  subroutine write_residual_history(first_stage, last_stage)
    !< This writes out the change in residual to a file for convergence history monitoring.

    class(fluid_t), intent(in) :: first_stage
    class(fluid_t), intent(in) :: last_stage

    real(rk) :: rho_diff   !< difference in the rho residual
    real(rk) :: rho_u_diff !< difference in the rhou residual
    real(rk) :: rho_v_diff !< difference in the rhov residual
    real(rk) :: rho_E_diff !< difference in the rhoE residual
    integer(ik) :: io
    integer(ik) :: i, j, ilo, jlo, ihi, jhi

    call debug_print('Running fluid_t%write_residual_history()', __FILE__, __LINE__)

    ilo = lbound(last_stage%rho, dim=1) + n_ghost_layers
    ihi = ubound(last_stage%rho, dim=1) - n_ghost_layers
    jlo = lbound(last_stage%rho, dim=2) + n_ghost_layers
    jhi = ubound(last_stage%rho, dim=2) - n_ghost_layers

    rho_diff = maxval(abs(last_stage%rho(ilo:ihi, jlo:jhi) - first_stage%rho(ilo:ihi, jlo:jhi)))
    rho_u_diff = maxval(abs(last_stage%rho_u(ilo:ihi, jlo:jhi) - first_stage%rho_u(ilo:ihi, jlo:jhi)))
    rho_v_diff = maxval(abs(last_stage%rho_v(ilo:ihi, jlo:jhi) - first_stage%rho_v(ilo:ihi, jlo:jhi)))
    rho_E_diff = maxval(abs(last_stage%rho_E(ilo:ihi, jlo:jhi) - first_stage%rho_E(ilo:ihi, jlo:jhi)))

    open(newunit=io, file=trim(first_stage%residual_hist_file), status='old', position="append")
    write(io, '(i0, ",", 5(es16.6, ","))') first_stage%iteration, first_stage%time * t_0, &
      rho_diff, rho_u_diff, rho_v_diff, rho_E_diff
    close(io)

  end subroutine write_residual_history
end module mod_fluid
