! MIT License
! Copyright (c) 2020 Sam Miller
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

module mod_global_fluid
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use, intrinsic :: ieee_arithmetic

  use mod_error, only: ALL_OK, NEG_DENSITY, NEG_PRESSURE, NANS_FOUND, error_msg
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
  use mod_local_fluid, only: local_fluid_t

  implicit none

  private
  public :: global_fluid_t

  logical, parameter :: filter_small_mach = .false.

  type :: global_fluid_t
    !< Fluid solver physics package

    private ! make all private by default

    real(rk), dimension(:, :), codimension[:], allocatable, public :: rho    !< (i, j); density (conserved & primitive)
    real(rk), dimension(:, :), codimension[:], allocatable, public :: rho_u  !< (i, j); density * x-velocity (conserved)
    real(rk), dimension(:, :), codimension[:], allocatable, public :: rho_v  !< (i, j); density * y-velocity (conserved)
    real(rk), dimension(:, :), codimension[:], allocatable, public :: rho_E  !< (i, j); density * total energy (conserved)
    real(rk), dimension(:, :), codimension[:], allocatable, public :: u      !< (i, j); x-velocity (primitive)
    real(rk), dimension(:, :), codimension[:], allocatable, public :: v      !< (i, j); y-velocity (primitive)
    real(rk), dimension(:, :), codimension[:], allocatable, public :: p      !< (i, j); pressure (primitive)
    real(rk), dimension(:, :), codimension[:], allocatable, public :: cs     !< (i, j); sound speed
    real(rk), dimension(:, :), codimension[:], allocatable, public :: mach_u !< (i, j); mach number in the x-direction
    real(rk), dimension(:, :), codimension[:], allocatable, public :: mach_v !< (i, j); mach number in the y-direction

    ! Intel compiler alignment hints
    !dir$ attributes align:__ALIGNBYTES__ :: rho
    !dir$ attributes align:__ALIGNBYTES__ :: rho_u
    !dir$ attributes align:__ALIGNBYTES__ :: rho_v
    !dir$ attributes align:__ALIGNBYTES__ :: rho_E
    !dir$ attributes align:__ALIGNBYTES__ :: u
    !dir$ attributes align:__ALIGNBYTES__ :: v
    !dir$ attributes align:__ALIGNBYTES__ :: p
    !dir$ attributes align:__ALIGNBYTES__ :: cs
    !dir$ attributes align:__ALIGNBYTES__ :: mach_u
    !dir$ attributes align:__ALIGNBYTES__ :: mach_v

    class(flux_solver_t), allocatable :: solver !< solver scheme used to flux quantities at cell interfaces

    ! Coarray image indexing
    integer(ik), dimension(4) :: neighbors !< neighboring images
    integer(ik), dimension(4) :: node_dims_global = 0 !< (ilo, ihi, jlo, jhi); global node-based dimensions (w/ ghost)
    integer(ik), dimension(4) :: node_dims_img = 0    !< (ilo, ihi, jlo, jhi); local image node-based dimensions (w/ ghost)
    integer(ik), dimension(4) :: cell_dims_global = 0 !< (ilo, ihi, jlo, jhi); global node-based dimensions (w/ ghost)
    integer(ik), dimension(4) :: cell_dims_img = 0    !< (ilo, ihi, jlo, jhi); local image node-based dimensions (w/ ghost)

    logical :: img_on_ilo_bc = .false. !< is this image on the boundary for ilo?
    logical :: img_on_ihi_bc = .false. !< is this image on the boundary for ihi?
    logical :: img_on_jlo_bc = .false. !< is this image on the boundary for jlo?
    logical :: img_on_jhi_bc = .false. !< is this image on the boundary for jhi?

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
    procedure :: synchronize_edges
    procedure :: calculate_derived_quantities
    procedure :: sanity_check
    procedure :: ssp_rk2
    procedure :: ssp_rk3
    procedure :: gather
    procedure :: get_image_partitioning

    ! Operators
    procedure, pass(lhs), public :: add_local_fluid
    procedure, pass(lhs), public :: subtract_local_fluid
    procedure, pass(lhs), public :: local_fluid_mul_real
    procedure, pass(rhs), public :: real_mul_local_fluid
    procedure, pass(lhs), public :: assign_local_fluid

    ! Public methods
    procedure, public :: initialize
    procedure, public :: set_time
    procedure, public :: integrate
    procedure, public :: t => time_derivative
    procedure, public :: force_finalization

    ! Finalizer
    final :: finalize

    ! Map operators to corresponding procedures
    generic :: operator(+) => add_local_fluid
    generic :: operator(-) => subtract_local_fluid
    generic :: operator(*) => real_mul_local_fluid, local_fluid_mul_real
    generic :: assignment(=) => assign_local_fluid
  end type global_fluid_t

contains

  ! function new_fluid(input, grid) result(fluid)
  !   !< Fluid constructor
  !   class(input_t), intent(in) :: input
  !   class(grid_t), intent(in) :: grid
  !   type(global_fluid_t), pointer :: fluid

  !   allocate(fluid)
  !   call fluid%initialize(input, grid)
  ! end function new_fluid

  subroutine initialize(self, input, grid)
    class(global_fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in) :: grid
    class(flux_solver_t), pointer :: solver => null()

    integer(ik) :: alloc_status, i, j, ilo, ihi, jlo, jhi, io

    alloc_status = 0
    call debug_print('Calling global_fluid_t%initialize()', __FILE__, __LINE__)

    if(.not. scale_factors_set) then
      call error_msg(module='mod_fluid', class='global_fluid_t', procedure='initialize', &
            message="Global non-dimensional scale factors haven't been set yet. These need to be set before fluid initialization", &
                     file_name=__FILE__, line_number=__LINE__)
    end if

    associate(imin => grid%ilo_bc_cell, &
              imax => grid%ihi_bc_cell, &
              jmin => grid%jlo_bc_cell, &
              jmax => grid%jhi_bc_cell)

      allocate(self%rho(imin:imax, jmin:jmax)[*])
      allocate(self%u(imin:imax, jmin:jmax)[*])
      allocate(self%v(imin:imax, jmin:jmax)[*])
      allocate(self%p(imin:imax, jmin:jmax)[*])
      allocate(self%rho_u(imin:imax, jmin:jmax)[*])
      allocate(self%rho_v(imin:imax, jmin:jmax)[*])
      allocate(self%rho_E(imin:imax, jmin:jmax)[*])
      allocate(self%cs(imin:imax, jmin:jmax)[*])
      allocate(self%mach_u(imin:imax, jmin:jmax)[*])
      allocate(self%mach_v(imin:imax, jmin:jmax)[*])
    end associate

    self%smooth_residuals = input%smooth_residuals

    self%time_integration_scheme = trim(input%time_integration_strategy)

    select case(trim(input%flux_solver))
    case('FVLEG')
      allocate(fvleg_solver_t :: solver)
    case('AUSM+-u', 'AUSM+-up', 'AUSM+-up_all_speed')
      error stop "There are issues in the AUSM+ solver for now; exiting..."
      allocate(ausm_plus_solver_t :: solver)
    case('M-AUSMPW+')
      allocate(m_ausmpw_plus_solver_t :: solver)
    case('AUSMPW+')
      allocate(ausmpw_plus_solver_t :: solver)
    case('SLAU', 'SLAU2', 'SD-SLAU', 'SD-SLAU2')
      allocate(slau_solver_t :: solver)
    case default
      call error_msg(module='mod_fluid', class='global_fluid_t', procedure='initialize', &
                     message="Invalid flux solver. It must be one of the following: "// &
                     "['FVLEG', 'AUSM+-u','AUSM+-a','AUSM+-up','AUSM+-up_all_speed', "// &
                     "'AUSMPW+', 'M-AUSMPW+', 'SLAU', 'SLAU2', 'SD-SLAU', 'SD-SLAU2'], "// &
                     "the input was: '"//trim(input%flux_solver)//"'", &
                     file_name=__FILE__, line_number=__LINE__)
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
    class(global_fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    type(hdf5_file) :: h5
    logical :: file_exists
    character(:), allocatable :: filename
    character(32) :: str_buff = ''

    real(rk), dimension(:, :), allocatable :: density
    real(rk), dimension(:, :), allocatable :: x_velocity
    real(rk), dimension(:, :), allocatable :: y_velocity
    real(rk), dimension(:, :), allocatable :: pressure

    call debug_print('Initializing global_fluid_t from hdf5', __FILE__, __LINE__)

    if(input%restart_from_file) then
      filename = trim(input%restart_file)
    else
      filename = trim(input%initial_condition_file)
    end if

    file_exists = .false.
    inquire(file=filename, exist=file_exists)

    if(.not. file_exists) then
      call error_msg(module='mod_fluid', class='global_fluid_t', procedure='initialize_from_hdf5', &
                     message='File not found: "'//filename//'"', file_name=__FILE__, line_number=__LINE__)
    end if

    call h5%initialize(filename=filename, status='old', action='r')

    call h5%get('/density', density)
    if(input%restart_from_file) then
      call h5%readattr('/density', 'units', str_buff)
      select case(trim(str_buff))
      case('g/cc')
      case default
        call error_msg(module='mod_fluid', class='global_fluid_t', procedure='initialize_from_hdf5', &
                 message="Unknown density units in .h5 file. Acceptable units are 'g/cc'", file_name=__FILE__, line_number=__LINE__)
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
        call error_msg(module='mod_fluid', class='global_fluid_t', procedure='initialize_from_hdf5', &
     message="Unknown velocity units in .h5 file. Acceptable units are 'km/s' and 'cm/s'", file_name=__FILE__, line_number=__LINE__)
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
        call error_msg(module='mod_fluid', class='global_fluid_t', procedure='initialize_from_hdf5', &
       message="Unknown pressure units in .h5 file. Acceptable units are 'barye', 'Mbar'", file_name=__FILE__, line_number=__LINE__)
      end select
    end if

    call h5%finalize()

    ! Non-dimensionalize
    density = density / rho_0
    x_velocity = x_velocity / v_0
    y_velocity = y_velocity / v_0
    pressure = pressure / p_0

    if(minval(pressure) < tiny(1.0_rk)) then
      call error_msg(module='mod_fluid', class='global_fluid_t', procedure='initialize_from_hdf5', &
                     message="Some (or all) of the pressure array is ~0", file_name=__FILE__, line_number=__LINE__)
    end if

    if(minval(density) < tiny(1.0_rk)) then
      call error_msg(module='mod_fluid', class='global_fluid_t', procedure='initialize_from_hdf5', &
                     message="Some (or all) of the density array is ~0", file_name=__FILE__, line_number=__LINE__)
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
    class(global_fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    real(rk) :: density, x_velocity, y_velocity, pressure

    ! Non-dimensionalize
    density = input%init_density / rho_0
    x_velocity = input%init_x_velocity / v_0
    y_velocity = input%init_y_velocity / v_0
    pressure = input%init_pressure / p_0

    call debug_print('Initializing global_fluid_t from .ini', __FILE__, __LINE__)
    write(*, '(a,4(f0.3, 1x))') 'Initializing with [rho,u,v,p]: ', &
      input%init_density, input%init_x_velocity, input%init_y_velocity, input%init_pressure

    self%rho = density
    self%u = x_velocity
    self%v = y_velocity
    self%p = pressure
    call eos%primitive_to_conserved(rho=self%rho, u=self%u, v=self%v, p=self%p, &
                                    rho_u=self%rho_u, rho_v=self%rho_v, rho_E=self%rho_E)

    if(near_zero(input%init_pressure)) then
      call error_msg(module='mod_fluid', class='global_fluid_t', procedure='initialize_from_ini', &
                     message="Some (or all) of the pressure array is ~0", file_name=__FILE__, line_number=__LINE__)
    end if

    if(near_zero(input%init_density)) then
      call error_msg(module='mod_fluid', class='global_fluid_t', procedure='initialize_from_ini', &
                     message="Some (or all) of the density array is ~0", file_name=__FILE__, line_number=__LINE__)
    end if

  end subroutine initialize_from_ini

  subroutine force_finalization(self)
    class(global_fluid_t), intent(inout) :: self

    call debug_print('Running global_fluid_t%force_finalization()', __FILE__, __LINE__)
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
    type(global_fluid_t), intent(inout) :: self

    call debug_print('Running global_fluid_t%finalize()', __FILE__, __LINE__)
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
    class(global_fluid_t), intent(inout) :: self
    real(rk), intent(in) :: time          !< simulation time
    real(rk), intent(in) :: dt            !< time-step
    integer(ik), intent(in) :: iteration  !< iteration

    call debug_print('Running global_fluid_t%set_time()', __FILE__, __LINE__)
    self%time = time
    self%iteration = iteration
    self%dt = dt
  end subroutine set_time

  subroutine integrate(self, dt, grid, error_code)
    !< Integrate in time
    class(global_fluid_t), intent(inout) :: self
    real(rk), intent(in) :: dt !< time step
    class(grid_t), intent(in) :: grid !< grid class - the solver needs grid topology
    integer(ik), intent(out) :: error_code

    call debug_print('Running global_fluid_t%integerate()', __FILE__, __LINE__)
    self%time = self%time + dt
    self%dt = dt

    select case(trim(self%time_integration_scheme))
    case('ssp_rk2')
      call self%ssp_rk2(grid, error_code)
    case('ssp_rk3')
      call self%ssp_rk3(grid, error_code)
    case default
      call error_msg(module='mod_fluid', class='global_fluid_t', procedure='assign_fluid', &
                     message="Unknown time integration scheme", file_name=__FILE__, line_number=__LINE__)
    end select
  end subroutine integrate

  function time_derivative(self, grid, stage) result(d_dt)
    !< Implementation of the time derivative

    ! Inputs/Output
    class(global_fluid_t), intent(inout) :: self
    class(grid_t), intent(in) :: grid     !< grid class - the solver needs grid topology
    type(local_fluid_t), allocatable :: d_dt    !< dU/dt
    integer(ik), intent(in) :: stage      !< stage in the time integration scheme

    ! Locals
    integer(ik), dimension(2) :: lbounds

    call debug_print('Running global_fluid_t%time_derivative()', __FILE__, __LINE__)

    lbounds = lbound(self%rho)
    allocate(d_dt)
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

    class(global_fluid_t), intent(inout) :: self
    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi

    call debug_print('Running global_fluid_t%calculate_derived_quantities()', __FILE__, __LINE__)
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

  function subtract_local_fluid(lhs, rhs) result(difference)
    !< Implementation of the (-) operator for the fluid type
    class(global_fluid_t), intent(in) :: lhs
    class(local_fluid_t), intent(in) :: rhs
    type(local_fluid_t), allocatable :: difference
  end function subtract_local_fluid

  function add_local_fluid(lhs, rhs) result(sum)
    !< Implementation of the (+) operator for the fluid type
    class(global_fluid_t), intent(in) :: lhs
    class(local_fluid_t), intent(in) :: rhs
    type(local_fluid_t), allocatable :: sum
  end function add_local_fluid

  function local_fluid_mul_real(lhs, rhs) result(product)
    !< Implementation of the fluid * real operation
    class(global_fluid_t), intent(in) :: lhs
    real(rk), intent(in) :: rhs
    type(local_fluid_t), allocatable :: product
  end function local_fluid_mul_real

  function real_mul_local_fluid(lhs, rhs) result(product)
    !< Implementation of the real * fluid operation
    class(global_fluid_t), intent(in) :: rhs
    real(rk), intent(in) :: lhs
    type(local_fluid_t), allocatable :: product
    integer(ik) :: alloc_status
  end function real_mul_local_fluid

  subroutine assign_local_fluid(lhs, rhs)
    !< Implementation of the (=) operator for the fluid type. e.g. lhs = rhs
    class(global_fluid_t), intent(inout) :: lhs
    type(local_fluid_t), intent(in) :: rhs
    integer(ik) :: error_code
    integer(ik) :: alloc_status
  end subroutine assign_local_fluid

  subroutine sanity_check(self, error_code)
    !< Run checks on the conserved variables. Density and pressure need to be > 0. No NaNs or Inifinite numbers either.
    class(global_fluid_t), intent(in) :: self
    integer(ik) :: i, j, ilo, jlo, ihi, jhi
    logical :: invalid_numbers, negative_numbers
    character(len=64) :: err_message = ''
    integer(ik), intent(out) :: error_code

    call debug_print('Running global_fluid_t%sanity_check()', __FILE__, __LINE__)

    error_code = ALL_OK

    negative_numbers = .false.
    invalid_numbers = .false.

    ilo = lbound(self%rho, dim=1) + n_ghost_layers
    ihi = ubound(self%rho, dim=1) - n_ghost_layers
    jlo = lbound(self%rho, dim=2) + n_ghost_layers
    jhi = ubound(self%rho, dim=2) - n_ghost_layers

    !$omp parallel default(none) &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp shared(self, error_code, err_message)
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(self%rho(i, j) < 0.0_rk) then
          write(err_message, '(a, i0, ", ", i0, a)') "Negative density at (", i, j, ")"
          call error_msg(module='mod_fluid', class='global_fluid_t', procedure='sanity_check', &
                         message=err_message, file_name=__FILE__, line_number=__LINE__, error_stop=.false.)
          error_code = NEG_DENSITY
        end if
      end do
    end do
    !$omp end do nowait

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(self%p(i, j) < 0.0_rk) then
          write(err_message, '(a, i0, ", ", i0, a)') "Negative pressure found at (", i, j, ")"
          call error_msg(module='mod_fluid', class='global_fluid_t', procedure='sanity_check', &
                         message=err_message, file_name=__FILE__, line_number=__LINE__, error_stop=.false.)
          error_code = NEG_PRESSURE
        end if
      end do
    end do
    !$omp end do nowait

    ! NaN checks

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%rho(i, j))) then
          write(err_message, '(a, i0, ", ", i0, a)') "NaN density found at (", i, j, ")"
          call error_msg(module='mod_fluid', class='global_fluid_t', procedure='sanity_check', &
                         message=err_message, file_name=__FILE__, line_number=__LINE__, error_stop=.false.)
          error_code = NANS_FOUND
        end if
      end do
    end do
    !$omp end do nowait

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%u(i, j))) then
          write(err_message, '(a, i0, ", ", i0, a)') "NaN x-velocity found at (", i, j, ")"
          call error_msg(module='mod_fluid', class='global_fluid_t', procedure='sanity_check', &
                         message=err_message, file_name=__FILE__, line_number=__LINE__, error_stop=.false.)
          error_code = NANS_FOUND
        end if
      end do
    end do
    !$omp end do nowait
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%v(i, j))) then
          write(err_message, '(a, i0, ", ", i0, a)') "NaN y-velocity found at (", i, j, ")"
          call error_msg(module='mod_fluid', class='global_fluid_t', procedure='sanity_check', &
                         message=err_message, file_name=__FILE__, line_number=__LINE__, error_stop=.false.)
          error_code = NANS_FOUND
        end if
      end do
    end do
    !$omp end do nowait
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        if(ieee_is_nan(self%p(i, j))) then
          write(err_message, '(a, i0, ", ", i0, a)') "NaN pressure found at (", i, j, ")"
          call error_msg(module='mod_fluid', class='global_fluid_t', procedure='sanity_check', &
                         message=err_message, file_name=__FILE__, line_number=__LINE__, error_stop=.false.)
          error_code = NANS_FOUND
        end if
      end do
    end do
    !$omp end do nowait

    !$omp end parallel
  end subroutine sanity_check

  subroutine ssp_rk3(U, grid, error_code)
    !< Strong-stability preserving Runge-Kutta 3rd order
    class(global_fluid_t), intent(inout) :: U
    class(grid_t), intent(in) :: grid
    integer(ik), intent(out) :: error_code

    ! type(local_fluid_t), allocatable :: U_1 !< first stage
    ! type(local_fluid_t), allocatable :: U_2 !< second stage
    ! type(local_fluid_t), allocatable :: R !< hist
    ! real(rk) :: dt

    ! call debug_print('Running global_fluid_t%ssp_rk3()', __FILE__, __LINE__)

    ! dt = U%dt
    ! allocate(U_1, source=U)
    ! allocate(U_2, source=U)
    ! allocate(R, source=U)

    ! ! 1st stage
    ! U_1 = U + U%t(grid, stage=1) * dt
    ! call U_1%residual_smoother()

    ! ! 2nd stage
    ! U_2 = U * (3.0_rk / 4.0_rk) &
    !       + U_1 * (1.0_rk / 4.0_rk) &
    !       + U_1%t(grid, stage=2) * ((1.0_rk / 4.0_rk) * dt)
    ! call U_2%residual_smoother()

    ! ! Final stage
    ! U = U * (1.0_rk / 3.0_rk) &
    !     + U_2 * (2.0_rk / 3.0_rk) &
    !     + U_2%t(grid, stage=3) * ((2.0_rk / 3.0_rk) * dt)
    ! call U%residual_smoother()
    ! call U%sanity_check(error_code)
    ! ! Convergence history
    ! call write_residual_history(first_stage=U_1, last_stage=U)

    ! deallocate(R)
    ! deallocate(U_1)
    ! deallocate(U_2)

  end subroutine ssp_rk3

  subroutine ssp_rk2(U, grid, error_code)
    !< Strong-stability preserving Runge-Kutta 2nd order
    class(global_fluid_t), intent(inout) :: U
    class(grid_t), intent(in) :: grid
    integer(ik), intent(out) :: error_code

    type(local_fluid_t), allocatable :: U_1 !< first stage
    type(local_fluid_t), allocatable :: R !< hist

    ! call debug_print('Running global_fluid_t%rk2()', __FILE__, __LINE__)

    ! allocate(U_1, source=U)
    ! allocate(R, source=U)

    ! ! 1st stage
    ! associate(dt => U%dt)
    !   call debug_print(new_line('a')//'Running global_fluid_t%ssp_rk2_t() 1st stage'//new_line('a'), __FILE__, __LINE__)
    !   U_1 = U + U%t(grid, stage=1) * dt
    !   call U_1%residual_smoother()

    !   ! Final stage
    !   call debug_print(new_line('a')//'Running global_fluid_t%ssp_rk2_t() 2nd stage'//new_line('a'), __FILE__, __LINE__)
    !   U = U * 0.5_rk + U_1 * 0.5_rk + &
    !       U_1%t(grid, stage=2) * (0.5_rk * dt)
    !   call U%residual_smoother()
    !   call U%sanity_check(error_code)

    !   ! ! Convergence history
    !   call write_residual_history(first_stage=U_1, last_stage=U)
    ! end associate

    ! deallocate(R)
    ! deallocate(U_1)

  end subroutine ssp_rk2

  subroutine synchronize_edges(self)
    class(global_fluid_t), intent(inout) :: self
  end subroutine synchronize_edges

  subroutine write_residual_history(first_stage, last_stage)
    !< This writes out the change in residual to a file for convergence history monitoring.

    class(global_fluid_t), intent(in) :: first_stage
    class(global_fluid_t), intent(in) :: last_stage

    real(rk) :: rho_diff   !< difference in the rho residual
    real(rk) :: rho_u_diff !< difference in the rhou residual
    real(rk) :: rho_v_diff !< difference in the rhov residual
    real(rk) :: rho_E_diff !< difference in the rhoE residual
    integer(ik) :: io
    integer(ik) :: i, j, ilo, jlo, ihi, jhi

    call debug_print('Running global_fluid_t%write_residual_history()', __FILE__, __LINE__)

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

  function gather(self, image) result(gathered_data)
    !< Gather all the conserved vars onto one image for I/O

    class(global_fluid_t), intent(in) :: self
    integer(ik), intent(in) :: image !< image number
    real(rk), dimension(4, &
                        self%cell_dims_global(1):self%cell_dims_global(2), &
                        self%cell_dims_global(3):self%cell_dims_global(4)) :: gathered_data

    real(rk), dimension(:, :, :), allocatable :: gather_coarray[:]

    ! associate(ilo => self%cell_dims_global(1), ihi => self%cell_dims_global(2), &
    !           jlo => self%cell_dims_global(3), jhi => self%cell_dims_global(4))
    !   allocate(gather_coarray(4, ilo:ihi, jlo:jhi)[*])
    ! end associate

    ! associate(ilo => self%cell_dims_img(1), ihi => self%cell_dims_img(2), &
    !           jlo => self%cell_dims_img(3), jhi => self%cell_dims_img(4))

    !   gather_coarray(4, ilo:ihi, jlo:jhi)[image] = self%conserved_vars(4, ilo:ihi, jlo:jhi)
    !   rho
    !   rho_u
    !   rho_v
    !   rho_E
    !   u
    !   v
    !   p
    !   cs
    !   mach_u
    !   mach_v

    !   sync all
    !   if(this_image() == image) gathered_data = gather_coarray
    ! end associate

    deallocate(gather_coarray)

  end function gather

  subroutine get_image_partitioning(self, partitioning)
    !< For the sake of visualization, make an array that has the image number for
    !< each cell
    class(global_fluid_t), intent(in) :: self
    integer(ik), dimension(self%cell_dims_global(1):self%cell_dims_global(2), &
                           self%cell_dims_global(3):self%cell_dims_global(4)), intent(out) :: partitioning
    integer(ik), dimension(:, :), allocatable :: gather_coarray[:]

    associate(ilo => self%cell_dims_global(1), ihi => self%cell_dims_global(2), &
              jlo => self%cell_dims_global(3), jhi => self%cell_dims_global(4))
      allocate(gather_coarray(ilo:ihi, jlo:jhi)[*])
    end associate

    associate(ilo => self%cell_dims_img(1), ihi => self%cell_dims_img(2), &
              jlo => self%cell_dims_img(3), jhi => self%cell_dims_img(4))
      gather_coarray(ilo:ihi, jlo:jhi)[this_image()] = this_image()
      sync all
    end associate

    if(this_image() == 1) partitioning = gather_coarray
  end subroutine get_image_partitioning
end module mod_global_fluid
