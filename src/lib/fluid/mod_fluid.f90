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

module mod_fluid
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use, intrinsic :: ieee_arithmetic

  use mod_field, only: field_2d_t, field_2d
  use mod_error, only: ALL_OK, NEG_DENSITY, NEG_PRESSURE, NANS_FOUND, error_msg
  use mod_globals, only: enable_debug_print, debug_print, print_evolved_cell_data, print_recon_data, n_ghost_layers
  use mod_nondimensionalization, only: scale_factors_set, rho_0, v_0, p_0, e_0, t_0
  use mod_floating_point_utils, only: near_zero, nearly_equal, neumaier_sum, neumaier_sum_2, neumaier_sum_3, neumaier_sum_4
  use mod_functional, only: operator(.sort.)
  use mod_units
  use mod_distinguisher, only: distinguish
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_bc_factory, only: bc_factory
  use mod_grid_block, only: grid_block_t
  use mod_grid_block_2d, only: grid_block_2d_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use hdf5_interface, only: hdf5_file
  use mod_flux_solver, only: flux_solver_t
  ! use mod_ausm_plus_solver, only: ausm_plus_solver_t
  ! use mod_fvleg_solver, only: fvleg_solver_t
  ! use mod_m_ausmpw_plus_solver, only: m_ausmpw_plus_solver_t
  use mod_ausmpw_plus_solver, only: ausmpw_plus_solver_t
  ! use mod_slau_solver, only: slau_solver_t

  implicit none

  private
  public :: fluid_t, new_fluid

  logical, parameter :: filter_small_mach = .false.

  type :: fluid_t
    !< Fluid solver physics package

    private ! make all private by default

    type(field_2d_t), public :: rho    !< (i, j); Conserved quantities
    type(field_2d_t), public :: rho_u  !< (i, j); Conserved quantities
    type(field_2d_t), public :: rho_v  !< (i, j); Conserved quantities
    type(field_2d_t), public :: rho_E  !< (i, j); Conserved quantities
    type(field_2d_t), public :: u      !< (i, j); Conserved quantities
    type(field_2d_t), public :: v      !< (i, j); Conserved quantities
    type(field_2d_t), public :: p      !< (i, j); Conserved quantities
    type(field_2d_t), public :: cs     !< (i, j); Conserved quantities
    type(field_2d_t), public :: mach_u   !< (i, j); Conserved quantities
    type(field_2d_t), public :: mach_v   !< (i, j); Conserved quantities

    class(flux_solver_t), allocatable :: solver !< solver scheme used to flux quantities at cell interfaces

    ! Time variables
    character(len=10) :: time_integration_scheme = 'ssp_rk2'
    real(rk) :: time = 0.0_rk !< current simulation time
    real(rk) :: dt = 0.0_rk   !< time step
    real(rk) :: cfl = 0.0_rk  !< Courant–Friedrichs–Lewy condition (CFL)
    integer(ik) :: iteration = 0 !< current iteration number

    logical, public :: prim_vars_updated = .false.
    logical :: smooth_residuals = .true.

    ! Residual history
    character(len=32) :: residual_hist_file = 'residual_hist.csv'
    logical :: residual_hist_header_written = .false.

    class(boundary_condition_t), allocatable :: bc_plus_x
    class(boundary_condition_t), allocatable :: bc_plus_y
    class(boundary_condition_t), allocatable :: bc_minus_x
    class(boundary_condition_t), allocatable :: bc_minus_y

  contains
    ! Private methods
    private
    procedure :: initialize_from_ini
    procedure :: initialize_from_hdf5
    procedure :: residual_smoother
    procedure :: calculate_derived_quantities
    procedure :: sanity_check
    procedure :: apply_bc
    procedure :: ssp_rk2
    procedure :: ssp_rk3
    procedure :: get_continuity_sensor

    ! Operators
    procedure, pass(lhs), public :: add_fluid
    procedure, pass(lhs), public :: subtract_fluid
    procedure, pass(lhs), public :: fluid_mul_real
    procedure, pass(rhs), public :: real_mul_fluid
    procedure, pass(lhs), public :: assign_fluid

    ! Public methods
    procedure, public :: initialize
    procedure, public :: set_time
    procedure, public :: get_timestep
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
    class(grid_block_t), intent(in) :: grid
    type(fluid_t), pointer :: fluid

    select type(grid)
    class is (grid_block_2d_t)
      allocate(fluid)
      call fluid%initialize(input, grid)
    end select
  end function new_fluid

  subroutine initialize(self, input, grid)
    class(fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_block_2d_t), intent(in) :: grid
    class(flux_solver_t), pointer :: solver => null()
    class(boundary_condition_t), pointer :: bc => null()

    integer(ik) :: alloc_status, i, j, ilo, ihi, jlo, jhi, io

    alloc_status = 0
    if(enable_debug_print) call debug_print('Calling fluid_t%initialize()', __FILE__, __LINE__)

    if(.not. scale_factors_set) then
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='initialize', &
                     message="Global non-dimensional scale factors haven't been set yet. "// &
                     "These need to be set before fluid initialization", &
                     file_name=__FILE__, line_number=__LINE__)
    end if

    self%rho = field_2d(name='rho', long_name='Density', &
                        descrip='Cell Density', units='g/cm^3', &
                        global_dims=grid%cell_dim_global, n_halo_cells=input%n_ghost_layers)

    self%rho_u = field_2d(name='rhou', long_name='rhou', descrip='Cell Conserved quantity (Density * X-Velocity)', &
                          units='g cm/cm^2 s', &
                          global_dims=grid%cell_dim_global, n_halo_cells=input%n_ghost_layers)

    self%rho_v = field_2d(name='rhov', long_name='rhov', descrip='Cell Conserved quantity (Density * Y-Velocity)', &
                          units='g cm/cm^2 s', &
                          global_dims=grid%cell_dim_global, n_halo_cells=input%n_ghost_layers)

    self%rho_E = field_2d(name='rhoE', long_name='rhoE', descrip='Cell Conserved quantity (Density * Total Energy)', &
                          units='g erg / cm^3', &
                          global_dims=grid%cell_dim_global, n_halo_cells=input%n_ghost_layers)

    self%u = field_2d(name='u', long_name='X Velocity', descrip='Cell X-Velocity', units='cm/s', &
                      global_dims=grid%cell_dim_global, n_halo_cells=input%n_ghost_layers)

    self%v = field_2d(name='v', long_name='Y Velocity', descrip='Cell Y-Velocity', units='cm/s', &
                      global_dims=grid%cell_dim_global, n_halo_cells=input%n_ghost_layers)

    self%p = field_2d(name='p', long_name='Pressure', descrip='Cell Pressure', units='barye', &
                      global_dims=grid%cell_dim_global, n_halo_cells=input%n_ghost_layers)

    self%cs = field_2d(name='cs', long_name='Sound Speed', descrip='Cell Sound Speed', units='cm/s', &
                       global_dims=grid%cell_dim_global, n_halo_cells=input%n_ghost_layers)

    self%mach_u = field_2d(name='mach_u', long_name='Mach X', descrip='Cell Mach number in x-direction', &
                           units='dimensionless', &
                           global_dims=grid%cell_dim_global, n_halo_cells=input%n_ghost_layers)

    self%mach_v = field_2d(name='mach_v', long_name='Mach Y', descrip='Cell Mach number in y-direction', &
                           units='dimensionless', &
                           global_dims=grid%cell_dim_global, n_halo_cells=input%n_ghost_layers)

    self%smooth_residuals = input%smooth_residuals

    self%time_integration_scheme = trim(input%time_integration_strategy)

    select case(trim(input%flux_solver))
      ! case('FVLEG')
      !   allocate(fvleg_solver_t :: solver)
      ! case('AUSM+-u', 'AUSM+-up', 'AUSM+-up_all_speed')
      !   error stop "There are issues in the AUSM+ solver for now; exiting..."
      !   allocate(ausm_plus_solver_t :: solver)
      ! case('M-AUSMPW+')
      !   allocate(m_ausmpw_plus_solver_t :: solver)
    case('AUSMPW+')
      allocate(ausmpw_plus_solver_t :: solver)
      ! case('SLAU', 'SLAU2', 'SD-SLAU', 'SD-SLAU2')
      !   allocate(slau_solver_t :: solver)
    case default
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='initialize', &
                     message="Invalid flux solver. It must be one of the following: "// &
                     "['FVLEG', 'AUSM+-u','AUSM+-a','AUSM+-up','AUSM+-up_all_speed', "// &
                     "'AUSMPW+', 'M-AUSMPW+', 'SLAU', 'SLAU2', 'SD-SLAU', 'SD-SLAU2'], "// &
                     "the input was: '"//trim(input%flux_solver)//"'", &
                     file_name=__FILE__, line_number=__LINE__)
    end select

    call solver%initialize(input)
    allocate(self%solver, source=solver)
    deallocate(solver)

    ! Set boundary conditions
    bc => bc_factory(bc_type=input%plus_x_bc, location='+x', input=input, grid=grid)
    allocate(self%bc_plus_x, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate bc_plus_x"
    deallocate(bc)

    bc => bc_factory(bc_type=input%plus_y_bc, location='+y', input=input, grid=grid)
    allocate(self%bc_plus_y, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate bc_plus_y"
    deallocate(bc)

    bc => bc_factory(bc_type=input%minus_x_bc, location='-x', input=input, grid=grid)
    allocate(self%bc_minus_x, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate bc_minus_x"
    deallocate(bc)

    bc => bc_factory(bc_type=input%minus_y_bc, location='-y', input=input, grid=grid)
    allocate(self%bc_minus_y, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate bc_minus_y"
    deallocate(bc)

    open(newunit=io, file=trim(self%residual_hist_file), status='replace')
    write(io, '(a)') 'iteration,time,rho,rho_u,rho_v,rho_E'
    close(io)

    if(input%read_init_cond_from_file .or. input%restart_from_file) then
      call self%initialize_from_hdf5(input)
    else
      call self%initialize_from_ini(input)
    end if

    call eos%sound_speed(p=self%p, rho=self%rho, cs=self%cs)
    self%mach_v = self%v / self%cs
    self%mach_u = self%u / self%cs
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

    if(enable_debug_print) call debug_print('Initializing fluid_t from hdf5', __FILE__, __LINE__)

    if(input%restart_from_file) then
      filename = trim(input%restart_file)
    else
      filename = trim(input%initial_condition_file)
    end if

    file_exists = .false.
    inquire(file=filename, exist=file_exists)

    if(.not. file_exists) then
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='initialize_from_hdf5', &
                     message='File not found: "'//filename//'"', file_name=__FILE__, line_number=__LINE__)
    end if

    call h5%initialize(filename=filename, status='old', action='r')

    call h5%get('/density', density)
    call h5%readattr('/density', 'units', str_buff)
    ! self%rho%units = trim(str_buff)

    if(input%restart_from_file) then
      select case(trim(str_buff))
      case('g/cc')
      case default
        call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='initialize_from_hdf5', &
                       message="Unknown density units in .h5 file. Acceptable units are 'g/cc'", &
                       file_name=__FILE__, line_number=__LINE__)
      end select
    end if

    call h5%get('/x_velocity', x_velocity)
    call h5%get('/y_velocity', y_velocity)
    call h5%readattr('/x_velocity', 'units', str_buff)
    ! self%u%units = trim(str_buff)
    ! self%v%units = trim(str_buff)

    if(input%restart_from_file) then
      select case(trim(str_buff))
      case('km/s')
        x_velocity = x_velocity * km_per_sec_to_cm_per_sec
        y_velocity = y_velocity * km_per_sec_to_cm_per_sec
      case('cm/s')
      case default
        call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='initialize_from_hdf5', &
                       message="Unknown velocity units in .h5 file. Acceptable units are 'km/s' and 'cm/s'", &
                       file_name=__FILE__, line_number=__LINE__)
      end select
    end if

    call h5%get('/pressure', pressure)
    call h5%readattr('/pressure', 'units', str_buff)

    if(input%restart_from_file) then
      select case(trim(str_buff))
      case('barye')
      case('Mbar')
        pressure = pressure * mega_bar_to_barye
      case default
        call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='initialize_from_hdf5', &
                       message="Unknown pressure units in .h5 file. Acceptable units are 'barye', 'Mbar'", &
                       file_name=__FILE__, line_number=__LINE__)
      end select
    end if

    call h5%finalize()

    associate(ilo => self%rho%lbounds(1), ihi => self%rho%ubounds(1), &
              jlo => self%rho%lbounds(2), jhi => self%rho%ubounds(2))
      self%rho%data(ilo:ihi, jlo:jhi) = density(ilo:ihi, jlo:jhi)
      self%u%data(ilo:ihi, jlo:jhi) = x_velocity(ilo:ihi, jlo:jhi)
      self%v%data(ilo:ihi, jlo:jhi) = y_velocity(ilo:ihi, jlo:jhi)
      self%p%data(ilo:ihi, jlo:jhi) = pressure(ilo:ihi, jlo:jhi)
    end associate

    call self%rho%make_non_dimensional(non_dim_factor=rho_0)
    call self%u%make_non_dimensional(non_dim_factor=v_0)
    call self%v%make_non_dimensional(non_dim_factor=v_0)
    call self%p%make_non_dimensional(non_dim_factor=p_0)

    if(minval(pressure) < tiny(1.0_rk)) then
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='initialize_from_hdf5', &
                     message="Some (or all) of the pressure array is ~0", file_name=__FILE__, line_number=__LINE__)
    end if

    if(minval(density) < tiny(1.0_rk)) then
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='initialize_from_hdf5', &
                     message="Some (or all) of the density array is ~0", file_name=__FILE__, line_number=__LINE__)
    end if

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

    error stop "Fix fluid_t%initialize_from_ini"
    ! ! Non-dimensionalize
    ! density = input%init_density / rho_0
    ! x_velocity = input%init_x_velocity / v_0
    ! y_velocity = input%init_y_velocity / v_0
    ! pressure = input%init_pressure / p_0

    ! if (enable_debug_print) call debug_print('Initializing fluid_t from .ini', __FILE__, __LINE__)
    ! write(*, '(a,4(f0.3, 1x))') 'Initializing with [rho,u,v,p]: ', &
    !   input%init_density, input%init_x_velocity, input%init_y_velocity, input%init_pressure

    ! self%rho = density
    ! self%u = x_velocity
    ! self%v = y_velocity
    ! self%p = pressure
    ! call eos%primitive_to_conserved(rho=self%rho, u=self%u, v=self%v, p=self%p, &
    !                                 rho_u=self%rho_u, rho_v=self%rho_v, rho_E=self%rho_E)

    ! if(near_zero(input%init_pressure)) then
    !   call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='initialize_from_ini', &
    !                  message="Some (or all) of the pressure array is ~0", file_name=__FILE__, line_number=__LINE__)
    ! end if

    ! if(near_zero(input%init_density)) then
    !   call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='initialize_from_ini', &
    !                  message="Some (or all) of the density array is ~0", file_name=__FILE__, line_number=__LINE__)
    ! end if

  end subroutine initialize_from_ini

  subroutine force_finalization(self)
    class(fluid_t), intent(inout) :: self

    if(enable_debug_print) call debug_print('Running fluid_t%force_finalization()', __FILE__, __LINE__)
    ! if(allocated(self%rho)) deallocate(self%rho)
    ! if(allocated(self%u)) deallocate(self%u)
    ! if(allocated(self%v)) deallocate(self%v)
    ! if(allocated(self%p)) deallocate(self%p)
    ! if(allocated(self%rho_u)) deallocate(self%rho_u)
    ! if(allocated(self%rho_v)) deallocate(self%rho_v)
    ! if(allocated(self%rho_E)) deallocate(self%rho_E)
    ! if(allocated(self%cs)) deallocate(self%cs)
    ! if(allocated(self%mach_u)) deallocate(self%mach_u)
    ! if(allocated(self%mach_v)) deallocate(self%mach_v)
    ! if(allocated(self%time_integrator)) deallocate(self%time_integrator)
  end subroutine force_finalization

  subroutine finalize(self)
    type(fluid_t), intent(inout) :: self

    if(enable_debug_print) call debug_print('Running fluid_t%finalize()', __FILE__, __LINE__)
    ! if(allocated(self%rho)) deallocate(self%rho)
    ! if(allocated(self%u)) deallocate(self%u)
    ! if(allocated(self%v)) deallocate(self%v)
    ! if(allocated(self%p)) deallocate(self%p)
    ! if(allocated(self%rho_u)) deallocate(self%rho_u)
    ! if(allocated(self%rho_v)) deallocate(self%rho_v)
    ! if(allocated(self%rho_E)) deallocate(self%rho_E)
    ! if(allocated(self%cs)) deallocate(self%cs)
    ! if(allocated(self%mach_u)) deallocate(self%mach_u)
    ! if(allocated(self%mach_v)) deallocate(self%mach_v)
    ! if(allocated(self%solver)) deallocate(self%solver)
  end subroutine finalize

  subroutine set_time(self, time, iteration)
    !< Set the time statistics
    class(fluid_t), intent(inout) :: self
    real(rk), intent(in) :: time          !< simulation time
    integer(ik), intent(in) :: iteration  !< iteration

    if(enable_debug_print) call debug_print('Running fluid_t%set_time()', __FILE__, __LINE__)
    self%time = time
    self%iteration = iteration
  end subroutine set_time

  subroutine integrate(self, dt, grid, error_code)
    !< Integrate in time
    class(fluid_t), intent(inout) :: self
    real(rk), intent(in) :: dt !< time step
    class(grid_block_t), intent(in) :: grid !< grid class - the solver needs grid topology
    integer(ik), intent(out) :: error_code
    
    if(enable_debug_print) call debug_print('Running fluid_t%integerate()', __FILE__, __LINE__)
    self%time = self%time + dt
    self%dt = dt

    select type(grid)
    class is (grid_block_2d_t)
      select case(trim(self%time_integration_scheme))
      case('ssp_rk2')
        call self%ssp_rk2(grid, error_code)
      case('ssp_rk3')
        call self%ssp_rk3(grid, error_code)
      case default
        call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='assign_fluid', &
                      message="Unknown time integration scheme", file_name=__FILE__, line_number=__LINE__)
      end select
    end select
  end subroutine integrate

  function time_derivative(self, grid, stage) result(d_dt)
    !< Implementation of the time derivative

    ! Inputs/Output
    class(fluid_t), intent(inout) :: self
    class(grid_block_2d_t), intent(in) :: grid  !< grid class - the solver needs grid topology
    type(fluid_t), allocatable :: d_dt          !< dU/dt
    integer(ik), intent(in) :: stage            !< stage in the time integration scheme

    if(enable_debug_print) call debug_print('Running fluid_t%time_derivative()', __FILE__, __LINE__)

    allocate(d_dt, source=self)
    call self%solver%solve(grid=grid, &
                           rho=self%rho, u=self%u, v=self%v, p=self%p, &
                           d_rho_dt=d_dt%rho, &
                           d_rho_u_dt=d_dt%rho_u, &
                           d_rho_v_dt=d_dt%rho_v, &
                           d_rho_E_dt=d_dt%rho_E)

  end function time_derivative

  real(rk) function get_timestep(self, grid) result(delta_t)
    class(fluid_t), intent(inout) :: self
    class(grid_block_t), intent(in) :: grid
    real(rk), save :: coarray_max_delta_t[*] !< the max allowable timestep on each subdomain/image
    integer(ik) :: ierr
    integer(ik) :: ilo, ihi, jlo, jhi ! fluid lo/hi indicies
    integer(ik) :: g_ilo, g_ihi, g_jlo, g_jhi ! grid lo/hi indices
    character(len=200) :: err_msg

    real(rk), allocatable, save, dimension(:, :) :: dx, dy

    g_ilo = grid%cell_lbounds(1)
    g_ihi = grid%cell_ubounds(1)

    g_jlo = grid%cell_lbounds(2)
    g_jhi = grid%cell_ubounds(2)

    if (.not. allocated(dx)) allocate(dx(g_ilo:g_ihi, g_jlo:g_jhi))
    if (.not. allocated(dy)) allocate(dy(g_ilo:g_ihi, g_jlo:g_jhi))
  
    select type(grid)
    class is (grid_block_2d_t)
      dy = grid%dy(g_ilo:g_ihi, g_jlo:g_jhi)
      dx = grid%dy(g_ilo:g_ihi, g_jlo:g_jhi)
    end select

    err_msg = ''

    ilo = self%u%lbounds(1)
    ihi = self%u%ubounds(1)

    jlo = self%u%lbounds(2)
    jhi = self%u%ubounds(2)
    print*,'fluid: ', ilo, ihi, jlo, jhi
    print*,'grid : ', g_ilo, g_ihi, g_jlo, g_jhi

    coarray_max_delta_t = minval(self%cfl / &
                       (((abs(self%u%data(ilo:ihi, jlo:jhi)) + &
                          self%cs%data(ilo:ihi, jlo:jhi)) / dx) + &
                        ((abs(self%v%data(ilo:ihi, jlo:jhi)) + &
                          self%cs%data(ilo:ihi, jlo:jhi)) / dy)))
    
    sync all
    ! Get the minimum timestep across all the images and save it on each image
    call co_min(coarray_max_delta_t, stat=ierr)
    if (ierr /= 0) then
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='get_timestep', &
                     message="Unable to run the co_min() to get min timestep across images: err_msg=" // trim(err_msg), &
                     file_name=__FILE__, line_number=__LINE__)
    end if
    delta_t = coarray_max_delta_t
  end function get_timestep

  subroutine calculate_derived_quantities(self)
    !< Find derived quantities like sound speed, mach number, primitive variables

    class(fluid_t), intent(inout) :: self
    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi

    if(enable_debug_print) call debug_print('Running fluid_t%calculate_derived_quantities()', __FILE__, __LINE__)
    call eos%conserved_to_primitive(rho=self%rho, rho_u=self%rho_u, &
                                    rho_v=self%rho_v, rho_E=self%rho_E, &
                                    u=self%u, v=self%v, p=self%p)
    call eos%sound_speed(p=self%p, rho=self%rho, cs=self%cs)

    self%mach_u = self%u / self%cs
    self%mach_v = self%v / self%cs
    self%prim_vars_updated = .true.
  end subroutine calculate_derived_quantities

  subroutine get_continuity_sensor(self)
    !< Run a distinguishing step that determines which regions in the domain are smooth (continuous),
    !< or discontinuous (linear or non-linear)
    class(fluid_t), intent(inout) :: self

    ! call distinguish(lbounds=lbound(self%rho), rho=self%rho, u=self%u, v=self%v, p=self%p, continuity_sensor=self%continuous_sensor)
  end subroutine get_continuity_sensor

  subroutine residual_smoother(self)
    class(fluid_t), intent(inout) :: self

    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: i, j
    real(rk), parameter :: EPS = 5e-14_rk

    ! if (enable_debug_print) call debug_print('Running fluid_t%residual_smoother()', __FILE__, __LINE__)
    ! if(self%smooth_residuals) then
    !   ilo = lbound(self%rho, dim=1) + 1
    !   ihi = ubound(self%rho, dim=1) - 1
    !   jlo = lbound(self%rho, dim=2) + 1
    !   jhi = ubound(self%rho, dim=2) - 1

    !   !$omp parallel default(none) &
    !   !$omp shared(self) &
    !   !$omp firstprivate(ilo,ihi,jlo,jhi) &
    !   !$omp private(i,j)
    !   !$omp do
    !   do j = jlo, jhi
    !     do i = ilo, ihi
    !       if(abs(self%rho(i, j) - self%rho(i - 1, j)) < EPS) self%rho(i, j) = self%rho(i - 1, j)
    !       if(abs(self%rho_u(i, j) - self%rho_u(i - 1, j)) < EPS) self%rho_u(i, j) = self%rho_u(i - 1, j)
    !       if(abs(self%rho_v(i, j) - self%rho_v(i - 1, j)) < EPS) self%rho_v(i, j) = self%rho_v(i - 1, j)
    !       if(abs(self%rho_E(i, j) - self%rho_E(i - 1, j)) < EPS) self%rho_E(i, j) = self%rho_E(i - 1, j)
    !     end do
    !   end do
    !   !$omp end do
    !   !$omp do
    !   do j = jlo, jhi
    !     do i = ilo, ihi
    !       if(abs(self%rho(i, j) - self%rho(i, j - 1)) < EPS) self%rho(i, j) = self%rho(i, j - 1)
    !       if(abs(self%rho_u(i, j) - self%rho_u(i, j - 1)) < EPS) self%rho_u(i, j) = self%rho_u(i, j - 1)
    !       if(abs(self%rho_v(i, j) - self%rho_v(i, j - 1)) < EPS) self%rho_v(i, j) = self%rho_v(i, j - 1)
    !       if(abs(self%rho_E(i, j) - self%rho_E(i, j - 1)) < EPS) self%rho_E(i, j) = self%rho_E(i, j - 1)
    !     end do
    !   end do
    !   !$omp end do
    !   !$omp end parallel
    ! end if

  end subroutine residual_smoother

  function subtract_fluid(lhs, rhs) result(difference)
    !< Implementation of the (-) operator for the fluid type
    class(fluid_t), intent(in) :: lhs
    class(fluid_t), intent(in) :: rhs
    type(fluid_t), allocatable :: difference

    if(enable_debug_print) call debug_print('Running fluid_t%subtract_fluid()', __FILE__, __LINE__)

    allocate(difference, source=lhs)
    difference%rho = lhs%rho - rhs%rho
    difference%rho_u = lhs%rho_u - rhs%rho_u
    difference%rho_v = lhs%rho_v - rhs%rho_v
    difference%rho_E = lhs%rho_E - rhs%rho_E
    difference%prim_vars_updated = .false.
  end function subtract_fluid

  function add_fluid(lhs, rhs) result(sum)
    !< Implementation of the (+) operator for the fluid type
    class(fluid_t), intent(in) :: lhs
    class(fluid_t), intent(in) :: rhs
    type(fluid_t), allocatable :: sum

    if(enable_debug_print) call debug_print('Running fluid_t%add_fluid()', __FILE__, __LINE__)

    allocate(sum, source=lhs)
    sum%rho = lhs%rho + rhs%rho
    sum%rho_u = lhs%rho_u + rhs%rho_u
    sum%rho_v = lhs%rho_v + rhs%rho_v
    sum%rho_E = lhs%rho_E + rhs%rho_E
    sum%prim_vars_updated = .false.
  end function add_fluid

  function fluid_mul_real(lhs, rhs) result(product)
    !< Implementation of the fluid * real operation
    class(fluid_t), intent(in) :: lhs
    real(rk), intent(in) :: rhs
    type(fluid_t), allocatable :: product

    integer(ik) :: alloc_status

    if(enable_debug_print) call debug_print('Running fluid_t%fluid_mul_real()', __FILE__, __LINE__)

    allocate(product, source=lhs, stat=alloc_status)
    if(alloc_status /= 0) then
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='fluid_mul_real', &
                     message="Unable to allocate lhs%solver", file_name=__FILE__, line_number=__LINE__)
    end if

    product%rho = lhs%rho * rhs
    product%rho_u = lhs%rho_u * rhs
    product%rho_v = lhs%rho_v * rhs
    product%rho_E = lhs%rho_E * rhs

    product%prim_vars_updated = .false.
  end function fluid_mul_real

  function real_mul_fluid(lhs, rhs) result(product)
    !< Implementation of the real * fluid operation
    class(fluid_t), intent(in) :: rhs
    real(rk), intent(in) :: lhs
    type(fluid_t), allocatable :: product
    integer(ik) :: alloc_status

    if(enable_debug_print) call debug_print('Running fluid_t%real_mul_fluid()', __FILE__, __LINE__)

    allocate(product, source=rhs, stat=alloc_status)
    if(alloc_status /= 0) then
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='real_mul_fluid', &
                     message="Unable to allocate product", file_name=__FILE__, line_number=__LINE__)
    end if

    product%rho = lhs * rhs%rho
    product%rho_u = lhs * rhs%rho_u
    product%rho_v = lhs * rhs%rho_v
    product%rho_E = lhs * rhs%rho_E
    product%prim_vars_updated = .false.
  end function real_mul_fluid

  subroutine assign_fluid(lhs, rhs)
    !< Implementation of the (=) operator for the fluid type. e.g. lhs = rhs
    class(fluid_t), intent(inout) :: lhs
    type(fluid_t), intent(in) :: rhs
    integer(ik) :: alloc_status

    if(enable_debug_print) call debug_print('Running fluid_t%assign_fluid()', __FILE__, __LINE__)
    lhs%residual_hist_header_written = rhs%residual_hist_header_written
    lhs%rho = rhs%rho
    lhs%rho_u = rhs%rho_u
    lhs%rho_v = rhs%rho_v
    lhs%rho_E = rhs%rho_E
    if(allocated(lhs%solver)) deallocate(lhs%solver)

    allocate(lhs%solver, source=rhs%solver, stat=alloc_status)
    if(alloc_status /= 0) then
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='assign_fluid', &
                     message="Unable to allocate lhs%solver", file_name=__FILE__, line_number=__LINE__)
    end if

    lhs%time_integration_scheme = rhs%time_integration_scheme
    lhs%time = rhs%time
    lhs%dt = rhs%dt
    lhs%iteration = rhs%iteration
    lhs%prim_vars_updated = rhs%prim_vars_updated
    lhs%smooth_residuals = rhs%smooth_residuals
    lhs%residual_hist_file = rhs%residual_hist_file
    lhs%residual_hist_header_written = rhs%residual_hist_header_written
    call lhs%calculate_derived_quantities()
  end subroutine assign_fluid

  subroutine sanity_check(self, error_code)
    !< Run checks on the conserved variables. Density and pressure need to be > 0. No NaNs or Inifinite numbers either.
    class(fluid_t), intent(in) :: self
    logical :: invalid_numbers, negative_numbers
    integer(ik), intent(out) :: error_code

    if(enable_debug_print) call debug_print('Running fluid_t%sanity_check()', __FILE__, __LINE__)

    error_code = ALL_OK

    invalid_numbers = self%rho%has_nans()
    if(invalid_numbers) then
      error_code = NANS_FOUND
      return
    end if

    invalid_numbers = self%u%has_nans()
    if(invalid_numbers) then
      error_code = NANS_FOUND
      return
    end if

    invalid_numbers = self%v%has_nans()
    if(invalid_numbers) then
      error_code = NANS_FOUND
      return
    end if

    invalid_numbers = self%p%has_nans()
    if(invalid_numbers) then
      error_code = NANS_FOUND
      return
    end if

    negative_numbers = self%p%has_negatives()
    if(negative_numbers) then
      error_code = NEG_PRESSURE
      return
    end if

    negative_numbers = self%rho%has_negatives()
    if(negative_numbers) then
      error_code = NEG_DENSITY
      return
    end if

  end subroutine sanity_check

  subroutine ssp_rk3(U, grid, error_code)
    !< Strong-stability preserving Runge-Kutta 3rd order
    class(fluid_t), intent(inout) :: U
    class(grid_block_2d_t), intent(in) :: grid

    type(fluid_t), allocatable :: U_1 !< first stage
    type(fluid_t), allocatable :: U_2 !< second stage
    integer(ik), intent(out) :: error_code
    real(rk) :: dt

    if(enable_debug_print) call debug_print('Running fluid_t%ssp_rk3()', __FILE__, __LINE__)

    dt = U%dt
    allocate(U_1, source=U)
    allocate(U_2, source=U)

    ! 1st stage
    if(enable_debug_print) call debug_print(new_line('a')//'Running fluid_t%ssp_rk3() 1st stage'//new_line('a'), __FILE__, __LINE__)
    call U%apply_bc()
    U_1 = U + U%t(grid, stage=1) * dt

    ! 2nd stage
    if(enable_debug_print) call debug_print(new_line('a')//'Running fluid_t%ssp_rk3() 2nd stage'//new_line('a'), __FILE__, __LINE__)
    call U_1%apply_bc()
    U_2 = U * (3.0_rk / 4.0_rk) &
          + U_1 * (1.0_rk / 4.0_rk) &
          + U_1%t(grid, stage=2) * ((1.0_rk / 4.0_rk) * dt)

    ! Final stage
    if(enable_debug_print) call debug_print(new_line('a')//'Running fluid_t%ssp_rk2() 3rd stage'//new_line('a'), __FILE__, __LINE__)
    call U_2%apply_bc()
    U = U * (1.0_rk / 3.0_rk) &
        + U_2 * (2.0_rk / 3.0_rk) &
        + U_2%t(grid, stage=3) * ((2.0_rk / 3.0_rk) * dt)

    call U%sanity_check(error_code)

    ! Convergence history
    call write_residual_history(first_stage=U_1, last_stage=U)

    deallocate(U_1)
    deallocate(U_2)

  end subroutine ssp_rk3

  subroutine ssp_rk2(U, grid, error_code)
    !< Strong-stability preserving Runge-Kutta 2nd order
    class(fluid_t), intent(inout) :: U
    class(grid_block_2d_t), intent(in) :: grid

    type(fluid_t), allocatable :: U_1 !< first stage
    integer(ik), intent(out) :: error_code
    real(rk) :: dt

    if(enable_debug_print) call debug_print('Running fluid_t%rk2()', __FILE__, __LINE__)

    dt = U%dt
    allocate(U_1, source=U)

    ! 1st stage
    if(enable_debug_print) then
      call debug_print(new_line('a')//'Running fluid_t%ssp_rk2_t() 1st stage'//new_line('a'), &
                       __FILE__, __LINE__)
    end if
    call U%apply_bc()
    U_1 = U + U%t(grid, stage=1) * dt

    ! Final stage
    if(enable_debug_print) then
      call debug_print(new_line('a')//'Running fluid_t%ssp_rk2_t() 2nd stage'//new_line('a'), &
                       __FILE__, __LINE__)
    end if
    call U_1%apply_bc()
    U = U * 0.5_rk + U_1 * 0.5_rk + U_1%t(grid, stage=2) * (0.5_rk * dt)

    call U%sanity_check(error_code)

    ! Convergence history
    call write_residual_history(first_stage=U_1, last_stage=U)

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
    integer(ik) :: ilo, jlo, ihi, jhi

    if(enable_debug_print) call debug_print('Running fluid_t%write_residual_history()', __FILE__, __LINE__)

    ilo = last_stage%rho%lbounds(1)
    ihi = last_stage%rho%ubounds(1)
    jlo = last_stage%rho%lbounds(2)
    jhi = last_stage%rho%ubounds(2)

    rho_diff = maxval(abs(last_stage%rho%data(ilo:ihi, jlo:jhi) - first_stage%rho%data(ilo:ihi, jlo:jhi)))
    rho_u_diff = maxval(abs(last_stage%rho_u%data(ilo:ihi, jlo:jhi) - first_stage%rho_u%data(ilo:ihi, jlo:jhi)))
    rho_v_diff = maxval(abs(last_stage%rho_v%data(ilo:ihi, jlo:jhi) - first_stage%rho_v%data(ilo:ihi, jlo:jhi)))
    rho_E_diff = maxval(abs(last_stage%rho_E%data(ilo:ihi, jlo:jhi) - first_stage%rho_E%data(ilo:ihi, jlo:jhi)))

    open(newunit=io, file=trim(first_stage%residual_hist_file), status='old', position="append")
    write(io, '(i0, ",", 5(es16.6, ","))') first_stage%iteration, first_stage%time * t_0, &
      rho_diff, rho_u_diff, rho_v_diff, rho_E_diff
    close(io)

  end subroutine write_residual_history

  subroutine apply_bc(self)
    !< Apply the boundary conditions
    class(fluid_t), intent(inout) :: self
    integer(ik) :: priority
    integer(ik) :: max_priority_bc !< highest goes first

    if(enable_debug_print) call debug_print('Running flux_solver_t%apply()', __FILE__, __LINE__)

    max_priority_bc = max(self%bc_plus_x%priority, self%bc_plus_y%priority, &
                          self%bc_minus_x%priority, self%bc_minus_y%priority)

    call self%bc_plus_x%set_time(time=self%time)
    call self%bc_minus_x%set_time(time=self%time)
    call self%bc_plus_y%set_time(time=self%time)
    call self%bc_minus_y%set_time(time=self%time)

    do priority = max_priority_bc, 0, -1

      if(self%bc_plus_x%priority == priority) then
        call self%bc_plus_x%apply(rho=self%rho, u=self%u, v=self%v, p=self%p)
      end if

      if(self%bc_plus_y%priority == priority) then
        call self%bc_plus_y%apply(rho=self%rho, u=self%u, v=self%v, p=self%p)
      end if

      if(self%bc_minus_x%priority == priority) then
        call self%bc_minus_x%apply(rho=self%rho, u=self%u, v=self%v, p=self%p)
      end if

      if(self%bc_minus_y%priority == priority) then
        call self%bc_minus_y%apply(rho=self%rho, u=self%u, v=self%v, p=self%p)
      end if

    end do

  end subroutine apply_bc

end module mod_fluid
