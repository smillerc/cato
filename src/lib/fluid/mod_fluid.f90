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
  !< Summary: Provide
  !< Date: 08/24/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !<      [1] D. Rouson, J. Xia, and X. Xu, "Scientific Software Design: The Object-Oriented Way"
  !<      [2] J. Blazek, "Computational Fluid Dynamics: Principles and Applications"
  !<      [3] S. Ruuth, R. Spiteri, "High-Order Strong-Stability-Preserving Runge–Kutta Methods with Downwind-Biased Spatial Discretizations",
  !<          SIAM Journal of Numerical Analysis, Vol. 42, No. 3, pp. 974–996, https://doi.org/10.1137/S0036142902419284

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
  use collectives, only: min_to_all, max_to_all
  ! use mod_ausm_plus_solver, only: ausm_plus_solver_t
  ! use mod_fvleg_solver, only: fvleg_solver_t
  use mod_m_ausmpw_plus_solver, only: m_ausmpw_plus_solver_t
  use mod_ausmpw_plus_solver, only: ausmpw_plus_solver_t
  ! use mod_slau_solver, only: slau_solver_t

  use mod_source, only: source_t

  implicit none

  private
  public :: fluid_t, new_fluid

  logical, parameter :: filter_small_mach = .false.

  type :: fluid_t
    !< Fluid solver physics package

    private ! make all private by default

    type(input_t), public :: input
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

    character(len=32) :: flux_solver_type = ''
    class(flux_solver_t), allocatable :: solver !< solver scheme used to flux quantities at cell interfaces

    ! Time variables
    character(len=10) :: time_integration_scheme = 'ssp_rk_2_2'
    real(rk) :: time = 0.0_rk !< current simulation time
    real(rk) :: dt = 0.0_rk   !< time step
    real(rk) :: cfl = 0.0_rk  !< Courant–Friedrichs–Lewy condition (CFL)
    integer(ik) :: iteration = 0 !< current iteration number

    logical, public :: prim_vars_updated = .false.

    ! Residual history
    character(len=32) :: residual_hist_file = 'residual_hist.csv'
    logical :: residual_hist_header_written = .false.

  contains
    ! Private methods
    private
    procedure :: initialize_from_hdf5
    procedure :: calculate_derived_quantities
    procedure :: sanity_check
    procedure :: ssp_rk_2_2
    procedure :: ssp_rk_3_3
    procedure :: ssp_rk_4_3
    procedure :: sync_fields

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

    ! Finalizer
    final :: finalize

    ! Map operators to corresponding procedures
    generic :: operator(+) => add_fluid
    generic :: operator(-) => subtract_fluid
    generic :: operator(*) => real_mul_fluid, fluid_mul_real
    generic :: assignment(=) => assign_fluid
  endtype fluid_t

contains

  function new_fluid(input, grid, time) result(fluid)
    !< Fluid constructor
    class(input_t), intent(in) :: input
    class(grid_block_t), intent(in) :: grid
    type(fluid_t), pointer :: fluid
    real(rk), intent(in) :: time

    select type(grid)
    class is(grid_block_2d_t)
      allocate(fluid)
      call fluid%initialize(input, grid, time)
    endselect
  endfunction new_fluid

  subroutine initialize(self, input, grid, time)
    class(fluid_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_block_2d_t), intent(in) :: grid
    real(rk), intent(in) :: time
    class(flux_solver_t), pointer :: solver => null()

    integer(ik) :: alloc_status, io

    self%input = input
    self%time = time
    self%cfl = input%cfl

    alloc_status = 0
    if(enable_debug_print) call debug_print('Calling fluid_t%initialize()', __FILE__, __LINE__)

    if(.not. scale_factors_set) then
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='initialize', &
                     message="Global non-dimensional scale factors haven't been set yet. "// &
                     "These need to be set before fluid initialization", &
                     file_name=__FILE__, line_number=__LINE__)
    endif

    self%rho = field_2d(name='rho', long_name='Density', &
                        descrip='Cell Density', units='g/cm^3', &
                        global_dims=grid%global_dims, n_halo_cells=input%n_ghost_layers)

    self%rho_u = field_2d(name='rhou', long_name='rhou', descrip='Cell Conserved quantity (Density * X-Velocity)', &
                          units='g cm/cm^2 s', &
                          global_dims=grid%global_dims, n_halo_cells=input%n_ghost_layers)

    self%rho_v = field_2d(name='rhov', long_name='rhov', descrip='Cell Conserved quantity (Density * Y-Velocity)', &
                          units='g cm/cm^2 s', &
                          global_dims=grid%global_dims, n_halo_cells=input%n_ghost_layers)

    self%rho_E = field_2d(name='rhoE', long_name='rhoE', descrip='Cell Conserved quantity (Density * Total Energy)', &
                          units='g erg / cm^3', &
                          global_dims=grid%global_dims, n_halo_cells=input%n_ghost_layers)

    self%u = field_2d(name='u', long_name='X Velocity', descrip='Cell X-Velocity', units='cm/s', &
                      global_dims=grid%global_dims, n_halo_cells=input%n_ghost_layers)

    self%v = field_2d(name='v', long_name='Y Velocity', descrip='Cell Y-Velocity', units='cm/s', &
                      global_dims=grid%global_dims, n_halo_cells=input%n_ghost_layers)

    self%p = field_2d(name='p', long_name='Pressure', descrip='Cell Pressure', units='barye', &
                      global_dims=grid%global_dims, n_halo_cells=input%n_ghost_layers)

    self%cs = field_2d(name='cs', long_name='Sound Speed', descrip='Cell Sound Speed', units='cm/s', &
                       global_dims=grid%global_dims, n_halo_cells=input%n_ghost_layers)

    self%mach_u = field_2d(name='mach_u', long_name='Mach X', descrip='Cell Mach number in x-direction', &
                           units='dimensionless', &
                           global_dims=grid%global_dims, n_halo_cells=input%n_ghost_layers)

    self%mach_v = field_2d(name='mach_v', long_name='Mach Y', descrip='Cell Mach number in y-direction', &
                           units='dimensionless', &
                           global_dims=grid%global_dims, n_halo_cells=input%n_ghost_layers)

    self%time_integration_scheme = trim(input%time_integration_strategy)
    self%flux_solver_type = trim(input%flux_solver)

    select case(trim(input%flux_solver))
    case('M-AUSMPW+')
      allocate(m_ausmpw_plus_solver_t :: solver)
    case('AUSMPW+')
      allocate(ausmpw_plus_solver_t :: solver)
    case default
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='initialize', &
                     message="Invalid flux solver. It must be one of the following: "// &
                     "['AUSMPW+', 'M-AUSMPW+'], "// &
                     "the input was: '"//trim(input%flux_solver)//"'", &
                     file_name=__FILE__, line_number=__LINE__)
    endselect

    call solver%initialize(input, self%time)
    allocate(self%solver, source=solver)
    deallocate(solver)

    open(newunit=io, file=trim(self%residual_hist_file), status='replace')
    write(io, '(a)') 'iteration,time,rho,rho_u,rho_v,rho_E'
    close(io)

    call self%initialize_from_hdf5(input)
    call eos%sound_speed(p=self%p, rho=self%rho, cs=self%cs)
    self%mach_v = self%v / self%cs
    self%mach_u = self%u / self%cs
    self%prim_vars_updated = .true.
    call self%sync_fields()
  endsubroutine initialize

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
    endif

    file_exists = .false.
    inquire(file=filename, exist=file_exists)

    if(.not. file_exists) then
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='initialize_from_hdf5', &
                     message='File not found: "'//filename//'"', file_name=__FILE__, line_number=__LINE__)
    endif

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
      endselect
    endif

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
      endselect
    endif

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
      endselect
    endif

    call h5%finalize()

    associate(ilo => self%rho%lbounds(1), ihi => self%rho%ubounds(1), &
              jlo => self%rho%lbounds(2), jhi => self%rho%ubounds(2), nh => self%rho%n_halo_cells)
      self%rho%data(ilo:ihi, jlo:jhi) = density(ilo + nh:ihi + nh, jlo + nh:jhi + nh)
      self%u%data(ilo:ihi, jlo:jhi) = x_velocity(ilo + nh:ihi + nh, jlo + nh:jhi + nh)
      self%v%data(ilo:ihi, jlo:jhi) = y_velocity(ilo + nh:ihi + nh, jlo + nh:jhi + nh)
      self%p%data(ilo:ihi, jlo:jhi) = pressure(ilo + nh:ihi + nh, jlo + nh:jhi + nh)
    endassociate

    call self%rho%make_non_dimensional(non_dim_factor=rho_0)
    call self%u%make_non_dimensional(non_dim_factor=v_0)
    call self%v%make_non_dimensional(non_dim_factor=v_0)
    call self%p%make_non_dimensional(non_dim_factor=p_0)

    if(minval(pressure) < tiny(1.0_rk)) then
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='initialize_from_hdf5', &
                     message="Some (or all) of the pressure array is ~0", file_name=__FILE__, line_number=__LINE__)
    endif

    if(minval(density) < tiny(1.0_rk)) then
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='initialize_from_hdf5', &
                     message="Some (or all) of the density array is ~0", file_name=__FILE__, line_number=__LINE__)
    endif

    call eos%primitive_to_conserved(rho=self%rho, u=self%u, v=self%v, p=self%p, &
                                    rho_u=self%rho_u, rho_v=self%rho_v, rho_E=self%rho_E)

    if(this_image() == 1) then
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
    endif

  endsubroutine initialize_from_hdf5

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
    if(allocated(self%solver)) deallocate(self%solver)
  endsubroutine finalize

  subroutine set_time(self, time, iteration)
    !< Set the time statistics
    class(fluid_t), intent(inout) :: self
    real(rk), intent(in) :: time          !< simulation time
    integer(ik), intent(in) :: iteration  !< iteration

    if(enable_debug_print) call debug_print('Running fluid_t%set_time()', __FILE__, __LINE__)
    self%time = time
    self%iteration = iteration
  endsubroutine set_time

  subroutine integrate(self, dt, grid, source_term, error_code)
    !< Integrate in time
    class(fluid_t), intent(inout) :: self
    real(rk), intent(in) :: dt !< time step
    class(grid_block_t), intent(in) :: grid !< grid class - the solver needs grid topology
    class(source_t), allocatable, intent(in) :: source_term

    integer(ik), intent(out) :: error_code

    if(enable_debug_print) call debug_print('Running fluid_t%integerate()', __FILE__, __LINE__)
    self%time = self%time + dt
    self%dt = dt

    select case(trim(self%time_integration_scheme))
    case('ssp_rk2')
      call self%ssp_rk_2_2(grid=grid, source_term=source_term, error_code=error_code)
    case('ssp_rk3', 'ssp_rk33')
      call self%ssp_rk_3_3(grid=grid, source_term=source_term, error_code=error_code)
    case('ssp_rk43')
      call self%ssp_rk_4_3(grid=grid, source_term=source_term, error_code=error_code)
    case default
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='assign_fluid', &
                     message="Unknown time integration scheme", file_name=__FILE__, line_number=__LINE__)
    endselect
  endsubroutine integrate

  subroutine sync_fields(self)
    class(fluid_t), intent(inout) :: self
    call self%rho%sync_edges()
    call self%u%sync_edges()
    call self%v%sync_edges()
    call self%p%sync_edges()
  endsubroutine sync_fields

  function time_derivative(self, grid, source_term) result(d_dt)
    !< Implementation of the time derivative
    class(fluid_t), intent(inout) :: self
    class(grid_block_t), intent(in) :: grid  !< grid class - the solver needs grid topology
    type(fluid_t), allocatable :: d_dt          !< dU/dt
    class(source_t), allocatable, intent(in) :: source_term

    real(rk), dimension(:, :), allocatable ::   d_rho_dt  !< d/dt of the density field
    real(rk), dimension(:, :), allocatable :: d_rho_u_dt  !< d/dt of the rhou field
    real(rk), dimension(:, :), allocatable :: d_rho_v_dt  !< d/dt of the rhov field
    real(rk), dimension(:, :), allocatable :: d_rho_E_dt  !< d/dt of the rhoE field

    if(enable_debug_print) call debug_print('Running fluid_t%time_derivative()', __FILE__, __LINE__)

    allocate(d_dt, source=self)

    select type(grid)
    class is(grid_block_2d_t)
      call self%sync_fields()
      call self%solver%solve(dt=self%dt, grid=grid, &
                             rho=self%rho, u=self%u, v=self%v, p=self%p, &
                             d_rho_dt=d_rho_dt, &
                             d_rho_u_dt=d_rho_u_dt, &
                             d_rho_v_dt=d_rho_v_dt, &
                             d_rho_E_dt=d_rho_E_dt)

      d_dt%rho%data = d_rho_dt
      d_dt%rho_u%data = d_rho_u_dt
      d_dt%rho_v%data = d_rho_v_dt
      d_dt%rho_E%data = d_rho_E_dt

      deallocate(d_rho_dt)
      deallocate(d_rho_u_dt)
      deallocate(d_rho_v_dt)
      deallocate(d_rho_E_dt)

      if(allocated(source_term)) then
        if(self%time <= source_term%max_time) then
          if(source_term%source_type == 'energy') then
            d_dt%rho_E = source_term%data + d_dt%rho_E
          endif
        endif
      endif
    endselect
  endfunction time_derivative

  real(rk) function get_timestep(self, grid) result(delta_t)
    class(fluid_t), intent(inout) :: self
    class(grid_block_t), intent(in) :: grid
    real(rk) :: min_delta_t !< the max allowable timestep on each subdomain/image
    integer(ik) :: ilo, ihi, jlo, jhi ! fluid lo/hi indicies
    integer(ik) :: g_ilo, g_ihi, g_jlo, g_jhi ! grid lo/hi indices
    character(len=200) :: err_msg

    real(rk), allocatable, save, dimension(:, :) :: dx, dy

    ! This seems silly to have to do, but the master class requires
    ! that the grid is a grid_block_t type, even though this fluid class
    ! will always use a grid_block_2d_t type.
    select type(grid)
    class is(grid_block_2d_t)
      g_ilo = grid%lbounds(1)
      g_ihi = grid%ubounds(1)
      g_jlo = grid%lbounds(2)
      g_jhi = grid%ubounds(2)

      if(.not. allocated(dx)) allocate(dx(g_ilo:g_ihi, g_jlo:g_jhi))
      if(.not. allocated(dy)) allocate(dy(g_ilo:g_ihi, g_jlo:g_jhi))

      dx = grid%dx(g_ilo:g_ihi, g_jlo:g_jhi)
      dy = grid%dy(g_ilo:g_ihi, g_jlo:g_jhi)
    endselect

    err_msg = ''

    ilo = self%u%lbounds(1)
    ihi = self%u%ubounds(1)

    jlo = self%u%lbounds(2)
    jhi = self%u%ubounds(2)

    ! I would have put this in a cleaner associate block, but GFortran+OpenCoarrays bugs out on this
    min_delta_t = minval(self%cfl / &
                         (((abs(self%u%data(ilo:ihi, jlo:jhi)) + &
                            self%cs%data(ilo:ihi, jlo:jhi)) / dx) + &
                          ((abs(self%v%data(ilo:ihi, jlo:jhi)) + &
                            self%cs%data(ilo:ihi, jlo:jhi)) / dy)))

    min_delta_t = min_to_all(min_delta_t)
    delta_t = min_delta_t
  endfunction get_timestep

  subroutine calculate_derived_quantities(self)
    !< Find derived quantities like sound speed, mach number, primitive variables

    class(fluid_t), intent(inout) :: self

    if(enable_debug_print) call debug_print('Running fluid_t%calculate_derived_quantities()', __FILE__, __LINE__)
    call eos%conserved_to_primitive(rho=self%rho, rho_u=self%rho_u, &
                                    rho_v=self%rho_v, rho_E=self%rho_E, &
                                    u=self%u, v=self%v, p=self%p)
    call self%rho%check_for_negatives()
    call self%p%check_for_negatives()

    call eos%sound_speed(p=self%p, rho=self%rho, cs=self%cs)

    self%mach_u = self%u / self%cs
    self%mach_v = self%v / self%cs

    ! Filter out the super small Mach numbers
    where(abs(self%mach_v%data) < 1e-13_rk)
      self%v%data = 0.0_rk
      self%mach_v%data = 0.0_rk
    endwhere

    where(abs(self%mach_u%data) < 1e-13_rk)
      self%u%data = 0.0_rk
      self%mach_u%data = 0.0_rk
    endwhere

    self%prim_vars_updated = .true.
  endsubroutine calculate_derived_quantities

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
  endfunction subtract_fluid

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
  endfunction add_fluid

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
    endif

    product%rho = lhs%rho * rhs
    product%rho_u = lhs%rho_u * rhs
    product%rho_v = lhs%rho_v * rhs
    product%rho_E = lhs%rho_E * rhs

    product%prim_vars_updated = .false.
  endfunction fluid_mul_real

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
    endif

    product%rho = lhs * rhs%rho
    product%rho_u = lhs * rhs%rho_u
    product%rho_v = lhs * rhs%rho_v
    product%rho_E = lhs * rhs%rho_E
    product%prim_vars_updated = .false.
  endfunction real_mul_fluid

  subroutine assign_fluid(lhs, rhs)
    !< Implementation of the (=) operator for the fluid type. e.g. lhs = rhs
    class(fluid_t), intent(inout) :: lhs
    class(fluid_t), intent(in) :: rhs
    character(len=50) :: err_msg
    class(flux_solver_t), pointer :: solver => null()
    err_msg = ''

    if(enable_debug_print) call debug_print('Running fluid_t%assign_fluid()', __FILE__, __LINE__)

    lhs%input = rhs%input
    lhs%residual_hist_header_written = rhs%residual_hist_header_written
    lhs%rho = rhs%rho
    lhs%rho_u = rhs%rho_u
    lhs%rho_v = rhs%rho_v
    lhs%rho_E = rhs%rho_E

    select case(trim(rhs%flux_solver_type))
    case('M-AUSMPW+')
      allocate(m_ausmpw_plus_solver_t :: solver)
    case('AUSMPW+')
      allocate(ausmpw_plus_solver_t :: solver)
    case default
      call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='assign_fluid', &
                     message="Invalid flux solver. It must be one of the following: "// &
                     "['AUSMPW+', 'M-AUSMPW+'], "// &
                     "the input was: '"//trim(rhs%flux_solver_type)//"'", &
                     file_name=__FILE__, line_number=__LINE__)
    endselect

    if(allocated(lhs%solver)) deallocate(lhs%solver)
    call solver%initialize(rhs%input, rhs%time)
    allocate(lhs%solver, source=solver)
    deallocate(solver)

    lhs%time_integration_scheme = rhs%time_integration_scheme
    lhs%time = rhs%time
    lhs%dt = rhs%dt
    lhs%iteration = rhs%iteration
    lhs%prim_vars_updated = rhs%prim_vars_updated
    lhs%residual_hist_file = rhs%residual_hist_file
    lhs%residual_hist_header_written = rhs%residual_hist_header_written

    call lhs%rho%check_for_nans()
    call lhs%u%check_for_nans()
    call lhs%v%check_for_nans()
    call lhs%p%check_for_nans()

    call lhs%calculate_derived_quantities()
  endsubroutine assign_fluid

  subroutine sanity_check(self, error_code)
    !< Run checks on the conserved variables. Density and pressure need to be > 0. No NaNs or Inifinite numbers either.
    class(fluid_t), intent(in) :: self
    integer(ik), intent(out) :: error_code

    error_code = 0

    if(enable_debug_print) call debug_print('Running fluid_t%sanity_check()', __FILE__, __LINE__)
    call self%rho%check_for_negatives()
    call self%p%check_for_negatives()

  endsubroutine sanity_check

  subroutine ssp_rk_3_3(U, source_term, grid, error_code)
    !< Strong-stability preserving Runge-Kutta 3-step, 3rd order time integration. See Ref [3]
    class(fluid_t), intent(inout) :: U
    class(grid_block_t), intent(in) :: grid
    class(source_t), allocatable, intent(in) :: source_term

    type(fluid_t), allocatable :: U_1 !< first stage
    type(fluid_t), allocatable :: U_2 !< second stage
    integer(ik), intent(out) :: error_code
    real(rk) :: dt

    call debug_print('Running fluid_t%ssp_rk_3_3()', __FILE__, __LINE__)

    dt = U%dt

    ! 1st stage
    allocate(U_1, source=U)
    U_1 = U + U%t(grid=grid, source_term=source_term) * dt

    ! 2nd stage
    allocate(U_2, source=U)
    U_2 = U * (3.0_rk / 4.0_rk) &
          + U_1 * (1.0_rk / 4.0_rk) &
          + U_1%t(grid=grid, source_term=source_term) * ((1.0_rk / 4.0_rk) * dt)

    ! Final stage
    U = U * (1.0_rk / 3.0_rk) &
        + U_2 * (2.0_rk / 3.0_rk) &
        + U_2%t(grid=grid, source_term=source_term) * ((2.0_rk / 3.0_rk) * dt)
    call U%sanity_check(error_code)
    call U%calculate_derived_quantities()

    ! Convergence history
    call write_residual_history(first_stage=U_1, last_stage=U)

    deallocate(U_1)
    deallocate(U_2)

  endsubroutine ssp_rk_3_3

  subroutine ssp_rk_4_3(U, source_term, grid, error_code)
    !< Strong-stability preserving Runge-Kutta 4-step, 3rd order time integration. See Ref [3]. According to the
    !< reference, the increase in stage number is more than offset by the allowable increase in CFL number.

    class(fluid_t), intent(inout) :: U
    class(grid_block_t), intent(in) :: grid
    class(source_t), allocatable, intent(in) :: source_term

    type(fluid_t), allocatable :: U_1 !< first stage
    type(fluid_t), allocatable :: U_2 !< second stage
    type(fluid_t), allocatable :: U_3 !< third stage
    integer(ik), intent(out) :: error_code
    real(rk) :: dt
    real(rk), parameter :: one_third = 1.0_rk / 3.0_rk
    real(rk), parameter :: one_sixth = 1.0_rk / 6.0_rk
    real(rk), parameter :: two_thirds = 2.0_rk / 3.0_rk

    call debug_print('Running fluid_t%ssp_rk_4_3()', __FILE__, __LINE__)

    dt = U%dt

    ! 1st stage
    allocate(U_1, source=U)
    U_1 = U + U%t(grid=grid, source_term=source_term) * 0.5_rk * dt

    ! 2nd stage
    allocate(U_2, source=U)
    U_2 = U_1 + U_1%t(grid=grid, source_term=source_term) * 0.5_rk * dt

    ! 3rd stage
    allocate(U_3, source=U)
    U_3 = (U * two_thirds) + (U_2 * one_third) + (U_2%t(grid=grid, source_term=source_term) * one_sixth * dt)

    ! Final stage
    U = U_3 + (U_3%t(grid=grid, source_term=source_term) * 0.5_rk * dt)
    call U%sanity_check(error_code)
    call U%calculate_derived_quantities()

    ! Convergence history
    call write_residual_history(first_stage=U_1, last_stage=U)

    deallocate(U_1)
    deallocate(U_2)
    deallocate(U_3)

  endsubroutine ssp_rk_4_3

  subroutine ssp_rk_2_2(U, source_term, grid, error_code)
    !< Strong-stability preserving Runge-Kutta 2-stage, 2nd order time integration. See Ref [3]
    class(fluid_t), intent(inout) :: U
    class(grid_block_t), intent(in) :: grid
    class(source_t), allocatable, intent(in) :: source_term

    type(fluid_t), allocatable :: U_1 !< first stage
    integer(ik), intent(out) :: error_code
    real(rk) :: dt

    if(enable_debug_print) call debug_print('Running fluid_t%rk2()', __FILE__, __LINE__)

    dt = U%dt
    allocate(U_1, source=U)

    ! 1st stage
    U_1 = U + U%t(grid=grid, source_term=source_term) * dt

    ! Final stage
    U = U * 0.5_rk + U_1 * 0.5_rk + &
        U_1%t(grid=grid, source_term=source_term) * (0.5_rk * dt)
    call U%sanity_check(error_code)
    call U%calculate_derived_quantities()

    ! ! Convergence history
    call write_residual_history(first_stage=U_1, last_stage=U)

    deallocate(U_1)

  endsubroutine ssp_rk_2_2

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
    rho_diff = max_to_all(rho_diff)
    rho_u_diff = maxval(abs(last_stage%rho_u%data(ilo:ihi, jlo:jhi) - first_stage%rho_u%data(ilo:ihi, jlo:jhi)))
    rho_u_diff = max_to_all(rho_u_diff)
    rho_v_diff = maxval(abs(last_stage%rho_v%data(ilo:ihi, jlo:jhi) - first_stage%rho_v%data(ilo:ihi, jlo:jhi)))
    rho_v_diff = max_to_all(rho_v_diff)
    rho_E_diff = maxval(abs(last_stage%rho_E%data(ilo:ihi, jlo:jhi) - first_stage%rho_E%data(ilo:ihi, jlo:jhi)))
    rho_E_diff = max_to_all(rho_E_diff)

    open(newunit=io, file=trim(first_stage%residual_hist_file), status='old', position="append")
    write(io, '(i0, ",", 5(es16.6, ","))') first_stage%iteration, first_stage%time * t_0, &
      rho_diff, rho_u_diff, rho_v_diff, rho_E_diff
    close(io)

  endsubroutine write_residual_history
endmodule mod_fluid
