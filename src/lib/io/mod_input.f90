module mod_input

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, output_unit
  use cfgio_mod, only: cfg_t, parse_cfg

  implicit none

  private
  public :: input_t

  type :: input_t
    ! general
    character(:), allocatable :: title  !< Name of the simulation

    ! grid
    character(len=32) :: grid_type = '2d_regular'  !< Structure/layout of the grid, e.g. '2d_regular'
    real(rk) :: xmin = 0.0_rk   !< Minimum extent of the grid in x (ignored for .h5 initial grids)
    real(rk) :: xmax = 0.0_rk   !< Maximum extent of the grid in x (ignored for .h5 initial grids)
    real(rk) :: ymin = 0.0_rk   !< Minimum extent of the grid in y (ignored for .h5 initial grids)
    real(rk) :: ymax = 0.0_rk   !< Maximum extent of the grid in y (ignored for .h5 initial grids)
    integer(ik) :: ni_nodes = 0 !< # of i nodes (not including ghost) (ignored for .h5 initial grids)
    integer(ik) :: nj_nodes = 0 !< # of j nodes (not including ghost) (ignored for .h5 initial grids)

    ! initial conditions
    character(:), allocatable :: initial_condition_file
    logical :: read_init_cond_from_file = .false.
    real(rk) :: init_x_velocity = 0.0_rk !< Initial condition for the x velocity field (useful only for testing)
    real(rk) :: init_y_velocity = 0.0_rk !< Initial condition for the y velocity field (useful only for testing)
    real(rk) :: init_density = 0.0_rk    !< Initial condition for the density field (useful only for testing)
    real(rk) :: init_pressure = 0.0_rk   !< Initial condition for the pressure field (useful only for testing)

    ! boundary conditions
    character(len=32) ::  plus_x_bc = 'periodic' !< Boundary condition at +x
    character(len=32) :: minus_x_bc = 'periodic' !< Boundary condition at -x
    character(len=32) ::  plus_y_bc = 'periodic' !< Boundary condition at +y
    character(len=32) :: minus_y_bc = 'periodic' !< Boundary condition at -y

    ! io
    character(:), allocatable :: contour_io_format !< e.g. 'xdmf'

    ! timing
    real(rk) :: max_time = 1.0_rk
    real(rk) :: initial_delta_t = 0.0_rk
    real(rk) :: contour_interval_dt = 0.5_rk
    character(:), allocatable :: time_integration_strategy !< How is time integration handled? e.g. 'rk2', 'rk4', etc.

    ! physics
    real(rk) :: polytropic_index = 5.0_rk / 3.0_rk !< e.g. gamma for the simulated gas

    ! finite volume scheme specifics
    character(len=32) :: reconstruction_type !< How are the cells being reconstructed
    real(rk) :: tau = 1.0e-5_rk !< time increment for FVEG and FVLEG schemes
    character(:), allocatable :: slope_limiter

  contains
    procedure, public :: initialize
    procedure, public :: read_from_ini
  end type input_t

contains

  subroutine initialize(self, ni, nj, xmin, xmax, ymin, ymax)
    class(input_t), intent(inout) :: self
    integer(ik), intent(in) :: ni, nj
    real(rk), intent(in) :: xmin, xmax, ymin, ymax

    self%ni_nodes = ni
    self%nj_nodes = nj
    self%xmin = xmin
    self%xmax = xmax
    self%ymin = ymin
    self%ymax = ymax
    self%reconstruction_type = 'piecewise_linear'
    self%slope_limiter = 'sun_ren_09'

    self%plus_x_bc = 'periodic'
    self%minus_x_bc = 'periodic'
    self%plus_y_bc = 'periodic'
    self%minus_y_bc = 'periodic'

  end subroutine initialize

  subroutine read_from_ini(self, filename)

    character(len=*), intent(in) :: filename
    class(input_t), intent(inout) :: self

    character(len=120) :: char_buffer
    type(cfg_t) :: cfg
    logical :: file_exists

    file_exists = .false.

    ! if(this_image() == 1) then
    write(output_unit, '(a)') 'Reading input file: '//trim(filename)
    ! end if

    inquire(file=filename, EXIST=file_exists)
    if(.not. file_exists) error stop "Input .ini file not found"

    cfg = parse_cfg(trim(filename))

    ! General
    call cfg%get("general", "title", char_buffer)
    self%title = trim(char_buffer)

    ! Time
    call cfg%get("time", "max_time", self%max_time)
    call cfg%get("time", "initial_delta_t", self%initial_delta_t)
    call cfg%get("time", "integration_strategy", char_buffer)
    self%time_integration_strategy = trim(char_buffer)

    ! Physics
    call cfg%get("physics", "polytropic_index", self%polytropic_index)

    ! Scheme
    call cfg%get("scheme", "tau", self%tau)

    call cfg%get("scheme", "reconstruction_type", char_buffer, 'piecewise_linear')
    self%reconstruction_type = trim(char_buffer)

    call cfg%get("scheme", "slope_limiter", char_buffer, 'sun_ren_09')
    self%slope_limiter = trim(char_buffer)

    ! Initial Conditions
    call cfg%get("initial_conditions", "read_from_file", self%read_init_cond_from_file, .true.)
    if(self%read_init_cond_from_file) then
      call cfg%get("initial_conditions", "initial_condition_file", char_buffer)
      self%initial_condition_file = trim(char_buffer)

    end if

    call cfg%get("initial_conditions", "init_density", self%init_density, 0.0_rk)
    call cfg%get("initial_conditions", "init_x_velocity", self%init_x_velocity, 0.0_rk)
    call cfg%get("initial_conditions", "init_y_velocity", self%init_y_velocity, 0.0_rk)
    call cfg%get("initial_conditions", "init_pressure", self%init_pressure, 0.0_rk)

    ! Grid
    call cfg%get("grid", "grid_type", char_buffer, '2d_regular')
    self%grid_type = trim(char_buffer)

    if(.not. self%read_init_cond_from_file) then
      call cfg%get("grid", "ni_nodes", self%ni_nodes)
      call cfg%get("grid", "xmin", self%xmin)
      call cfg%get("grid", "xmax", self%xmax)

      call cfg%get("grid", "nj_nodes", self%nj_nodes)
      call cfg%get("grid", "ymin", self%ymin)
      call cfg%get("grid", "ymax", self%ymax)
    end if

    ! Boundary conditions
    call cfg%get("boundary_conditions", "plus_x", char_buffer, 'periodic')
    self%plus_x_bc = trim(char_buffer)
    call cfg%get("boundary_conditions", "minus_x", char_buffer, 'periodic')
    self%minus_x_bc = trim(char_buffer)

    call cfg%get("boundary_conditions", "plus_y", char_buffer, 'periodic')
    self%plus_y_bc = trim(char_buffer)
    call cfg%get("boundary_conditions", "minus_y", char_buffer, 'periodic')
    self%minus_y_bc = trim(char_buffer)

    ! Input/Output
    call cfg%get("io", "format", char_buffer, 'xdmf')
    self%contour_io_format = trim(char_buffer)

    call cfg%get("io", "contour_interval_dt", self%contour_interval_dt, 0.1_rk)

  end subroutine read_from_ini

end module mod_input
