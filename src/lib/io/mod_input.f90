module mod_input

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, output_unit
  use cfgio_mod, only: cfg_t, parse_cfg

  implicit none

  private
  public :: input_t

  type :: input_t
    ! general
    character(:), allocatable :: title  !< Name of the simulation
    character(:), allocatable :: unit_system !< Type of units to output (CGS, ICF, MKS, etc...)

    ! reference state
    real(rk) :: reference_pressure = 1.0_rk
    real(rk) :: reference_density = 1.0_rk

    ! grid
    character(len=32) :: grid_type = '2d_regular'  !< Structure/layout of the grid, e.g. '2d_regular'
    real(rk) :: xmin = 0.0_rk   !< Minimum extent of the grid in x (ignored for .h5 initial grids)
    real(rk) :: xmax = 0.0_rk   !< Maximum extent of the grid in x (ignored for .h5 initial grids)
    real(rk) :: ymin = 0.0_rk   !< Minimum extent of the grid in y (ignored for .h5 initial grids)
    real(rk) :: ymax = 0.0_rk   !< Maximum extent of the grid in y (ignored for .h5 initial grids)
    integer(ik) :: ni_nodes = 0 !< # of i nodes (not including ghost) (ignored for .h5 initial grids)
    integer(ik) :: nj_nodes = 0 !< # of j nodes (not including ghost) (ignored for .h5 initial grids)
    integer(ik) :: n_ghost_layers = 1 !< # of ghost layers to use (based on spatial reconstruction order)

    ! initial conditions
    character(:), allocatable :: initial_condition_file
    logical :: read_init_cond_from_file = .false.
    real(rk) :: init_x_velocity = 0.0_rk !< Initial condition for the x velocity field (useful only for testing)
    real(rk) :: init_y_velocity = 0.0_rk !< Initial condition for the y velocity field (useful only for testing)
    real(rk) :: init_density = 0.0_rk    !< Initial condition for the density field (useful only for testing)
    real(rk) :: init_pressure = 0.0_rk   !< Initial condition for the pressure field (useful only for testing)

    ! restart files
    logical :: restart_from_file = .false.
    character(:), allocatable :: restart_file

    ! boundary conditions
    character(:), allocatable :: bc_pressure_input_file
    logical :: apply_constant_bc_pressure = .false.
    real(rk) :: constant_bc_pressure_value = 0.0_rk
    real(rk) :: bc_density = 5e-3_rk !< density input for a pressure boundary condition
    real(rk) :: bc_pressure_scale_factor = 1.0_rk
    character(len=32) ::  plus_x_bc = 'periodic' !< Boundary condition at +x
    character(len=32) :: minus_x_bc = 'periodic' !< Boundary condition at -x
    character(len=32) ::  plus_y_bc = 'periodic' !< Boundary condition at +y
    character(len=32) :: minus_y_bc = 'periodic' !< Boundary condition at -y

    ! source inputs (i.e. energy source injection at a grid location)
    logical :: enable_source_terms = .false.
    character(len=32) :: source_term_type = 'energy'
    logical :: apply_constant_source = .false.
    character(:), allocatable :: source_file
    real(rk) :: constant_source_value = 0.0_rk
    real(rk) :: source_scale_factor = 1.0_rk
    integer(ik) :: source_ilo = 0
    integer(ik) :: source_ihi = 0
    integer(ik) :: source_jlo = 0
    integer(ik) :: source_jhi = 0

    ! io
    character(:), allocatable :: contour_io_format !< e.g. 'xdmf'
    logical :: append_date_to_result_folder = .false.
    logical :: plot_reconstruction_states = .false.
    logical :: plot_reference_states = .false.
    logical :: plot_evolved_states = .false.
    logical :: plot_64bit = .false.
    logical :: plot_ghost_cells = .false.

    ! timing
    real(rk) :: max_time = 1.0_rk
    real(rk) :: cfl = 0.1_rk ! Courant–Friedrichs–Lewy condition
    real(rk) :: initial_delta_t = 0.0_rk
    real(rk) :: contour_interval_dt = 0.5_rk
    integer(ik) :: max_iterations = huge(1)
    character(:), allocatable :: time_integration_strategy !< How is time integration handled? e.g. 'rk2', 'rk4', etc.
    logical :: smooth_residuals = .true.

    ! physics
    real(rk) :: polytropic_index = 5.0_rk / 3.0_rk !< e.g. gamma for the simulated gas

    ! finite volume scheme specifics
    character(len=32) :: evolution_operator_type = 'fvleg'        !< How are the cells being reconstructed
    character(len=32) :: cell_reconstruction = 'piecewise_linear' !< How are the cells being reconstructed
    character(len=32) :: gradient_scheme = 'green_gauss'          !< How is the gradient estimated?
    character(len=32) :: edge_interpolation_scheme = 'TVD2'       !< How are the edge values interpolated?

    real(rk) :: tau = 1.0e-5_rk !< time increment for FVEG and FVLEG schemes
    character(:), allocatable :: limiter

    ! debug
    logical :: plot_limiters = .false.
    logical :: plot_gradients = .false.

  contains
    procedure, public :: initialize
    procedure, public :: read_from_ini
    procedure, public :: display_config
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
    self%cell_reconstruction = 'piecewise_linear'
    self%limiter = 'minmod'

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
    integer(ik) :: required_n_ghost_layers = 0

    file_exists = .false.

    ! if(this_image() == 1) then
    write(output_unit, '(a)') 'Reading input file: '//trim(filename)
    ! end if

    inquire(file=filename, exist=file_exists)
    if(.not. file_exists) error stop "Input .ini file not found"

    cfg = parse_cfg(trim(filename))

    ! General
    call cfg%get("general", "title", char_buffer)
    self%title = trim(char_buffer)

    call cfg%get("general", "units", char_buffer, 'cgs')
    self%unit_system = trim(char_buffer)

    ! Reference state
    call cfg%get("reference_state", "reference_pressure", self%reference_pressure)
    call cfg%get("reference_state", "reference_density", self%reference_density)

    ! Time
    call cfg%get("time", "max_time", self%max_time)
    call cfg%get("time", "initial_delta_t", self%initial_delta_t, 0.0_rk)
    call cfg%get("time", "cfl", self%cfl, 0.1_rk)
    call cfg%get("time", "integration_strategy", char_buffer)
    call cfg%get("time", "max_iterations", self%max_iterations, huge(1))
    self%time_integration_strategy = trim(char_buffer)

    ! Physics
    call cfg%get("physics", "polytropic_index", self%polytropic_index)

    ! Scheme
    call cfg%get("scheme", "smooth_residuals", self%smooth_residuals, .true.)
    call cfg%get("scheme", "cell_reconstruction", char_buffer, 'piecewise_linear')
    self%cell_reconstruction = trim(char_buffer)

    call cfg%get("scheme", "gradient_scheme", char_buffer, 'green_gauss')
    self%gradient_scheme = trim(char_buffer)

    call cfg%get("scheme", "edge_interpolation_scheme", char_buffer, 'TVD2')
    self%edge_interpolation_scheme = trim(char_buffer)

    call cfg%get("scheme", "limiter", char_buffer, 'minmod')
    self%limiter = trim(char_buffer)

    ! Restart files
    call cfg%get("restart", "restart_from_file", self%restart_from_file, .false.)

    write(*, '(a, l2)') 'Reading from restart?', self%restart_from_file

    if(self%restart_from_file) then
      call cfg%get("restart", "restart_file", char_buffer)
      self%restart_file = trim(char_buffer)

      inquire(file=self%restart_file, exist=file_exists)
      if(.not. file_exists) then
        write(*, '(a)') 'Restart file not found: "'//trim(self%restart_file)//'"'
        error stop "Restart file not found"
      end if
    end if

    ! Initial Conditions
    if(.not. self%restart_from_file) then
      call cfg%get("initial_conditions", "read_from_file", self%read_init_cond_from_file, .true.)

      if(self%read_init_cond_from_file) then
        call cfg%get("initial_conditions", "initial_condition_file", char_buffer)
        self%initial_condition_file = trim(char_buffer)

        inquire(file=self%initial_condition_file, exist=file_exists)
        if(.not. file_exists) then
          write(*, '(a)') 'Initial conditions file not found: "'//trim(self%initial_condition_file)//'"'
          error stop "Initial conditions file not found"
        end if

      end if

      call cfg%get("initial_conditions", "init_density", self%init_density, 0.0_rk)
      call cfg%get("initial_conditions", "init_x_velocity", self%init_x_velocity, 0.0_rk)
      call cfg%get("initial_conditions", "init_y_velocity", self%init_y_velocity, 0.0_rk)
      call cfg%get("initial_conditions", "init_pressure", self%init_pressure, 0.0_rk)
    end if

    ! Grid
    select case(trim(self%edge_interpolation_scheme))
    case('TVD2')
      required_n_ghost_layers = 1
      call cfg%get("grid", "n_ghost_layers", self%n_ghost_layers, required_n_ghost_layers)
    case('TVD3', 'TVD5', 'MLP3', 'MLP5')
      required_n_ghost_layers = 2
      call cfg%get("grid", "n_ghost_layers", self%n_ghost_layers, required_n_ghost_layers)
    case default
      error stop "Unknown edge interpolation scheme, must be one of the following: "// &
        "'TVD2', 'TVD3', 'TVD5', 'MLP3', or 'MLP5'"
    end select

    if(self%n_ghost_layers /= required_n_ghost_layers) then
      error stop "The number of required ghost cell layers doesn't match the edge interpolation order"
    end if

    call cfg%get("grid", "grid_type", char_buffer, '2d_regular')
    self%grid_type = trim(char_buffer)

    if(.not. self%read_init_cond_from_file .and. .not. self%restart_from_file) then
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

    if(self%plus_x_bc == 'pressure_input' .or. &
       self%minus_x_bc == 'pressure_input' .or. &
       self%plus_y_bc == 'pressure_input' .or. &
       self%minus_y_bc == 'pressure_input') then

      call cfg%get("boundary_conditions", "apply_constant_bc_pressure", self%apply_constant_bc_pressure, .false.)
      call cfg%get("boundary_conditions", "bc_density", self%bc_density)

      if(self%apply_constant_bc_pressure) then
        call cfg%get("boundary_conditions", "constant_bc_pressure_value", self%constant_bc_pressure_value)
      else
        call cfg%get("boundary_conditions", "bc_pressure_input_file", char_buffer)
        self%bc_pressure_input_file = trim(char_buffer)
      end if

      call cfg%get("boundary_conditions", "bc_pressure_scale_factor", self%bc_pressure_scale_factor, 1.0_rk)
    end if

    ! Source terms
    call cfg%get('source_terms', 'enable_source_terms', self%enable_source_terms, .false.)

    if(self%enable_source_terms) then
      call cfg%get('source_terms', 'apply_constant_source', self%apply_constant_source, .false.)
      call cfg%get('source_terms', 'source_file', char_buffer)
      self%source_file = trim(char_buffer)

      call cfg%get('source_terms', 'source_term_type', char_buffer)
      self%source_term_type = trim(char_buffer)

      call cfg%get('source_terms', 'constant_source_value', self%constant_source_value, 0.0_rk)
      call cfg%get('source_terms', 'source_scale_factor', self%source_scale_factor, 1.0_rk)
      call cfg%get('source_terms', 'ilo', self%source_ilo, 0)
      call cfg%get('source_terms', 'ihi', self%source_ihi, 0)
      call cfg%get('source_terms', 'jlo', self%source_jlo, 0)
      call cfg%get('source_terms', 'jhi', self%source_jhi, 0)

      if(self%source_ilo == 0 .and. &
         self%source_ihi == 0 .and. &
         self%source_jlo == 0 .and. &
         self%source_jhi == 0) then
        error stop "All of the (i,j) ranges in source_terms are 0"
      end if

      if(self%source_ilo > self%source_ihi .and. self%source_ihi /= -1) then
        error stop "Error: Invalid source application range; source_terms%ilo > source_terms%ihi"
      else if(self%source_jlo > self%source_jhi .and. self%source_jhi /= -1) then
        error stop "Error: Invalid source application range; source_terms%jlo > source_terms%jhi"
      end if
    end if

    ! Input/Output
    call cfg%get("io", "format", char_buffer, 'xdmf')
    self%contour_io_format = trim(char_buffer)

    call cfg%get("io", "contour_interval_dt", self%contour_interval_dt, 0.1_rk)
    call cfg%get("io", "append_date_to_result_folder", self%append_date_to_result_folder, .false.)
    call cfg%get("io", "plot_reconstruction_states", self%plot_reconstruction_states, .false.)
    call cfg%get("io", "plot_reference_states", self%plot_reference_states, .false.)
    call cfg%get("io", "plot_evolved_states", self%plot_evolved_states, .false.)
    call cfg%get("io", "plot_ghost_cells", self%plot_ghost_cells, .false.)
    call cfg%get("io", "plot_64bit", self%plot_64bit, .false.)

  end subroutine read_from_ini

  subroutine display_config(self)
    class(input_t), intent(in) :: self

    write(*, '(a)') "Input settings:"
    write(*, '(a)') "==============="
    write(*, *)
    write(*, '(a)') "[general]"
    write(*, '(3(a))') "title = '", trim(self%title), "'"
    write(*, '(3(a))') "unit_system = '", trim(self%unit_system), "'"

    write(*, *)
    write(*, '(a)') "[time]"
    write(*, '(a, es9.3)') "max_time = ", self%max_time
    write(*, '(a, f5.3)') "cfl = ", self%cfl
    write(*, '(a, es9.3)') "initial_delta_t = ", self%initial_delta_t
    write(*, '(a, i0)') "max_iterations = ", self%max_iterations
    write(*, '(3(a))') "time_integration_strategy = '", trim(self%unit_system), "'"

    write(*, *)
    write(*, '(a)') "[grid]"
    write(*, '(a, a)') "grid_type = ", self%grid_type
    write(*, '(a, es9.3)') "xmin = ", self%xmin
    write(*, '(a, es9.3)') "xmax = ", self%xmax
    write(*, '(a, es9.3)') "ymin = ", self%ymin
    write(*, '(a, es9.3)') "ymax = ", self%ymax
    write(*, '(a, i0)') "ni_nodes = ", self%ni_nodes
    write(*, '(a, i0)') "nj_nodes = ", self%nj_nodes
    write(*, '(a, i0)') "n_ghost_layers = ", self%n_ghost_layers

    write(*, *)
    write(*, '(a)') "[reference_state]"
    write(*, '(a, es9.3)') "reference_pressure = ", self%reference_pressure
    write(*, '(a, es9.3)') "reference_density = ", self%reference_density

    write(*, *)
    write(*, '(a)') "[initial_conditions]"
    write(*, '(a, a)') "initial_condition_file = ", self%initial_condition_file
    write(*, '(a, l1)') "read_init_cond_from_file = ", self%read_init_cond_from_file
    write(*, '(a, es9.3)') "init_x_velocity = ", self%init_x_velocity
    write(*, '(a, es9.3)') "init_y_velocity = ", self%init_y_velocity
    write(*, '(a, es9.3)') "init_density  = ", self%init_density
    write(*, '(a, es9.3)') "init_pressure = ", self%init_pressure

    write(*, *)
    write(*, '(a)') "[restart]"
    write(*, '(a, l1)') "restart_from_file = ", self%restart_from_file
    write(*, '(a, a)') "restart_file = ", self%restart_file

    write(*, *)
    write(*, '(a)') "[source_terms]"
    write(*, '(a, l1)') "enable_source_terms = ", self%enable_source_terms
    write(*, '(a, a)') "source_term_type = ", self%source_term_type
    write(*, '(a, l1)') "apply_constant_source = ", self%apply_constant_source
    write(*, '(a, a)') "source_file = ", self%source_file
    write(*, '(a, es9.3)') "constant_source_value = ", self%constant_source_value
    write(*, '(a, es9.3)') "source_scale_factor   = ", self%source_scale_factor
    write(*, '(a, i0)') "source_ilo = ", self%source_ilo
    write(*, '(a, i0)') "source_ihi = ", self%source_ihi
    write(*, '(a, i0)') "source_jlo = ", self%source_jlo
    write(*, '(a, i0)') "source_jhi = ", self%source_jhi

    write(*, *)
    write(*, '(a)') "[boundary_conditions]"
    write(*, '(a, a)') "bc_pressure_input_file = ", self%bc_pressure_input_file
    write(*, '(a, l1)') "apply_constant_bc_pressure = ", self%apply_constant_bc_pressure
    write(*, '(a, es9.3)') "constant_bc_pressure_value = ", self%constant_bc_pressure_value
    write(*, '(a, es9.3)') "bc_pressure_scale_factor = ", self%bc_pressure_scale_factor
    write(*, '(a, a)') "plus_x_bc = ", self%plus_x_bc
    write(*, '(a, a)') "minus_x_bc = ", self%minus_x_bc
    write(*, '(a, a)') "plus_y_bc  = ", self%plus_y_bc
    write(*, '(a, a)') "minus_y_bc = ", self%minus_y_bc

    write(*, *)
    write(*, '(a)') "[scheme]"
    write(*, '(a, l1)') "smooth_residuals = ", self%smooth_residuals
    write(*, '(a, a)') "evolution_operator_type = ", self%evolution_operator_type
    write(*, '(a, a)') "cell_reconstruction = ", self%cell_reconstruction
    write(*, '(a, a)') "edge_interpolation_scheme = ", self%edge_interpolation_scheme
    write(*, '(a, es9.3)') "tau = ", self%tau
    write(*, '(a, a)') "limiter = ", self%limiter

    write(*, *)
    write(*, '(a)') "[physics]"
    write(*, '(a, f5.3)') "polytropic_index = ", self%polytropic_index

    write(*, *)
    write(*, '(a)') "[io]"
    write(*, '(a, a)') "contour_io_format = ", self%contour_io_format
    write(*, '(a, es9.3)') "contour_interval_dt = ", self%contour_interval_dt
    write(*, '(a, l1)') "append_date_to_result_folder = ", self%append_date_to_result_folder
    write(*, '(a, l1)') "plot_reconstruction_states = ", self%plot_reconstruction_states
    write(*, '(a, l1)') "plot_reference_states = ", self%plot_reference_states
    write(*, '(a, l1)') "plot_evolved_states = ", self%plot_evolved_states
    write(*, '(a, l1)') "plot_64bit = ", self%plot_64bit
    write(*, '(a, l1)') "plot_ghost_cells = ", self%plot_ghost_cells

    write(*, '(a)') "==============="
    print *

  end subroutine display_config
end module mod_input
