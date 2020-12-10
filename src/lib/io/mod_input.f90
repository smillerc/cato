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

module mod_input

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, output_unit
  use mod_globals, only: debug_print
  use mod_error, only: error_msg
  use cfgio_mod, only: cfg_t, parse_cfg

  implicit none

  private
  public :: input_t

  type :: input_t
    ! general
    character(:), allocatable :: title  !< Name of the simulation
    character(:), allocatable :: unit_system !< Type of units to output (CGS, ICF, MKS, etc...)

    ! reference state
    logical :: non_dimensionalize = .true.
    real(rk) :: reference_pressure = 1.0_rk
    real(rk) :: reference_density = 1.0_rk
    real(rk) :: reference_mach = 1.0_rk

    ! grid
    character(len=32) :: grid_type = 'XY'  !< Structure/layout of the grid, e.g. '2d_regular'
    real(rk) :: xmin = 0.0_rk   !< Minimum extent of the grid in x (ignored for .h5 initial grids)
    real(rk) :: xmax = 0.0_rk   !< Maximum extent of the grid in x (ignored for .h5 initial grids)
    real(rk) :: ymin = 0.0_rk   !< Minimum extent of the grid in y (ignored for .h5 initial grids)
    real(rk) :: ymax = 0.0_rk   !< Maximum extent of the grid in y (ignored for .h5 initial grids)
    integer(ik) :: ni_nodes = 0 !< # of i nodes (not including ghost) (ignored for .h5 initial grids)
    integer(ik) :: nj_nodes = 0 !< # of j nodes (not including ghost) (ignored for .h5 initial grids)
    integer(ik) :: n_ghost_layers = 2 !< # of ghost layers to use (based on spatial reconstruction order)

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

    real(rk) :: source_xlo = 0.0_rk
    real(rk) :: source_xhi = 0.0_rk
    real(rk) :: source_ylo = 0.0_rk
    real(rk) :: source_yhi = 0.0_rk

    real(rk) :: source_center_x = 0.0_rk
    real(rk) :: source_center_y = 0.0_rk
    real(rk) :: source_gaussian_fwhm_x = 0.0_rk
    real(rk) :: source_gaussian_fwhm_y = 0.0_rk
    integer(ik) :: source_gaussian_order = 1

    character(:), allocatable :: source_geometry

    ! io
    character(:), allocatable :: contour_io_format !< e.g. 'xdmf'
    logical :: append_date_to_result_folder = .false.
    logical :: plot_reconstruction_states = .false.
    logical :: plot_reference_states = .false.
    logical :: plot_evolved_states = .false.
    logical :: plot_64bit = .true.
    logical :: plot_ghost_cells = .true.

    ! timing
    real(rk) :: max_time = 1.0_rk
    real(rk) :: cfl = 0.1_rk ! Courant–Friedrichs–Lewy condition
    real(rk) :: initial_delta_t = 0.0_rk
    logical :: use_constant_delta_t = .false.
    real(rk) :: contour_interval_dt = 0.5_rk
    integer(ik) :: max_iterations = huge(1)
    character(:), allocatable :: time_integration_strategy !< How is time integration handled? e.g. 'rk2', 'rk4', etc.
    logical :: smooth_residuals = .true.

    ! physics
    real(rk) :: polytropic_index = 5.0_rk / 3.0_rk !< e.g. gamma for the simulated gas

    ! finite volume scheme specifics
    character(len=32) :: flux_solver = 'AUSMPW+'          !< flux solver, ('FVLEG', 'AUSM+-up', etc.)
    character(len=32) :: spatial_reconstruction = 'MUSCL' !< (MUSCL, green_gauss) How are the edge values interpolated?
    logical :: apply_low_mach_fix = .true.                !< some flux solvers have this option
    character(len=32) :: limiter = 'minmod'               !< Flux limiter, e.g. minmod, superbee, TVD3, MLP3, MLP5, etc.

    ! AUSM solver specifics, see the AUSM solver packages for more details.
    ! These are only read in if the flux solver is from the AUSM family
    real(rk) :: ausm_beta = 1.0_rk / 8.0_rk             !< beta parameter
    real(rk) :: ausm_pressure_diffusion_coeff = 0.25_rk !< K_p; pressure diffusion coefficient
    real(rk) :: ausm_pressure_flux_coeff = 0.75_rk      !< K_u; pressure flux coefficient
    real(rk) :: ausm_sonic_point_sigma = 1.0_rk         !< sigma; another pressure diffusion coefficient

    real(rk) :: tau = 1.0e-5_rk !< time increment for FVEG and FVLEG schemes

    ! debug
    logical :: plot_limiters = .false.
    logical :: plot_gradients = .false.

  contains
    procedure, public :: initialize
    procedure, public :: read_from_ini
    procedure, public :: display_config
    final :: finalize
  endtype input_t

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
    self%spatial_reconstruction = 'MUSCL'
    self%limiter = 'minmod'

    self%plus_x_bc = 'periodic'
    self%minus_x_bc = 'periodic'
    self%plus_y_bc = 'periodic'
    self%minus_y_bc = 'periodic'

  endsubroutine initialize

  subroutine finalize(self)
    type(input_t), intent(inout) :: self

    call debug_print('Running input_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%time_integration_strategy)) deallocate(self%time_integration_strategy)
    if(allocated(self%contour_io_format)) deallocate(self%contour_io_format)
    if(allocated(self%source_file)) deallocate(self%source_file)
    if(allocated(self%bc_pressure_input_file)) deallocate(self%bc_pressure_input_file)
    if(allocated(self%restart_file)) deallocate(self%restart_file)
    if(allocated(self%initial_condition_file)) deallocate(self%initial_condition_file)
    if(allocated(self%title)) deallocate(self%title)
    if(allocated(self%unit_system)) deallocate(self%unit_system)
  endsubroutine

  subroutine read_from_ini(self, filename)

    character(len=*), intent(in) :: filename
    class(input_t), intent(inout) :: self

    character(len=120) :: char_buffer
    type(cfg_t) :: cfg
    logical :: file_exists
    integer(ik) :: required_n_ghost_layers = 2

    file_exists = .false.

    if(this_image() == 1) write(output_unit, '(a)') 'Reading input file: '//trim(filename)

    inquire(file=filename, exist=file_exists)
    if(.not. file_exists) then
      call error_msg(module_name='mod_input', class_name="input_t", procedure_name='read_from_ini', &
                     message="Input .ini file not found", &
                     file_name=__FILE__, line_number=__LINE__)
    endif

    cfg = parse_cfg(trim(filename))

    ! General
    call cfg%get("general", "title", char_buffer)
    self%title = trim(char_buffer)

    call cfg%get("general", "units", char_buffer, 'cgs')
    self%unit_system = trim(char_buffer)

    ! Reference state
    call cfg%get("reference_state", "non_dimensionalize", self%non_dimensionalize, .true.)

    if(self%non_dimensionalize) then
      call cfg%get("reference_state", "reference_pressure", self%reference_pressure)
      call cfg%get("reference_state", "reference_density", self%reference_density)
      call cfg%get("reference_state", "reference_mach", self%reference_mach, 1.0_rk)
    endif
    ! Time
    call cfg%get("time", "max_time", self%max_time)
    call cfg%get("time", "initial_delta_t", self%initial_delta_t, 0.0_rk)
    call cfg%get("time", "use_constant_delta_t", self%use_constant_delta_t, .false.)
    call cfg%get("time", "cfl", self%cfl, 0.1_rk)
    call cfg%get("time", "integration_strategy", char_buffer)
    call cfg%get("time", "max_iterations", self%max_iterations, huge(1))
    self%time_integration_strategy = trim(char_buffer)

    ! Physics
    call cfg%get("physics", "polytropic_index", self%polytropic_index)

    ! Scheme
    call cfg%get("scheme", "flux_solver", self%flux_solver)

    if(trim(self%flux_solver) == 'AUSM+-up') then
      call cfg%get("ausm", "beta", self%ausm_beta, 1.0_rk / 8.0_rk)
      call cfg%get("ausm", "pressure_diffusion_coeff", self%ausm_pressure_diffusion_coeff, 0.25_rk)
      call cfg%get("ausm", "pressure_flux_coeff", self%ausm_pressure_flux_coeff, 0.75_rk)
      call cfg%get("ausm", "sonic_point_sigma", self%ausm_sonic_point_sigma, 1.0_rk)
    endif

    call cfg%get("scheme", "smooth_residuals", self%smooth_residuals, .true.)

    ! Spatial reconstruction
    call cfg%get("scheme", "spatial_reconstruction", self%spatial_reconstruction, 'MUSCL')
    call cfg%get("scheme", "limiter", self%limiter, 'minmod')
    call cfg%get("scheme", "apply_low_mach_fix", self%apply_low_mach_fix, .true.)

    ! Restart files
    call cfg%get("restart", "restart_from_file", self%restart_from_file, .false.)

    if(this_image() == 1) write(*, '(a, l2)') 'Reading from restart?', self%restart_from_file

    if(self%restart_from_file) then
      call cfg%get("restart", "restart_file", char_buffer)
      self%restart_file = trim(char_buffer)

      inquire(file=self%restart_file, exist=file_exists)
      if(.not. file_exists) then
        call error_msg(module_name='mod_input', class_name="input_t", procedure_name='read_from_ini', &
                       message='Restart file not found: "'//trim(self%restart_file)//'"', &
                       file_name=__FILE__, line_number=__LINE__)
      endif
    endif

    ! Initial Conditions
    if(.not. self%restart_from_file) then
      call cfg%get("initial_conditions", "read_from_file", self%read_init_cond_from_file, .true.)

      if(self%read_init_cond_from_file) then
        call cfg%get("initial_conditions", "initial_condition_file", char_buffer)
        self%initial_condition_file = trim(char_buffer)

        inquire(file=self%initial_condition_file, exist=file_exists)
        if(.not. file_exists) then
          call error_msg(module_name='mod_input', class_name="input_t", procedure_name='read_from_ini', &
                         message='Initial conditions file not found: "'//trim(self%initial_condition_file)//'"', &
                         file_name=__FILE__, line_number=__LINE__)
        endif

      endif

      call cfg%get("initial_conditions", "init_density", self%init_density, 0.0_rk)
      call cfg%get("initial_conditions", "init_x_velocity", self%init_x_velocity, 0.0_rk)
      call cfg%get("initial_conditions", "init_y_velocity", self%init_y_velocity, 0.0_rk)
      call cfg%get("initial_conditions", "init_pressure", self%init_pressure, 0.0_rk)
    endif

    ! Grid
    select case(trim(self%limiter))
    case('minmod', 'superbee', 'van_leer', 'none', 'TVD2', 'TVD3', 'MLP3', 'eMLP3')
      required_n_ghost_layers = 2
      call cfg%get("grid", "n_ghost_layers", self%n_ghost_layers, required_n_ghost_layers)
    case('TVD5', 'MLP5', 'eMLP5')
      required_n_ghost_layers = 3
      call cfg%get("grid", "n_ghost_layers", self%n_ghost_layers, required_n_ghost_layers)
    case default
      call error_msg(module_name='mod_input', class_name="input_t", procedure_name='read_from_ini', &
                     message="Unknown edge interpolation scheme, must be one of the following: "// &
                     "'TVD2', 'TVD3', 'TVD5', 'MLP3', or 'MLP5'", &
                     file_name=__FILE__, line_number=__LINE__)
    endselect

    if(self%n_ghost_layers /= required_n_ghost_layers) then
      call error_msg(module_name='mod_input', class_name="input_t", procedure_name='read_from_ini', &
                     message="The number of required ghost cell layers doesn't match the edge interpolation order", &
                     file_name=__FILE__, line_number=__LINE__)
    endif

    call cfg%get("grid", "grid_type", char_buffer, '2d_regular')
    self%grid_type = trim(char_buffer)

    if(.not. self%read_init_cond_from_file .and. .not. self%restart_from_file) then
      call cfg%get("grid", "ni_nodes", self%ni_nodes)
      call cfg%get("grid", "xmin", self%xmin)
      call cfg%get("grid", "xmax", self%xmax)

      call cfg%get("grid", "nj_nodes", self%nj_nodes)
      call cfg%get("grid", "ymin", self%ymin)
      call cfg%get("grid", "ymax", self%ymax)
    endif

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
      call cfg%get("boundary_conditions", "bc_density", self%bc_density, 0.0_rk)

      if(self%apply_constant_bc_pressure) then
        call cfg%get("boundary_conditions", "constant_bc_pressure_value", self%constant_bc_pressure_value)
      else
        call cfg%get("boundary_conditions", "bc_pressure_input_file", char_buffer)
        self%bc_pressure_input_file = trim(char_buffer)
      endif

      call cfg%get("boundary_conditions", "bc_pressure_scale_factor", self%bc_pressure_scale_factor, 1.0_rk)
    endif

    ! Source terms
    call cfg%get('source_terms', 'enable_source_terms', self%enable_source_terms, .false.)

    if(self%enable_source_terms) then

      call cfg%get('source_terms', 'geometry', char_buffer)
      self%source_geometry = trim(char_buffer)

      call cfg%get('source_terms', 'apply_constant_source', self%apply_constant_source, .false.)
      call cfg%get('source_terms', 'source_file', char_buffer)
      self%source_file = trim(char_buffer)

      call cfg%get('source_terms', 'source_term_type', char_buffer)
      self%source_term_type = trim(char_buffer)

      call cfg%get('source_terms', 'constant_source_value', self%constant_source_value, 0.0_rk)
      call cfg%get('source_terms', 'source_scale_factor', self%source_scale_factor, 1.0_rk)

      select case(self%source_geometry)
      case('constant_xy')
        call cfg%get('source_terms', 'xlo', self%source_xlo)
        call cfg%get('source_terms', 'xhi', self%source_xhi)
        call cfg%get('source_terms', 'ylo', self%source_ylo)
        call cfg%get('source_terms', 'yhi', self%source_yhi)
      case('1d_gaussian', '2d_gaussian')
        call cfg%get('source_terms', 'center_x', self%source_center_x, 0.0_rk)
        call cfg%get('source_terms', 'center_y', self%source_center_y, 0.0_rk)
        call cfg%get('source_terms', 'gaussian_fwhm_x', self%source_gaussian_fwhm_x, 0.0_rk)
        call cfg%get('source_terms', 'gaussian_fwhm_y', self%source_gaussian_fwhm_y, 0.0_rk)
        call cfg%get('source_terms', 'gaussian_order', self%source_gaussian_order, 1)
      endselect
    endif

    ! Input/Output
    call cfg%get("io", "format", char_buffer, 'xdmf')
    self%contour_io_format = trim(char_buffer)

    call cfg%get("io", "contour_interval_dt", self%contour_interval_dt, 0.1_rk)
    call cfg%get("io", "append_date_to_result_folder", self%append_date_to_result_folder, .false.)
    call cfg%get("io", "plot_reconstruction_states", self%plot_reconstruction_states, .false.)
    call cfg%get("io", "plot_reference_states", self%plot_reference_states, .false.)
    call cfg%get("io", "plot_evolved_states", self%plot_evolved_states, .false.)
    call cfg%get("io", "plot_ghost_cells", self%plot_ghost_cells, .true.)
    call cfg%get("io", "plot_64bit", self%plot_64bit, .true.)

  endsubroutine read_from_ini

  subroutine display_config(self)
    class(input_t), intent(in) :: self

    if(this_image() == 1) then
      write(*, '(a)') "Input settings:"
      write(*, '(a)') "==============="
      write(*, *)
      write(*, '(a)') "[general]"
      write(*, '(3(a))') "title = '", trim(self%title), "'"
      write(*, '(3(a))') "unit_system = '", trim(self%unit_system), "'"

      write(*, *)
      write(*, '(a)') "[time]"
      write(*, '(a, es10.3)') "max_time = ", self%max_time
      write(*, '(a, f5.3)') "cfl = ", self%cfl
      write(*, '(a, l1)') "use_constant_delta_t = ", self%use_constant_delta_t
      write(*, '(a, es10.3)') "initial_delta_t = ", self%initial_delta_t
      write(*, '(a, i0)') "max_iterations = ", self%max_iterations
      write(*, '(3(a))') "time_integration_strategy = '", trim(self%time_integration_strategy), "'"

      write(*, *)
      write(*, '(a)') "[grid]"
      write(*, '(a, a)') "grid_type = ", self%grid_type
      write(*, '(a, es10.3)') "xmin = ", self%xmin
      write(*, '(a, es10.3)') "xmax = ", self%xmax
      write(*, '(a, es10.3)') "ymin = ", self%ymin
      write(*, '(a, es10.3)') "ymax = ", self%ymax
      write(*, '(a, i0)') "ni_nodes = ", self%ni_nodes
      write(*, '(a, i0)') "nj_nodes = ", self%nj_nodes
      write(*, '(a, i0)') "n_ghost_layers = ", self%n_ghost_layers

      write(*, *)
      write(*, '(a)') "[reference_state]"
      write(*, '(a, es10.3)') "reference_pressure = ", self%reference_pressure
      write(*, '(a, es10.3)') "reference_density = ", self%reference_density

      write(*, *)
      write(*, '(a)') "[initial_conditions]"
      write(*, '(a, a)') "initial_condition_file = ", self%initial_condition_file
      write(*, '(a, l1)') "read_init_cond_from_file = ", self%read_init_cond_from_file
      write(*, '(a, es10.3)') "init_x_velocity = ", self%init_x_velocity
      write(*, '(a, es10.3)') "init_y_velocity = ", self%init_y_velocity
      write(*, '(a, es10.3)') "init_density  = ", self%init_density
      write(*, '(a, es10.3)') "init_pressure = ", self%init_pressure

      write(*, *)
      write(*, '(a)') "[restart]"
      write(*, '(a, l1)') "restart_from_file = ", self%restart_from_file
      write(*, '(3(a))') "restart_file = '", trim(self%restart_file), "'"

      write(*, *)
      write(*, '(a)') "[source_terms]"
      write(*, '(a, l1)') "enable_source_terms = ", self%enable_source_terms
      write(*, '(a, a)') "source_term_type = ", self%source_term_type
      write(*, '(a, l1)') "apply_constant_source = ", self%apply_constant_source
      write(*, '(a, a)') "source_file = ", self%source_file
      write(*, '(a, es10.3)') "constant_source_value = ", self%constant_source_value
      write(*, '(a, es10.3)') "source_scale_factor   = ", self%source_scale_factor
      write(*, '(a, g0.3)') "source_ilo = ", self%source_xlo
      write(*, '(a, g0.3)') "source_ihi = ", self%source_xhi
      write(*, '(a, g0.3)') "source_jlo = ", self%source_ylo
      write(*, '(a, g0.3)') "source_jhi = ", self%source_yhi

      write(*, *)
      write(*, '(a)') "[boundary_conditions]"
      write(*, '(a, a)') "bc_pressure_input_file = ", self%bc_pressure_input_file
      write(*, '(a, l1)') "apply_constant_bc_pressure = ", self%apply_constant_bc_pressure
      write(*, '(a, es10.3)') "constant_bc_pressure_value = ", self%constant_bc_pressure_value
      write(*, '(a, es10.3)') "bc_pressure_scale_factor = ", self%bc_pressure_scale_factor
      write(*, '(a, a)') "plus_x_bc = ", self%plus_x_bc
      write(*, '(a, a)') "minus_x_bc = ", self%minus_x_bc
      write(*, '(a, a)') "plus_y_bc  = ", self%plus_y_bc
      write(*, '(a, a)') "minus_y_bc = ", self%minus_y_bc

      write(*, *)
      write(*, '(a)') "[scheme]"
      write(*, '(a, l1)') "smooth_residuals = ", self%smooth_residuals
      write(*, '(a, a)') "flux_solver = ", self%flux_solver
      write(*, '(a, a)') "spatial_reconstruction = ", self%spatial_reconstruction
      write(*, '(a, a)') "limiter = ", self%limiter
      write(*, '(a, es10.3)') "tau = ", self%tau

      write(*, *)
      write(*, '(a)') "[physics]"
      write(*, '(a, f5.3)') "polytropic_index = ", self%polytropic_index

      write(*, *)
      write(*, '(a)') "[io]"
      write(*, '(a, a)') "contour_io_format = ", self%contour_io_format
      write(*, '(a, es10.3)') "contour_interval_dt = ", self%contour_interval_dt
      write(*, '(a, l1)') "append_date_to_result_folder = ", self%append_date_to_result_folder
      write(*, '(a, l1)') "plot_reconstruction_states = ", self%plot_reconstruction_states
      write(*, '(a, l1)') "plot_reference_states = ", self%plot_reference_states
      write(*, '(a, l1)') "plot_evolved_states = ", self%plot_evolved_states
      write(*, '(a, l1)') "plot_64bit = ", self%plot_64bit
      write(*, '(a, l1)') "plot_ghost_cells = ", self%plot_ghost_cells

      write(*, '(a)') "==============="
      print *
    endif

  endsubroutine display_config
endmodule mod_input
