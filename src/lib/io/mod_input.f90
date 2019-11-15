module mod_input

  use iso_fortran_env, only: ik => int32, rk => real64
  use cfgio_mod, only: cfg_t, parse_cfg

  implicit none

  private
  public :: input_t

  type :: input_t
    ! general
    character(:), allocatable :: title

    ! grid
    logical :: read_from_hdf5 = .false.
    character(:), allocatable :: grid_type
    real(rk) :: xmin
    real(rk) :: xmax
    real(rk) :: ymin
    real(rk) :: ymax
    integer(ik) :: ni
    integer(ik) :: nj

    ! boundary
    character(:), allocatable ::  plus_x_bc
    character(:), allocatable :: minus_x_bc

    character(:), allocatable ::  plus_y_bc
    character(:), allocatable :: minus_y_bc

    character(:), allocatable ::  plus_z_bc
    character(:), allocatable :: minus_z_bc

    ! io
    character(:), allocatable :: contour_io_format

    ! timing
    real(rk) :: max_time = 1
    real(rk) :: contour_interval = 0.5_rk

    ! physics
    real(rk) :: polytropic_index = 5.0_rk / 3.0_rk

    ! finite volume scheme specifics
    character(:), allocatable :: reconstruction_type
    real(rk) :: tau = 0.0_rk !< time increment for FVEG and FVLEG schemes

  contains
    procedure, public :: initialize
    procedure, public :: read_from_ini
  end type input_t

contains
  subroutine initialize(self, ni, nj, xmin, xmax, ymin, ymax)
    class(input_t), intent(inout) :: self
    integer(ik), intent(in) :: ni, nj
    real(rk), intent(in) :: xmin, xmax, ymin, ymax

    self%ni = ni
    self%nj = nj
    self%xmin = xmin
    self%xmax = xmax
    self%ymin = ymin
    self%ymax = ymax

  end subroutine initialize

  subroutine read_from_ini(self, filename)

    character(len=*), intent(in) :: filename
    class(input_t), intent(inout) :: self

    character(len=120) :: char_buffer
    type(cfg_t) :: cfg
    logical :: file_exists

    file_exists = .false.

    if(this_image() == 1) then
      print *, 'Reading input file: '//trim(filename)
    end if

    inquire(file=filename, EXIST=file_exists)
    if(.not. file_exists) error stop "Error: Input .ini file not found"

    cfg = parse_cfg(trim(filename))

    ! General
    call cfg%get("general", "title", char_buffer)
    self%title = trim(char_buffer)
    call cfg%get("general", "contour_interval", self%contour_interval)

    ! Physics
    call cfg%get("physics", "polytropic_index", self%polytropic_index)

    ! FV scheme
    call cfg%get("scheme", "reconstruction_type", char_buffer, 'piecewise_linear')
    self%reconstruction_type = trim(char_buffer)

    ! Grid
    call cfg%get("grid", "read_from_hdf5", self%read_from_hdf5, .false.)
    call cfg%get("grid", "grid_type", char_buffer, '2d_regular')
    self%grid_type = trim(char_buffer)

    call cfg%get("grid", "xmin", self%xmin)
    call cfg%get("grid", "xmax", self%xmax)
    call cfg%get("grid", "ymin", self%ymin)
    call cfg%get("grid", "ymax", self%ymax)

    ! boundary conditions
    call cfg%get("boundary_conditions", "plus_x", char_buffer, 'periodic')
    self%plus_x_bc = trim(char_buffer)
    call cfg%get("boundary_conditions", "minus_x", char_buffer, 'periodic')
    self%minus_x_bc = trim(char_buffer)

    call cfg%get("boundary_conditions", "plus_y", char_buffer, 'periodic')
    self%plus_y_bc = trim(char_buffer)
    call cfg%get("boundary_conditions", "minus_y", char_buffer, 'periodic')
    self%minus_y_bc = trim(char_buffer)

    call cfg%get("boundary_conditions", "plus_z", char_buffer, 'periodic')
    self%plus_z_bc = trim(char_buffer)
    call cfg%get("boundary_conditions", "minus_z", char_buffer, 'periodic')
    self%minus_z_bc = trim(char_buffer)

    ! Input/Output
    call cfg%get("io", "format", char_buffer, 'xdmf')
    self%contour_io_format = trim(char_buffer)

  end subroutine read_from_ini

end module mod_input
