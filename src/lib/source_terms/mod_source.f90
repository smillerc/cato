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

module mod_source
  !< Define the base class for all source term classes. These facilitate injecting a
  !< "source term" into the domain, e.g. energy, pressure, etc...

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_out => output_unit
  use mod_units
  use math_constants, only: pi
  use linear_interpolation_module, only: linear_interp_1d
  use mod_input, only: input_t
  use mod_grid_block, only: grid_block_t
  use mod_grid_block_2d, only: grid_block_2d_t
  use mod_error, only: error_msg
  use mod_nondimensionalization, only: t_0, rho_0, l_0, v_0, p_0, e_0

  implicit none

  type :: source_t

    integer(ik) :: io_unit = 0
    character(len=:), allocatable :: source_type !< what type of source, e.g. energy, force, etc.?
    character(len=:), allocatable :: source_geometry !< Shape of the source, e.g. uniform, constant, or 1D/2D gaussian

    real(rk) :: dt = 0.0_rk   !< Current time step
    real(rk) :: time = 0.0_rk !< Current time
    real(rk) :: max_time = 0.0_rk !< Max time in source (e.g. stop after this)

    real(rk) :: multiplier = 1.0_rk
    real(rk) :: nondim_scale_factor = 1.0_rk

    ! Bounds
    integer(ik) :: ilo = 0 !< lower i index bound on the source region
    integer(ik) :: jlo = 0 !< lower j index bound on the source region
    integer(ik) :: ihi = 0 !< upper i index bound on the source region
    integer(ik) :: jhi = 0 !< upper j index bound on the source region

    real(rk) :: xlo = 0.0 !< lower x extent bound on the source region
    real(rk) :: ylo = 0.0 !< lower y extent bound on the source region
    real(rk) :: xhi = 0.0 !< upper x extent bound on the source region
    real(rk) :: yhi = 0.0 !< upper y extent bound on the source region

    ! For gaussian source terms
    real(rk) :: fwhm_x = 0.0 !< full-width half-max value in the x-direction
    real(rk) :: fwhm_y = 0.0 !< full-width half-max value in the y-direction
    integer(ik) :: gaussian_order = 0

    real(rk) :: x_center = 0.0_rk !< x-center of the source term (if any)
    real(rk) :: y_center = 0.0_rk !< y-center of the source term (if any)

    logical :: constant_source = .false.
    real(rk) :: constant_source_value = 0.0_rk
    logical :: constant_in_x = .false. !< is the source term 1D in the x-direction, e.g. only depends on x
    logical :: constant_in_y = .false. !< is the source term 1D in the y-direction, e.g. only depends on y

    ! Temporal inputs
    type(linear_interp_1d) :: temporal_source_input
    character(len=:), allocatable :: input_filename

    ! Field data
    real(rk), allocatable, dimension(:, :) :: data !< (i,j); actual source data to apply
    real(rk), pointer, dimension(:, :) :: centroid_x => null() !< (i,j); cell x centroid used to apply positional source data
    real(rk), pointer, dimension(:, :) :: centroid_y => null() !< (i,j); cell y centroid used to apply positional source data
    real(rk), pointer, dimension(:, :) :: volume => null()
  contains
    private
    procedure, public :: integrate
    procedure, public :: read_input_file
    procedure, public :: get_desired_source_input
    procedure, public :: get_source_field
  endtype

contains

  function new_source(input, grid, time)
    !< Source constructor
    type(source_t), pointer :: new_source
    class(input_t), intent(in) :: input
    class(grid_block_t), intent(in), target :: grid
    real(rk), intent(in) :: time
    allocate(new_source)

    new_source%source_type = 'energy'
    new_source%time = time
    new_source%multiplier = input%source_scale_factor
    new_source%nondim_scale_factor = e_0 / t_0
    new_source%source_geometry = trim(input%source_geometry)
    new_source%constant_source = input%apply_constant_source

    select type(grid)
    class is(grid_block_2d_t)
      new_source%volume => grid%volume

      if(.not. new_source%constant_source) then
        new_source%input_filename = trim(input%source_file)
        call new_source%read_input_file()
      else
        new_source%constant_source_value = input%constant_source_value
      endif

      select case(new_source%source_geometry)
      case('uniform')
        ! case('constant_xy', '1d_gaussian', '2d_gaussian')
      case('constant_xy', '1d_gaussian')
        ! Get the cell centroid data since the source is position dependant
        new_source%centroid_x => grid%centroid_x
        new_source%centroid_y => grid%centroid_y
      case default
        call error_msg(message="Unsupported geometry type, must be one of ['uniform', 'constant_xy', '1d_gaussian']", &
             module_name='mod_source', class_name='source_t', procedure_name='new_source', file_name=__FILE__, line_number=__LINE__)
      endselect

      ! Determine the extents of the source term (if any)
      select case(new_source%source_geometry)
      case('constant_xy')
        new_source%xlo = input%source_xlo / l_0
        new_source%xhi = input%source_xhi / l_0
        new_source%ylo = input%source_ylo / l_0
        new_source%yhi = input%source_yhi / l_0
      case('1d_gaussian', '2d_gaussian')
        new_source%x_center = input%source_center_x / l_0
        new_source%y_center = input%source_center_y / l_0
        new_source%fwhm_x = input%source_gaussian_fwhm_x / l_0
        new_source%fwhm_y = input%source_gaussian_fwhm_y / l_0
        new_source%gaussian_order = input%source_gaussian_order

        if(new_source%gaussian_order < 1) then
          call error_msg(message="Source term gaussian order is invalid (< 1)", &
             module_name='mod_source', class_name='source_t', procedure_name='new_source', file_name=__FILE__, line_number=__LINE__)
        endif

        if(new_source%fwhm_x > 0.0_rk .and. new_source%fwhm_y <= 0.0_rk) then
          new_source%constant_in_y = .true.
        else if(new_source%fwhm_y > 0.0_rk .and. new_source%fwhm_x <= 0.0_rk) then
          new_source%constant_in_x = .true.
        endif
      endselect

      allocate(new_source%data, mold=grid%centroid_x)
      new_source%data = 0.0_rk
    endselect
  endfunction new_source

  subroutine read_input_file(self)
    !< Read the input file and initialize the linear iterpolated object. This
    !< makes it easy to get the input at any given time
    class(source_t), intent(inout) :: self

    real(rk), dimension(:), allocatable :: time ! time in seconds
    real(rk), dimension(:), allocatable :: source_data ! source data in cgs units
    character(len=300) :: line_buffer
    logical :: has_header_line
    logical :: file_exists
    integer(ik) :: input_unit, io_status, interp_status
    integer(ik) :: line, nlines

    file_exists = .false.
    has_header_line = .false.

    inquire(file=trim(self%input_filename), exist=file_exists)

    if(.not. file_exists) then
      error stop 'Error in source_t%read_input_file(); input file not found, exiting...'
    endif

    open(newunit=input_unit, file=trim(self%input_filename))
    nlines = 0
    do
      read(input_unit, *, iostat=io_status) line_buffer
      if(nlines == 0) then
        line_buffer = adjustl(line_buffer)
        if(line_buffer(1:1) == '#') has_header_line = .true.
      endif

      if(io_status /= 0) exit
      nlines = nlines + 1
    enddo
    close(input_unit)

    open(newunit=input_unit, file=trim(self%input_filename))
    if(has_header_line) then
      read(input_unit, *, iostat=io_status) ! skip the first line
      nlines = nlines - 1
    endif

    allocate(time(nlines))
    allocate(source_data(nlines))
    time = 0.0_rk
    source_data = 0.0_rk

    do line = 1, nlines
      read(input_unit, *, iostat=io_status) time(line), source_data(line)
      if(io_status /= 0) exit
    enddo
    close(input_unit)

    time = time / t_0
    self%max_time = maxval(time)

    source_data = source_data / self%nondim_scale_factor

    ! Initialize the linear interpolated data object so we can query the pressure at any time
    call self%temporal_source_input%initialize(time, source_data, interp_status)
    if(interp_status /= 0) error stop "Error initializing pressure_source_t%temporal_source_input"

    deallocate(time)
    deallocate(source_data)
  endsubroutine read_input_file

  real(rk) function get_desired_source_input(self) result(source_value)
    !< Get the desired source term amount. This can vary with time, so this can
    !< use linear interpolation.

    class(source_t), intent(inout) :: self
    integer(ik) :: interp_stat

    if(self%constant_source) then
      source_value = self%constant_source_value
    else
      if(self%time <= self%max_time) then
        call self%temporal_source_input%evaluate(x=self%time, &
                                                 f=source_value, &
                                                 istat=interp_stat)

        if(interp_stat /= 0) then
          write(std_out, '(a)') "Unable to interpolate value within "// &
            "source_t%get_desired_source_input()"
          error stop "Unable to interpolate value within "// &
            "source_t%get_desired_source_input()"
        endif
      else
        source_value = 0.0_rk
      endif
    endif

    if(source_value > 0.0_rk) then
      write(*, '(a, es10.3)') 'Applying source term of: ', source_value * self%nondim_scale_factor
    endif

  endfunction get_desired_source_input

  subroutine integrate(self, dt)
    !< Create the source field to be passed to the fluid class or others
    class(source_t), intent(inout) :: self
    real(rk), intent(in) :: dt

    self%dt = dt
    self%time = self%time + dt

    call self%get_source_field()
  endsubroutine integrate

  subroutine get_source_field(self)
    !< For the given time find the source strength and distribution based on the geometry type
    class(source_t), intent(inout) :: self
    real(rk) :: source_value
    real(rk) :: gauss_amplitude, c
    real(rk), parameter :: fwhm_to_c = 1.0_rk / (2.0_rk * sqrt(2.0_rk * log(2.0)))

    source_value = self%get_desired_source_input() * self%multiplier
    self%data = 0.0_rk

    if(abs(source_value) > 0.0_rk) then
      select case(self%source_geometry)
      case('uniform')
        self%data = source_value
      case('constant_xy')
        where(self%centroid_x < self%xhi .and. self%centroid_x > self%xlo .and. &
              self%centroid_y < self%yhi .and. self%centroid_x > self%ylo)
          self%data = 1.0_rk
        endwhere
      case('1d_gaussian')
        if(self%constant_in_y) then
          associate(x => self%centroid_x, x0 => self%x_center, &
                    fwhm => self%fwhm_x, order => self%gaussian_order)

            c = fwhm * fwhm_to_c
            gauss_amplitude = source_value / (sqrt(2.0_rk * pi) * abs(c))
            self%data = gauss_amplitude * exp(-((((x - x0)**2) / ((2.0_rk * c)**2))))

          endassociate
        else
          error stop "not done yet"
          associate(y => self%centroid_y, y0 => self%y_center, &
                    fwhm => self%fwhm_y, order => self%gaussian_order)
            self%data = exp(-((((y - y0)**2) / fwhm**2)**order))
          endassociate
        endif
      case('2d_gaussian')

        associate(x => self%centroid_x, x0 => self%x_center, fwhm_x => self%fwhm_x, &
                  y => self%centroid_y, y0 => self%y_center, fwhm_y => self%fwhm_y, &
                  order => self%gaussian_order)
          self%data = exp(-(((((x - x0)**2) / fwhm_x**2) + &
                             (((y - y0)**2) / fwhm_y**2))**order))
        endassociate
      endselect

      ! Remove tiny numbers
      where(abs(self%data) < tiny(1.0_rk)) self%data = 0.0_rk

      ! self%data = self%data / maxval(self%data)

      ! write(*,'(a, 10(es16.6))') "min/max: ", minval(self%data), maxval(self%data)
      ! write(*,'(a, es16.6)') "source_value", source_value
      ! write(*,'(a, es16.6)') "sum(self%data)", sum(self%data)
      ! write(*,'(a, es16.6)') "sum(self%volume, mask=abs(self%data) > 0.0_rk)", sum(self%volume, mask=abs(self%data) > 0.0_rk)
      self%data = self%data / sum(self%volume, mask=abs(self%data) > 0.0_rk)
      ! write(*,'(a, es16.6)') "self%data", sum(self%data)
    endif
    ! error stop
  endsubroutine get_source_field

  subroutine finalize(self)
    type(source_t), intent(inout) :: self
    logical :: is_open = .false.

    if(allocated(self%data)) deallocate(self%data)
    if(allocated(self%source_type)) deallocate(self%source_type)
    if(allocated(self%source_geometry)) deallocate(self%source_geometry)
    if(allocated(self%input_filename)) deallocate(self%input_filename)
    if(associated(self%centroid_x)) deallocate(self%centroid_x)
    if(associated(self%centroid_y)) deallocate(self%centroid_y)
    if(associated(self%volume)) deallocate(self%volume)
  endsubroutine finalize
endmodule mod_source
