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
  use collectives, only: min_to_all, max_to_all
  use mod_field, only: field_2d_t, field_2d
  use math_constants, only: pi
  use linear_interpolation_module, only: linear_interp_1d
  use mod_input, only: input_t
  use mod_grid_block, only: grid_block_t
  use mod_grid_block_2d, only: grid_block_2d_t
  use mod_error, only: error_msg
  use mod_nondimensionalization

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

    logical :: this_image_deposits = .false.
    ! Temporal inputs
    type(linear_interp_1d) :: temporal_source_input
    character(len=:), allocatable :: input_filename

    ! Field data
    ! real(rk), allocatable, dimension(:, :) :: data !< (i,j); actual source data to apply
    type(field_2d_t) :: source    !< (i, j); Conserved quantities

    real(rk), pointer, dimension(:, :) :: centroid_x => null() !< (i,j); cell x centroid used to apply positional source data
    real(rk), pointer, dimension(:, :) :: centroid_y => null() !< (i,j); cell y centroid used to apply positional source data
    real(rk), pointer, dimension(:, :) :: volume => null()
    integer(ik) :: i_dep_center = 0
    integer(ik) :: j_dep_center = 0
    integer(ik), dimension(2) :: i_dep_range = 0 !< (imin, imax); i-index range for energy deposition
    integer(ik), dimension(2) :: j_dep_range = 0 !< (jmin, jmax); j-index range for energy deposition
  contains
    private
    procedure, public :: integrate
    procedure, public :: read_input_file
    procedure, public :: get_desired_source_input
    procedure, public :: get_source_field
    procedure, public :: find_deposition_location
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
    new_source%nondim_scale_factor = energy_to_nondim
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
        new_source%xlo = input%source_xlo * len_to_nondim
        new_source%xhi = input%source_xhi * len_to_nondim
        new_source%ylo = input%source_ylo * len_to_nondim
        new_source%yhi = input%source_yhi * len_to_nondim
      case('1d_gaussian', '2d_gaussian')
        new_source%x_center = input%source_center_x  * len_to_nondim
        new_source%y_center = input%source_center_y  * len_to_nondim
        new_source%fwhm_x = input%source_gaussian_fwhm_x  * len_to_nondim
        new_source%fwhm_y = input%source_gaussian_fwhm_y  * len_to_nondim
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

      new_source%source = field_2d(name='source_term', long_name='Source Term', &
                             descrip='Source Term', units='', &
                             global_dims=grid%global_dims, n_halo_cells=input%n_ghost_layers)
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

    time = time * t_to_nondim
    self%max_time = maxval(time)

    source_data = source_data * self%nondim_scale_factor

    ! Initialize the linear interpolated data object so we can query the pressure at any time
    call self%temporal_source_input%initialize(time, source_data, interp_status)
    if(interp_status /= 0) error stop "Error initializing source_t%temporal_source_input"

    deallocate(time)
    deallocate(source_data)
  endsubroutine read_input_file

  real(rk) function get_desired_source_input(self, t) result(source_value)
    !< Get the desired source term amount. This can vary with time, so this can
    !< use linear interpolation.

    class(source_t), intent(inout) :: self
    real(rk), intent(in) :: t !< time
    integer(ik) :: interp_stat

    if(self%constant_source) then
      source_value = self%constant_source_value
    else
      if(t <= self%max_time) then
        call self%temporal_source_input%evaluate(x=t, &
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


  endfunction get_desired_source_input

  subroutine integrate(self, dt, density)
    !< Create the source field to be passed to the fluid class or others
    class(source_t), intent(inout) :: self
    class(field_2d_t), intent(in) :: density
    real(rk), intent(in) :: dt

    self%dt = dt
    self%time = self%time + dt

    call self%find_deposition_location(density)
    call self%get_source_field(density)

    ! self%source%data = self%source%data / self%volume

  endsubroutine integrate

  subroutine find_deposition_location(self, density)
    class(source_t), intent(inout) :: self
    class(field_2d_t), intent(in) :: density
    integer(ik) :: i, ilo, ihi, j, jlo, jhi, i_dep
    integer(ik) :: imax, imin
    integer(ik), parameter :: I_WINDOW = 100 !< index window in which we look to deposit energy
    real(rk) :: x_center
    real(rk) :: critical_density !< where do we deposit the energy TODO: make this a user input

    critical_density = 0.1_rk * density_to_nondim ! 0.10 g/cc to nondimensionalize

    self%this_image_deposits = .false.
    imax = 0
    imin = 0
    x_center = 0.0_rk

    ihi = density%ihi
    ilo = density%ilo
    jhi = density%jhi
    jlo = density%jlo

    ! Scan in from the right; this is only for a "laser" from the +x boundary
    ! Each image does this, but not all actually deposit any energy
    do i = ihi, ilo, -1
      if (any(density%data(i, jlo:jhi) > critical_density)) then
        x_center = maxval(self%centroid_x(i,:))
        i_dep = i
        exit
      endif
    enddo

    self%x_center = max_to_all(x_center)
    i_dep = max_to_all(i_dep) ! broadcast the true i deposit index everywhere
    self%i_dep_range = [ilo, ihi] ! default to the entire range

    ! Check to see if the index window is in the current image's domain
    imin = max(i_dep - I_WINDOW, ilo)
    imax = min(i_dep + I_WINDOW, ihi)
    
    if (imin > ihi .or. imax < ilo) then
      self%this_image_deposits = .false.
    else
      self%this_image_deposits = .true.
      self%i_dep_range = [imin, imax]
    endif
  end subroutine

  subroutine get_source_field(self, density)
    !< For the given time find the source strength and distribution based on the geometry type
    class(source_t), intent(inout) :: self
    class(field_2d_t), intent(in) :: density
    ! real(rk) :: source_value
    real(rk) :: gauss_amplitude, c
    real(rk) :: p1, p2
    real(rk) :: e_input !< input energy
    real(rk), parameter :: fwhm_to_c = 1.0_rk / (2.0_rk * sqrt(2.0_rk * log(2.0_rk)))
    integer(ik) :: ilo, ihi, jlo, jhi

    ihi = density%ihi
    ilo = density%ilo
    jhi = density%jhi
    jlo = density%jlo


    ! The source term in the Euler equations is Q = [0, rho f_x, rho f_y, q_dot_h]
    ! q_dot_h is the heat transfer per unit mass, or energy per mass. In cgs this
    ! has units of erg/g. We are given a power from the input, which is erg/s
    e_input = self%get_desired_source_input(t=self%time) * self%multiplier
    ! p1 = self%get_desired_source_input(t=self%time) * self%multiplier
    ! p2 = self%get_desired_source_input(t=self%time + self%dt) * self%multiplier
    ! e_input = ((p1 + p2)* 0.5_rk) * self%dt
    self%source%data = 0.0_rk

    if(e_input > 0.0_rk) then
      if (this_image() == 1) write(*, '(2(a, es10.3))') 'Applying source term of [dim version]: ', &
      e_input / self%nondim_scale_factor, " at", self%x_center * len_to_dim
    endif
    
    if(abs(e_input) > 0.0_rk .and. self%this_image_deposits) then
      select case(self%source_geometry)
      case('uniform')
        self%source%data = e_input
      case('constant_xy')
        where(self%centroid_x < self%xhi .and. self%centroid_x > self%xlo .and. &
              self%centroid_y < self%yhi .and. self%centroid_x > self%ylo)
          self%source%data = 1.0_rk
        endwhere
      case('1d_gaussian')
        if(self%constant_in_y) then

          associate(imin => self%i_dep_range(1), imax => self%i_dep_range(2), &
                    x => self%centroid_x(self%i_dep_range(1):self%i_dep_range(2),:), x0 => self%x_center, &
                    fwhm => self%fwhm_x, order => self%gaussian_order)

            ! Make a gaussian
            self%source%data(imin:imax,:) = exp((-4.0_rk * log(2.0_rk) * (x - x0)**2) / fwhm**2)
            
            ! Filter out the small values at the fringes
            where(abs(self%source%data) < 1e-6_rk) self%source%data = 0.0_rk

            self%source%data(imin:imax,:) = self%source%data(imin:imax,:) * e_input

            !  / density%data(imin:imax,:) / self%volume(imin:imax,:)
            ! The source term is per unit mass, e.g. erg/gram -> erg * cm^3 / g / cm^3
            
          endassociate
        else
          error stop "not done yet"
          associate(y => self%centroid_y, y0 => self%y_center, &
                    fwhm => self%fwhm_y, order => self%gaussian_order)
            self%source%data = exp(-((((y - y0)**2) / fwhm**2)**order))
          endassociate
        endif
      case('2d_gaussian')

        associate(x => self%centroid_x, x0 => self%x_center, fwhm_x => self%fwhm_x, &
                  y => self%centroid_y, y0 => self%y_center, fwhm_y => self%fwhm_y, &
                  order => self%gaussian_order)
          self%source%data = exp(-(((((x - x0)**2) / fwhm_x**2) + &
                             (((y - y0)**2) / fwhm_y**2))**order))
        endassociate
      endselect
    endif
  endsubroutine get_source_field

  subroutine finalize(self)
    type(source_t), intent(inout) :: self
    logical :: is_open = .false.

    if(allocated(self%source_type)) deallocate(self%source_type)
    if(allocated(self%source_geometry)) deallocate(self%source_geometry)
    if(allocated(self%input_filename)) deallocate(self%input_filename)
    if(associated(self%centroid_x)) deallocate(self%centroid_x)
    if(associated(self%centroid_y)) deallocate(self%centroid_y)
    if(associated(self%volume)) deallocate(self%volume)
  endsubroutine finalize
endmodule mod_source
