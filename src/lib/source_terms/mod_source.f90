module mod_source
  !< Define the base class for all source term classes. These facilitate injecting a
  !< "source term" into the domain, e.g. energy, pressure, etc...

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_out => output_unit
  use mod_units
  use linear_interpolation_module, only: linear_interp_1d

  implicit none

  type, abstract :: source_t
    character(len=:), allocatable :: source_type
    real(rk) :: scale_factor = 1.0_rk !< User scale factor do scale energy up/down
    real(rk), private :: time = 0.0_rk ! Solution time (for time dependent bc's)
    real(rk), private :: max_time = 0.0_rk !< Max time in source (e.g. stop after this)
    integer(ik) :: ilo = 0 !< Index to apply source term at
    integer(ik) :: jlo = 0 !< Index to apply source term at
    integer(ik) :: ihi = 0 !< Index to apply source term at
    integer(ik) :: jhi = 0 !< Index to apply source term at
    integer(ik) :: io_unit !< file that the source value is written out to (time,source) pairs

    logical :: constant_source = .false.
    real(rk) :: source_input = 0.0_rk

    ! Temporal inputs
    type(linear_interp_1d) :: temporal_source_input
    character(len=100) :: input_filename

  contains
    procedure, public :: set_time
    procedure, public :: get_time
    procedure, public :: get_application_bounds
    procedure(apply_source), deferred :: apply_source
    procedure(copy_source), public, deferred :: copy
    procedure, public :: read_input_file
    procedure, public :: get_desired_source_value
    generic :: assignment(=) => copy
  end type

  abstract interface
    subroutine copy_source(out_source, in_source)
      import :: source_t
      class(source_t), intent(in) :: in_source
      class(source_t), intent(inout) :: out_source
    end subroutine copy_source

    subroutine apply_source(self, conserved_vars, lbounds, time)
      import :: source_t, rk, ik
      class(source_t), intent(inout) :: self
      integer(ik), dimension(3), intent(in) :: lbounds
      real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(inout) :: conserved_vars
      real(rk), intent(in) :: time
    end subroutine apply_source
  end interface
contains
  subroutine set_time(self, time)
    class(source_t), intent(inout) :: self
    real(rk), intent(in) :: time
    self%time = time
  end subroutine

  function get_time(self) result(time)
    class(source_t), intent(in) :: self
    real(rk) :: time
    time = self%time
  end function

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
    end if

    open(newunit=input_unit, file=trim(self%input_filename))
    nlines = 0
    do
      read(input_unit, *, iostat=io_status) line_buffer
      if(nlines == 0) then
        line_buffer = adjustl(line_buffer)
        if(line_buffer(1:1) == '#') has_header_line = .true.
      end if

      if(io_status /= 0) exit
      nlines = nlines + 1
    end do
    close(input_unit)

    open(newunit=input_unit, file=trim(self%input_filename))
    if(has_header_line) then
      read(input_unit, *, iostat=io_status) ! skip the first line
      nlines = nlines - 1
    end if

    allocate(time(nlines))
    allocate(source_data(nlines))
    time = 0.0_rk
    source_data = 0.0_rk

    do line = 1, nlines
      read(input_unit, *, iostat=io_status) time(line), source_data(line)
      if(io_status /= 0) exit
    end do
    close(input_unit)

    self%max_time = maxval(time)

    ! Initialize the linear interpolated data object so we can query the pressure at any time
    call self%temporal_source_input%initialize(time, source_data, interp_status)
    if(interp_status /= 0) error stop "Error initializing pressure_source_t%temporal_source_input"

    deallocate(time)
    deallocate(source_data)
  end subroutine read_input_file

  real(rk) function get_desired_source_value(self) result(source_value)
    !< Get the desired source term amount. This can vary with time, so this can
    !< use linear interpolation.

    class(source_t), intent(inout) :: self
    integer(ik) :: interp_stat
    real(rk) :: t

    t = self%get_time()
    if(self%constant_source) then
      source_value = self%source_input
    else
      if(t <= self%max_time) then
        call self%temporal_source_input%evaluate(x=t, &
                                                 f=source_value, &
                                                 istat=interp_stat)

        if(interp_stat /= 0) then
          write(std_out, '(a)') "Unable to interpolate value within "// &
            "source_t%get_desired_source_value()"
          error stop "Unable to interpolate value within "// &
            "source_t%get_desired_source_value()"
        end if
      else
        source_value = 0.0_rk
      end if
    end if

    source_value = source_value * self%scale_factor

    if(source_value > 0.0_rk) write(*, '(a, es10.3)') 'Applying source term of: ', source_value

    write(self%io_unit, '(2(es14.5, 1x))') t * io_time_units, source_value

  end function get_desired_source_value

  function get_application_bounds(self, lbounds, ubounds) result(ranges)
    !< Get the i and j ranges where the source is to be applied

    class(source_t), intent(in) :: self
    integer(ik), dimension(3), intent(in) :: lbounds
    integer(ik), dimension(3), intent(in) :: ubounds
    integer(ik), dimension(4) :: ranges !< [ilo, ihi, jlo, jhi]

    associate(ilo=>ranges(1), ihi=>ranges(2), &
              jlo=>ranges(3), jhi=>ranges(4))

      ilo = self%ilo
      ihi = self%ihi
      jlo = self%jlo
      jhi = self%jhi

      if(self%ilo == 0 .and. self%ihi == 0) then
        ! Apply across entire i range
        ilo = lbounds(2)
        ihi = ubounds(2)
      else if(self%jlo == 0 .and. self%jhi == 0) then
        ! Apply across entire j range
        jlo = lbounds(3)
        jhi = ubounds(3)
      end if

      ! For convienence, -1 just means go to the last index, ala python
      if(self%ihi == -1) ihi = ubounds(2)
      if(self%jhi == -1) jhi = ubounds(3)

      if(ilo > ihi .or. jlo > jhi) then
        write(*, '(a, 4(i0, 1x), a)') "Error: Invalid source appliation range [ilo,ihi,jlo,jhi]: [", ilo, ihi, jlo, jhi, ']'
        error stop "Error: Invalid source appliation range"
      end if

    end associate
  end function get_application_bounds

end module mod_source
