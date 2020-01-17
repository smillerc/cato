module mod_pressure_source
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64

  use mod_globals, only: debug_print
  use mod_source, only: source_t
  use mod_input, only: input_t
  use linear_interpolation_module, only: linear_interp_1d

  implicit none

  private
  public :: pressure_source_t, new_pressure_source

  type, extends(source_t) :: pressure_source_t
    logical :: constant_pressure = .false.
    real(rk) :: pressure_input = 0.0_rk

    ! Temporal inputs
    type(linear_interp_1d) :: temporal_pressure_input
    character(len=100) :: input_filename
  contains
    procedure :: apply_source => apply_pressure_source
    procedure :: copy => copy_pressure_source
    procedure, private :: read_pressure_input
  end type

contains

  function new_pressure_source(input) result(pressure_source)
    type(pressure_source_t), pointer :: pressure_source
    class(input_t), intent(in) :: input
    allocate(pressure_source)

    pressure_source%source_type = 'pressure'
    pressure_source%constant_pressure = input%apply_constant_source_pressure
    if(.not. pressure_source%constant_pressure) then
      pressure_source%input_filename = trim(input%pressure_source_file)
      call pressure_source%read_pressure_input()
    else
      pressure_source%pressure_input = input%constant_source_pressure_value
    end if

    if(input%pressure_source_ilo == 0 .and. &
       input%pressure_source_ihi == 0 .and. &
       input%pressure_source_jlo == 0 .and. &
       input%pressure_source_jhi == 0) then
      error stop "Error in pressure_source_t constructor; all (i,j) indicies are 0"
    end if

    pressure_source%ilo = input%pressure_source_ilo
    pressure_source%ihi = input%pressure_source_ihi
    pressure_source%jlo = input%pressure_source_jlo
    pressure_source%jhi = input%pressure_source_jhi

  end function

  subroutine read_pressure_input(self)
    !< Read the pressure input file and initialize the linear iterpolated object. This
    !< makes it easy to get the input pressure at any given time
    class(pressure_source_t), intent(inout) :: self

    real(rk), dimension(:), allocatable :: time_sec
    real(rk), dimension(:), allocatable :: pressure_barye
    character(len=300) :: line_buffer
    logical :: has_header_line
    logical :: file_exists
    integer(ik) :: input_unit, io_status, interp_status
    integer(ik) :: line, nlines

    file_exists = .false.
    has_header_line = .false.

    inquire(file=trim(self%input_filename), exist=file_exists)

    if(.not. file_exists) then
      error stop 'Error in pressure_source_t%read_pressure_input(); pressure input file not found, exiting...'
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

    allocate(time_sec(nlines))
    allocate(pressure_barye(nlines))
    time_sec = 0.0_rk
    pressure_barye = 0.0_rk

    do line = 1, nlines
      read(input_unit, *, iostat=io_status) time_sec(line), pressure_barye(line)
      if(io_status /= 0) exit
    end do
    close(input_unit)

    ! Initialize the linear interpolated data object so we can query the pressure at any time
    call self%temporal_pressure_input%initialize(time_sec, pressure_barye, interp_status)
    if(interp_status /= 0) error stop "Error initializing pressure_source_t%temporal_pressure_input"

    deallocate(time_sec)
    deallocate(pressure_barye)
  end subroutine

  subroutine copy_pressure_source(out_source, in_source)
    class(source_t), intent(in) :: in_source
    class(pressure_source_t), intent(inout) :: out_source

    select type(in_source)
    class is(pressure_source_t)
      out_source%constant_pressure = in_source%constant_pressure
      out_source%pressure_input = in_source%pressure_input
      out_source%temporal_pressure_input = in_source%temporal_pressure_input
      out_source%input_filename = in_source%input_filename
    class default
      error stop 'pressure_source_t%copy_pressure_source: unsupported in_source class'
    end select

  end subroutine

  subroutine apply_pressure_source(self, primitive_vars, time)
    class(pressure_source_t), intent(inout) :: self
    real(rk), dimension(:, 0:, 0:), intent(inout) :: primitive_vars
    real(rk), intent(in) :: time

    integer(ik) :: interp_stat
    real(rk) :: pressure

    if(self%constant_pressure) then
      pressure = self%pressure_input
    else
      call self%temporal_pressure_input%evaluate(x=time, f=pressure, istat=interp_stat)
      if(interp_stat /= 0) then
        error stop "Unable to interpolate pressure within "// &
          "pressure_source_t%apply_pressure_source()"
      end if
    end if

    if(pressure > 0.0_rk) then
      ! Apply across all i indices
      if(self%ilo == 0 .and. self%ihi == 0) then
        write(*, '(a, es14.5, 2(a, i0), a)') "Applying pressure = ", pressure, " [barye] at (i=:, j=", self%jlo, ":", self%jhi, ")"
        primitive_vars(4, :, self%jlo:self%jhi) = pressure
        ! primitive_vars(4, :, self%jlo:self%jhi) = primitive_vars(4, :, self%jlo:self%jhi) + pressure

        ! Apply across all j indicies
      else if(self%jlo == 0 .and. self%jhi == 0) then
        write(*, '(a, es14.5, 2(a, i0), a)') "Applying pressure = ", pressure, " [barye] at (i=", self%ilo, ":", self%ihi, ", j=:)"
        primitive_vars(4, self%ilo:self%ihi, :) = pressure
        ! primitive_vars(4, self%ilo:self%ihi, :) = primitive_vars(4, self%ilo:self%ihi, :) + pressure

        ! Apply within the index intervals
      else
        write(*, '(a, es14.5, 4(a, i0), a)') "Applying pressure = ", pressure, &
          " [barye] at (i=", self%ilo, ":", self%ihi, ", j=", self%jlo, ":", self%jhi, ")"
        primitive_vars(4, self%ilo:self%ihi, self%jlo:self%jhi) = pressure
        ! primitive_vars(4, self%ilo:self%ihi, self%jlo:self%jhi) = primitive_vars(4, self%ilo:self%ihi, self%jlo:self%jhi) + pressure
      end if
    end if

  end subroutine

end module mod_pressure_source
