module mod_energy_source
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64

  use mod_globals, only: debug_print
  use mod_source, only: source_t
  use mod_input, only: input_t
  use linear_interpolation_module, only: linear_interp_1d

  implicit none

  private
  public :: energy_source_t, new_energy_source

  type, extends(source_t) :: energy_source_t
    logical :: constant_energy = .false.
    real(rk) :: energy_input = 0.0_rk

    ! Temporal inputs
    type(linear_interp_1d) :: temporal_energy_input
    character(len=100) :: input_filename
  contains
    procedure :: apply_source => apply_energy_source
    procedure :: copy => copy_energy_source
    procedure, private :: read_energy_input
  end type

contains

  function new_energy_source(input) result(energy_source)
    type(energy_source_t), pointer :: energy_source
    class(input_t), intent(in) :: input
    allocate(energy_source)

    energy_source%source_type = 'energy'
    energy_source%constant_energy = input%apply_constant_source
    if(.not. energy_source%constant_energy) then
      energy_source%input_filename = trim(input%source_file)
      call energy_source%read_energy_input()
    else
      energy_source%energy_input = input%constant_source_value
    end if

    if(input%source_ilo == 0 .and. &
       input%source_ihi == 0 .and. &
       input%source_jlo == 0 .and. &
       input%source_jhi == 0) then
      error stop "Error in energy_source_t constructor; all (i,j) indicies are 0"
    end if

    energy_source%ilo = input%source_ilo
    energy_source%ihi = input%source_ihi
    energy_source%jlo = input%source_jlo
    energy_source%jhi = input%source_jhi

  end function

  subroutine read_energy_input(self)
    !< Read the energy input file and initialize the linear iterpolated object. This
    !< makes it easy to get the input energy at any given time
    class(energy_source_t), intent(inout) :: self

    real(rk), dimension(:), allocatable :: time_sec
    real(rk), dimension(:), allocatable :: energy_erg
    character(len=300) :: line_buffer
    logical :: has_header_line
    logical :: file_exists
    integer(ik) :: input_unit, io_status, interp_status
    integer(ik) :: line, nlines

    file_exists = .false.
    has_header_line = .false.

    inquire(file=trim(self%input_filename), exist=file_exists)

    if(.not. file_exists) then
      error stop 'Error in energy_source_t%read_energy_input(); energy input file not found, exiting...'
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
    allocate(energy_erg(nlines))
    time_sec = 0.0_rk
    energy_erg = 0.0_rk

    do line = 1, nlines
      read(input_unit, *, iostat=io_status) time_sec(line), energy_erg(line)
      if(io_status /= 0) exit
    end do
    close(input_unit)

    ! Initialize the linear interpolated data object so we can query the energy at any time
    call self%temporal_energy_input%initialize(time_sec, energy_erg, interp_status)
    if(interp_status /= 0) error stop "Error initializing energy_source_t%temporal_energy_input"

    deallocate(time_sec)
    deallocate(energy_erg)
  end subroutine

  subroutine copy_energy_source(out_source, in_source)
    class(source_t), intent(in) :: in_source
    class(energy_source_t), intent(inout) :: out_source

    select type(in_source)
    class is(energy_source_t)
      out_source%constant_energy = in_source%constant_energy
      out_source%energy_input = in_source%energy_input
      out_source%temporal_energy_input = in_source%temporal_energy_input
      out_source%input_filename = in_source%input_filename
    class default
      error stop 'energy_source_t%copy_energy_source: unsupported in_source class'
    end select

  end subroutine

  subroutine apply_energy_source(self, conserved_vars, time)
    ! //TODO: Fix this so that it deals with rho E not just E
    class(energy_source_t), intent(inout) :: self
    real(rk), dimension(:, 0:, 0:), intent(inout) :: conserved_vars
    real(rk), intent(in) :: time

    integer(ik) :: interp_stat
    real(rk) :: energy

    if(self%constant_energy) then
      energy = self%energy_input
    else
      call self%temporal_energy_input%evaluate(x=time, f=energy, istat=interp_stat)
      if(interp_stat /= 0) then
        error stop "Unable to interpolate energy within "// &
          "energy_source_t%apply_energy_source()"
      end if
    end if

    if(energy > 0.0_rk) then
      ! Apply across all i indices
      if(self%ilo == 0 .and. self%ihi == 0) then
        write(*, '(a, es14.5, 2(a, i0), a)') "Applying energy = ", energy, " [erg] at (i=:, j=", self%jlo, ":", self%jhi, ")"
        ! conserved_vars(4, :, self%jlo:self%jhi) = energy
        conserved_vars(4, :, self%jlo:self%jhi) = conserved_vars(4, :, self%jlo:self%jhi) + energy

        ! Apply across all j indicies
      else if(self%jlo == 0 .and. self%jhi == 0) then
        write(*, '(a, es14.5, 2(a, i0), a)') "Applying energy = ", energy, " [erg] at (i=", self%ilo, ":", self%ihi, ", j=:)"
        ! conserved_vars(4, self%ilo:self%ihi, :) = energy
        conserved_vars(4, self%ilo:self%ihi, :) = conserved_vars(4, self%ilo:self%ihi, :) + energy

        ! Apply within the index intervals
      else
        write(*, '(a, es14.5, 4(a, i0), a)') "Applying energy = ", energy, &
          " [erg] at (i=", self%ilo, ":", self%ihi, ", j=", self%jlo, ":", self%jhi, ")"
        ! conserved_vars(4, self%ilo:self%ihi, self%jlo:self%jhi) = energy
        conserved_vars(4, self%ilo:self%ihi, self%jlo:self%jhi) = conserved_vars(4, self%ilo:self%ihi, self%jlo:self%jhi) + energy
      end if
    end if

  end subroutine

end module mod_energy_source
