module mod_energy_source
  !< Define the class for injecting energy into the domain

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64

  use mod_globals, only: debug_print
  use math_constants, only: universal_gas_const
  use mod_source, only: source_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use linear_interpolation_module, only: linear_interp_1d

  implicit none

  private
  public :: energy_source_t, new_energy_source

  type, extends(source_t) :: energy_source_t
  contains
    procedure :: apply_source => apply_energy_source
    procedure :: copy => copy_energy_source
  end type

contains

  function new_energy_source(input) result(energy_source)
    type(energy_source_t), pointer :: energy_source
    class(input_t), intent(in) :: input
    allocate(energy_source)

    energy_source%source_type = 'energy'
    energy_source%constant_source = input%apply_constant_source
    if(.not. energy_source%constant_source) then
      energy_source%input_filename = trim(input%source_file)
      call energy_source%read_input_file()
    else
      energy_source%source_input = input%constant_source_value
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

  subroutine copy_energy_source(out_source, in_source)
    class(source_t), intent(in) :: in_source
    class(energy_source_t), intent(inout) :: out_source

    select type(in_source)
    class is(energy_source_t)
      out_source%constant_source = in_source%constant_source
      out_source%source_input = in_source%source_input
      out_source%temporal_source_input = in_source%temporal_source_input
      out_source%input_filename = in_source%input_filename
    class default
      error stop 'energy_source_t%copy_energy_source: unsupported in_source class'
    end select

  end subroutine copy_energy_source

  subroutine apply_energy_source(self, conserved_vars, lbounds, time)
    !< Inject energy into the conserved variables
    class(energy_source_t), intent(inout) :: self
    integer(ik), dimension(3), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(inout) :: conserved_vars
    real(rk), intent(in) :: time

    integer(ik) :: interp_stat
    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik), dimension(4) :: index_ranges
    real(rk) :: input_energy
    real(rk), dimension(:, :), allocatable :: new_energy

    call self%set_time(time)
    input_energy = self%get_desired_source_value()
    if(input_energy > 0.0_rk) then
      index_ranges = self%get_application_bounds(lbounds=lbound(conserved_vars), &
                                                 ubounds=ubound(conserved_vars))
      ilo = index_ranges(1)
      ihi = index_ranges(2)
      jlo = index_ranges(3)
      jhi = index_ranges(4)

      allocate(new_energy(ilo:ihi, jlo:jhi))
      ! Get the energy from the conserved variable array (convert from rho E to E)
      new_energy = conserved_vars(4, ilo:ihi, jlo:jhi) / conserved_vars(1, ilo:ihi, jlo:jhi)

      ! Add the new amount of desired energy
      new_energy = new_energy + input_energy

      conserved_vars(4, ilo:ihi, jlo:jhi) = new_energy * conserved_vars(1, ilo:ihi, jlo:jhi)
      deallocate(new_energy)
    end if

  end subroutine apply_energy_source

end module mod_energy_source
