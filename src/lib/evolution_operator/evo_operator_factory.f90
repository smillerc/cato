module mod_evo_operator_factory

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_input, only: input_t
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_grid, only: grid_t
  use mod_local_evo_operator, only: local_evo_operator_t

  implicit none

  private
  public :: evo_operator_factory

contains

  function evo_operator_factory(input, grid_target, recon_operator_target, reconstructed_state_target, lbounds) result(operator)
    class(input_t), intent(in) :: input
    class(abstract_evo_operator_t), pointer :: operator

    class(grid_t), intent(in), target :: grid_target
    class(abstract_reconstruction_t), intent(in), target :: recon_operator_target
    integer(ik), dimension(5), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):, &
                        lbounds(4):, lbounds(5):), intent(in), target :: reconstructed_state_target

    character(len=32) :: evo_type = ''

    evo_type = trim(input%evolution_operator_type)
    call debug_print('Making a "'//trim(evo_type)//'" evolution operator', __FILE__, __LINE__)

    select case(trim(evo_type))
    case('cato')
      allocate(local_evo_operator_t :: operator)
      call operator%initialize(input=input, grid_target=grid_target, &
                               recon_operator_target=recon_operator_target, &
                               reconstructed_state_target=reconstructed_state_target, &
                               lbounds=lbounds)
    case default
      ! if(this_image() == 1) then
      error stop "Error in evo_operator_factory, unrecognizable evolution type"
      ! end if
    end select
    ! call debug_print('Done', __FILE__, __LINE__)
  end function

end module mod_evo_operator_factory
