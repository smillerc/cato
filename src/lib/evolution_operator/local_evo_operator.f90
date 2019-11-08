module mod_local_evo_operator
  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_evo_operator, only: abstract_evo_operator_t

  implicit none

  private
  public :: local_evo_operator_t

  type, extends(abstract_evo_operator_t) :: local_evo_operator_t
  contains
    procedure, public :: initialize => init_local_evo
  end type local_evo_operator_t

contains

  subroutine init_local_evo(self)
    class(local_evo_operator_t), intent(inout) :: self
  end subroutine init_local_evo

end module mod_local_evo_operator
