module mod_abstract_evo_operator
  !< Provide a base type for the evolution operators

  use iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  private
  public :: abstract_evo_operator_t

  type, abstract :: abstract_evo_operator_t
    character(:), allocatable :: name
    integer(ik) :: error_code
    character(:), allocatable :: error_message
  contains
    procedure(initialize), deferred, public :: initialize
  end type abstract_evo_operator_t

  abstract interface
    subroutine initialize(self)
      import :: abstract_evo_operator_t
      class(abstract_evo_operator_t), intent(inout) :: self
    end subroutine initialize
  end interface
end module mod_abstract_evo_operator
