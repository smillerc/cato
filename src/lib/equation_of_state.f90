module mod_eos

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t
  implicit none

  private
  public :: set_equation_of_state, eos

  type :: eos_t
    private
    real(rk) :: gamma = 5.0_rk / 3.0_rk
  contains
    procedure, public :: get_gamma
    procedure, public :: set_gamma
  end type eos_t

  interface eos_t
    module procedure :: constructor
  end interface

  type(eos_t) :: eos !< singleton equation of state object

contains
  type(eos_t) function constructor(input) result(eq_of_state)
    !< Constructor for the equation of state type
    class(input_t), intent(in) :: input
    call eq_of_state%set_gamma(input%polytropic_index)
  end function

  subroutine set_gamma(self, gamma)
    class(eos_t), intent(inout) :: self
    real(rk), intent(in) :: gamma
    self%gamma = gamma
  end subroutine

  real(rk) pure function get_gamma(self) result(gamma)
    class(eos_t), intent(in) :: self
    gamma = self%gamma
  end function

  subroutine set_equation_of_state(input)
    class(input_t), intent(in) :: input
    eos = constructor(input)
  end subroutine

end module mod_eos
