module mod_floating_point_utils
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  private
  public :: near_zero, equal

  real(rk), parameter :: TINY_RK = tiny(1.0_rk)
  real(rk), parameter :: TINY_FACTOR = 5.0_rk

contains

  elemental function near_zero(number, epsilon) result(return_value)
    !< Test if a floading point is near 0. Use this instead of doing a
    !< `if (number == 0.0_rk)` check. Taken from "Modern Fortran: Style and Usage"
    !<  by N. Clerman and W. Spector in Chapter 13, Section 2, Rule 171

    real(rk), intent(in)  :: number
    real(rk), intent(in), optional :: epsilon
    logical :: return_value
    real(kind(epsilon)) :: local_epsilon

    local_epsilon = TINY_FACTOR * TINY_RK
    if(present(epsilon)) then
      if(abs(epsilon) >= TINY_RK) &
        local_epsilon = abs(epsilon)
    end if

    return_value = abs(number) < local_epsilon
  end function near_zero

  elemental function equal(lhs, rhs, epsilon) result(return_value)
    !< Test if two floating point numbers are equal

    real(rk), intent(in)  :: lhs
    real(rk), intent(in)  :: rhs
    real(rk), intent(in), optional :: epsilon
    logical :: return_value
    real(kind(epsilon)) :: local_epsilon

    local_epsilon = TINY_FACTOR * TINY_RK
    if(present(epsilon)) then
      if(abs(epsilon) >= TINY_RK) &
        local_epsilon = abs(epsilon)
    end if

    return_value = abs(lhs - rhs) < local_epsilon
  end function equal

end module mod_floating_point_utils
