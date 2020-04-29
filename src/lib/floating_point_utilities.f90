module mod_floating_point_utils
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  private
  public :: near_zero, equal, EPS

  real(rk), parameter :: TINY_RK = tiny(1.0_rk)
  real(rk), parameter :: TINY_FACTOR = 5.0_rk
  real(rk), parameter :: EPS = TINY_FACTOR * TINY_RK

  interface near_zero
    module procedure :: near_zero_base, near_zero_custom
  end interface

  interface equal
    module procedure :: equal_base, equal_custom
  end interface

contains

  elemental function near_zero_base(number) result(return_value)
    !< Test if a floading point is near 0. Use this instead of doing a
    !< `if (number == 0.0_rk)` check. Taken from "Modern Fortran: Style and Usage"
    !<  by N. Clerman and W. Spector in Chapter 13, Section 2, Rule 171
    real(rk), intent(in) :: number
    logical :: return_value

    return_value = abs(number) < EPS
  end function near_zero_base

  elemental function near_zero_custom(number, epsilon) result(return_value)
    !< Test if a floading point is near 0. Use this instead of doing a
    !< `if (number == 0.0_rk)` check. Taken from "Modern Fortran: Style and Usage"
    !<  by N. Clerman and W. Spector in Chapter 13, Section 2, Rule 171
    real(rk), intent(in) :: number
    real(rk), intent(in) :: epsilon
    logical :: return_value
    real(rk) :: local_epsilon

    if(abs(epsilon) >= TINY_RK) local_epsilon = abs(epsilon)
    return_value = abs(number) < local_epsilon

  end function near_zero_custom

  elemental function equal_base(lhs, rhs) result(return_value)
    !< Test if two floating point numbers are equal
    real(rk), intent(in)  :: lhs
    real(rk), intent(in)  :: rhs
    logical :: return_value

    return_value = abs(lhs - rhs) < EPS

  end function equal_base

  elemental function equal_custom(lhs, rhs, epsilon) result(return_value)
    !< Test if two floating point numbers are equal

    real(rk), intent(in)  :: lhs
    real(rk), intent(in)  :: rhs
    real(rk), intent(in) :: epsilon
    logical :: return_value
    real(rk) :: local_epsilon

    if(abs(epsilon) >= TINY_RK) local_epsilon = abs(epsilon)
    return_value = abs(lhs - rhs) < local_epsilon

  end function equal_custom

end module mod_floating_point_utils
