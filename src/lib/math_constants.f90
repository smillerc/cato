module math_constants

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  real(rk), parameter :: pi = 4.0_rk * atan(1.0_rk)

contains
  real(rk) elemental function rad2deg(rad) result(deg)
    real(rk), intent(in) :: rad
    deg = rad * 180.0_rk / pi
  end function

  real(rk) elemental function deg2rad(deg) result(rad)
    real(rk), intent(in) :: deg
    rad = deg * pi / 180.0_rk
  end function
end module math_constants
