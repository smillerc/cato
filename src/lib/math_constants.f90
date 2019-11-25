module math_constants

  use iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  real(rk), parameter :: pi = 4.0_rk * atan(1.0_rk)
end module math_constants
