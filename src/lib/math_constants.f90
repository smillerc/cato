module math_constants

  use iso_fortran_env, only : int32, real64

  implicit none
  
  real(real64), parameter :: pi = 4.0_real64*tan(1.0_real64)
end module math_constants