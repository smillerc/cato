module mod_error
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit

  implicit none

  integer(ik), parameter :: ALL_OK = 0
  integer(ik), parameter :: NEG_DENSITY = 1
  integer(ik), parameter :: NEG_PRESSURE = 2
  integer(ik), parameter :: NANS_FOUND = 3
end module mod_error
