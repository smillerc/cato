module mod_units
  !< Convert from CGS to some other unit

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  ! Time conversion
  real(rk), parameter :: pico_seconds = 1e-12_rk !< sec to ps conversion
  real(rk), parameter :: nano_seconds = 1e-9_rk !< sec to ps conversion

  ! Vel
  real(rk), parameter :: km_per_sec = 1e-5_rk !< cm/s to km/s conversion
  real(rk), parameter :: um_per_ns = 1e-5_rk !< cm/s to um/ns conversion

  ! Distance
  real(rk), parameter :: microns = 1e4_rk !< cm to microns conversion

contains

end module mod_units
