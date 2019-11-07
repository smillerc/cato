program fvleg

  use iso_fortran_env, only: int32, real64

  implicit none

  real(real64) :: t, t_stop, delta_t

  do while(t < t_stop)

    ! U = U

    t = t + delta_t
  end do

end program
