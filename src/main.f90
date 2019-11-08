program fvleg

  use iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  real(rk) :: t, t_stop, delta_t

  do while(t < t_stop)

    ! U = U

    t = t + delta_t
  end do

end program
