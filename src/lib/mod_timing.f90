module mod_timing
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  implicit none

  private
  public :: timer_t

  type :: timer_t
    real(rk) :: start_time = 0.0_rk
    real(rk) :: stop_time = 0.0_rk
    real(rk) :: total_time = 0.0_rk
  contains
    procedure :: start
    procedure :: stop
    procedure :: output_stats
  end type

contains

  subroutine start(self)
    class(timer_t), intent(inout) :: self
    call cpu_time(self%start_time)
  end subroutine

  subroutine stop(self)
    class(timer_t), intent(inout) :: self
    call cpu_time(self%stop_time)
    self%total_time = self%stop_time - self%start_time
  end subroutine

  subroutine output_stats(self)
    class(timer_t), intent(in) :: self
    write(*, '(a, es10.3)') "Total elapsed time [s]:", self%total_time
    write(*, '(a, es10.3)') "Total elapsed time [m]:", self%total_time / 60.0_rk
    write(*, '(a, es10.3)') "Total elapsed time [hr]:", self%total_time / 3600.0_rk
  end subroutine

end module mod_timing
