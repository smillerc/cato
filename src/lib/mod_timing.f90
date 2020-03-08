module mod_timing
  !< Define the type used for timing

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, int64
  implicit none

  private
  public :: timer_t

  type :: timer_t
    private
    integer(int64) :: count_start = 0
    integer(int64) :: count_end = 0
    integer(int64) :: count_rate = 0
    integer(int64) :: count_max = 0
    real(rk) :: counter_rate = 0.0_rk
    real(rk) :: start_cputime = 0.0_rk
    real(rk) :: end_cputime = 0.0_rk

    real(rk), public :: elapsed_cputime = 0.0_rk
    real(rk), public :: elapsed_walltime = 0.0_rk
  contains
    procedure :: start
    procedure :: stop
    procedure :: output_stats
  end type

contains

  subroutine start(self)
    class(timer_t), intent(inout) :: self
    call cpu_time(self%start_cputime)

    ! Initialize the clock
    call system_clock(count_rate=self%count_rate)
    call system_clock(count_max=self%count_max)
    self%counter_rate = real(self%count_rate, rk)

    ! Start the clock
    call system_clock(count=self%count_start)
  end subroutine

  subroutine stop(self)
    class(timer_t), intent(inout) :: self
    call cpu_time(self%end_cputime)
    call system_clock(count=self%count_end)
    self%elapsed_walltime = (self%count_end - self%count_start) / self%counter_rate
    self%elapsed_cputime = self%end_cputime - self%start_cputime
  end subroutine

  subroutine output_stats(self)
    class(timer_t), intent(in) :: self
    write(*, '(a, es10.3)') "Total elapsed wall time [s]:", self%elapsed_walltime
    write(*, '(a, es10.3)') "Total elapsed wall time [m]:", self%elapsed_walltime / 60.0_rk
    write(*, '(a, es10.3)') "Total elapsed wall time [hr]:", self%elapsed_walltime / 3600.0_rk
    write(*, '(a, es10.3)') "Total elapsed CPU time [s]:", self%elapsed_cputime
    write(*, '(a, es10.3)') "Total elapsed CPU time [m]:", self%elapsed_cputime / 60.0_rk
    write(*, '(a, es10.3)') "Total elapsed CPU time [hr]:", self%elapsed_cputime / 3600.0_rk
  end subroutine

end module mod_timing
