module mod_timing
  !< Define the type used for timing

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, int64
  use mod_finite_volume_schemes, only: finite_volume_scheme_t
  use mod_nondimensionalization, only: t_0
  use mod_fluid, only: fluid_t

  implicit none

  private
  public :: timer_t, get_timestep

  type :: timer_t
    private
    integer(ik) :: log_file
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
    procedure :: log_time
  end type

contains

  subroutine start(self)
    class(timer_t), intent(inout) :: self

    open(newunit=self%log_file, file='timing.csv')

    write(self%log_file, '(a)') 'iteration, elapsed_wall_time[sec], elapsed_cpu_time[sec], timestep[sec]'
    call cpu_time(self%start_cputime)

    ! Initialize the clock
    call system_clock(count_rate=self%count_rate)
    call system_clock(count_max=self%count_max)
    self%counter_rate = real(self%count_rate, rk)

    ! Start the clock
    call system_clock(count=self%count_start)
  end subroutine

  subroutine log_time(self, iteration, timestep)
    !< Keep a running log of the timings and save it to a csv file
    class(timer_t), intent(inout) :: self
    integer(ik), intent(in) :: iteration
    real(rk), intent(in) :: timestep

    integer(int64) :: count_end
    real(rk) :: elapsed_walltime
    real(rk) :: elapsed_cputime, end_cputim
    real(rk) :: end_cputime

    call cpu_time(end_cputime)
    call system_clock(count=count_end)

    elapsed_walltime = (count_end - self%count_start) / self%counter_rate
    elapsed_cputime = end_cputime - self%start_cputime

    ! header is 'iteration elapsed_wall_time[sec] elapsed_cpu_time[sec] timestep[sec]'
    write(self%log_file, '(i0, 3(", ", es14.4))') &
      iteration, elapsed_walltime, elapsed_cputime, timestep

  end subroutine log_time

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

  real(rk) function get_timestep(cfl, fv, fluid) result(delta_t)
    real(rk), intent(in) :: cfl
    class(finite_volume_scheme_t), intent(in) :: fv
    class(fluid_t), intent(in) :: fluid

    integer(ik) :: ilo, ihi, jlo, jhi

    ilo = fv%grid%ilo_cell
    ihi = fv%grid%ihi_cell
    jlo = fv%grid%jlo_cell
    jhi = fv%grid%jhi_cell

    associate(dx=>fv%grid%cell_dx, dy=>fv%grid%cell_dy)

      delta_t = minval(cfl / (((fluid%u(ilo:ihi, jlo:jhi) + fluid%cs(ilo:ihi, jlo:jhi)) / dx(ilo:ihi, jlo:jhi)) + &
                              ((fluid%v(ilo:ihi, jlo:jhi) + fluid%cs(ilo:ihi, jlo:jhi)) / dy(ilo:ihi, jlo:jhi))))
    end associate

  end function

end module mod_timing
