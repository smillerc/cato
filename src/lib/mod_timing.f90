module mod_timing
  !< Define the type used for timing

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, int64, std_out => output_unit
  use mod_master_puppeteer, only: master_puppeteer_t
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
    real(rk) :: elapsed_cputime
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

    integer(ik) :: io
    integer(ik) :: i, unit
    integer(ik), dimension(2) :: units

    open(newunit=io, file='timing_summary.yaml', status='replace')
    units = [io, std_out]

    do i = 1, size(units)
      unit = units(i)
      write(unit, '(a, es10.3)') "Total elapsed wall time [s]:", self%elapsed_walltime
      write(unit, '(a, es10.3)') "Total elapsed wall time [m]:", self%elapsed_walltime / 60.0_rk
      write(unit, '(a, es10.3)') "Total elapsed wall time [hr]:", self%elapsed_walltime / 3600.0_rk
      write(unit, '(a, es10.3)') "Total elapsed CPU time [s]:", self%elapsed_cputime
      write(unit, '(a, es10.3)') "Total elapsed CPU time [m]:", self%elapsed_cputime / 60.0_rk
      write(unit, '(a, es10.3)') "Total elapsed CPU time [hr]:", self%elapsed_cputime / 3600.0_rk
    end do
    close(io)
  end subroutine

  real(rk) function get_timestep(cfl, master) result(delta_t)
    real(rk), intent(in) :: cfl
    class(master_puppeteer_t), intent(in) :: master

    integer(ik) :: ilo, ihi, jlo, jhi

    ilo = master%grid%ilo_cell
    ihi = master%grid%ihi_cell
    jlo = master%grid%jlo_cell
    jhi = master%grid%jhi_cell

    if(.not. master%fluid%prim_vars_updated) error stop "Error fluid%prim_vars_updated is .false."
    associate(dx=>master%grid%cell_dx, dy=>master%grid%cell_dy)

      delta_t = minval(cfl / &
                       (((abs(master%fluid%u(ilo:ihi, jlo:jhi)) + master%fluid%cs(ilo:ihi, jlo:jhi)) / dx(ilo:ihi, jlo:jhi)) + &
                        ((abs(master%fluid%v(ilo:ihi, jlo:jhi)) + master%fluid%cs(ilo:ihi, jlo:jhi)) / dy(ilo:ihi, jlo:jhi))))
    end associate
    ! !$omp end workshare

  end function

end module mod_timing
