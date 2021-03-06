! MIT License
! Copyright (c) 2019 Sam Miller
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module mod_timing
  !< Summary: Provide implementations for timing and calculating the timestep
  !< Author: Sam Miller

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, int64, std_out => output_unit
  use mod_master_puppeteer, only: master_puppeteer_t
  use mod_nondimensionalization
  use mod_fluid, only: fluid_t

  implicit none

  private
  public :: timer_t

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
  endtype

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
  endsubroutine

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

  endsubroutine log_time

  subroutine stop(self)
    class(timer_t), intent(inout) :: self
    call cpu_time(self%end_cputime)
    call system_clock(count=self%count_end)
    self%elapsed_walltime = (self%count_end - self%count_start) / self%counter_rate
    self%elapsed_cputime = self%end_cputime - self%start_cputime
  endsubroutine

  subroutine output_stats(self)
    class(timer_t), intent(in) :: self

    integer(ik) :: io
    integer(ik) :: i, unit
    integer(ik), dimension(2) :: units

    if(this_image() == 1) then
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
      enddo
      close(io)
    endif
  endsubroutine
endmodule mod_timing
