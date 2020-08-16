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

program cato
#ifdef USE_OPENMP
  use omp_lib
#endif /* USE_OPENMP */

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_out => output_unit, std_error => error_unit
  use mod_error, only: ALL_OK, error_msg
  use mod_contour_writer, only: contour_writer_t
  use mod_globals, only: print_version_stats, open_debug_files
  use mod_units, only: set_output_unit_system, io_time_label, io_time_units
  use mod_input, only: input_t
  use mod_nondimensionalization, only: set_scale_factors, t_0
  use mod_timing, only: timer_t, get_timestep
  use mod_master_puppeteer, only: master_puppeteer_t, make_master
  use mod_eos, only: set_equation_of_state

  implicit none

  character(150) :: command_line_arg
  character(50) :: input_filename
  integer(ik) :: error_code = 0
  type(input_t) :: input
  class(master_puppeteer_t), pointer :: master

  type(contour_writer_t) :: contour_writer
  type(timer_t) :: timer
  real(rk) :: time = 0.0_rk
  real(rk) :: delta_t = 0.0_rk
  real(rk) :: next_output_time = 0.0_rk
  integer(ik) :: iteration = 0
  integer(ik) :: last_io_iteration = 0
  logical :: file_exists = .false.
  real(rk) :: contour_interval_dt, max_time

  open(std_error, file='std.err')

  ! ascii art for the heck of it :)
  if(this_image() == 1) then
    write(std_out, '(a)')
    write(std_out, '(a)') "  ______      ___   .___________.  ______  "
    write(std_out, '(a)') " /      |    /   \  |           | /  __  \ "
    write(std_out, '(a)') "|  ,----'   /  ^  \ `---|  |----`|  |  |  |"
    write(std_out, '(a)') "|  |       /  /_\  \    |  |     |  |  |  |"
    write(std_out, '(a)') "|  `----. /  _____  \   |  |     |  `--'  |"
    write(std_out, '(a)') " \______|/__/     \__\  |__|      \______/ "
    write(std_out, '(a)')
    write(std_out, '(a, i0, a)') "Running on ", num_images(), " images"
  end if

#ifdef USE_OPENMP
  write(std_out, '(a, i0, a)') "Running with ", omp_get_max_threads(), " OpenMP threads"
#endif /* USE_OPENMP */

  call open_debug_files()
  call print_version_stats()
  call get_command_argument(1, command_line_arg)
  input_filename = trim(command_line_arg)

  inquire(file=trim(input_filename), exist=file_exists)

  if(.not. file_exists) then
    call error_msg(module='cato', procedure='main', &
                   message="CATO input (input.ini) file not found", file_name=__FILE__, line_number=__LINE__)
  end if

  call input%read_from_ini(input_filename)
  call input%display_config()

  call set_equation_of_state(input)
  call set_output_unit_system(input%unit_system)

  master => make_master(input)
  contour_writer = contour_writer_t(input)

  ! Non-dimensionalize
  contour_interval_dt = input%contour_interval_dt / t_0
  max_time = input%max_time / t_0

  if(input%restart_from_file) then
    time = master%time / t_0
    next_output_time = time + contour_interval_dt
    iteration = master%iteration
  else
    next_output_time = next_output_time + contour_interval_dt
    call contour_writer%write_contour(master, time, iteration)
  end if

  if(this_image() == 1) then
    print *
    write(std_out, '(a)') '--------------------------------------------'
    write(std_out, '(a)') ' Starting time loop:'
    write(std_out, '(a)') '--------------------------------------------'
    print *
  end if
  error stop

  call timer%start()
  if(input%use_constant_delta_t) then
    delta_t = input%initial_delta_t / t_0
  else
    if(input%initial_delta_t > 0.0_rk) then
      delta_t = input%initial_delta_t / t_0
      write(std_out, '(a, es10.3)') "Starting timestep (given via input, not calculated):", delta_t
    else
      delta_t = 0.1_rk * get_timestep(cfl=input%cfl, master=master)
    endif
  end if

  call master%set_time(time, delta_t, iteration)

  ! if(input%restart_from_file) then
  write(std_out, '(a, es10.3)') "Starting time:", time
  write(std_out, '(a, es10.3)') "Starting timestep:", delta_t
  ! end if

  do while(time < max_time .and. iteration < input%max_iterations)

    write(std_out, '(2(a, es10.3), a, i0)') 'Time =', time * io_time_units * t_0, &
      ' '//trim(io_time_label)//', Delta t =', delta_t * t_0, ' s, Iteration: ', iteration

    ! call master%apply_source_terms(conserved_vars=U%conserved_vars, &
    !                            lbounds=bounds)
    ! Integrate in time
    call master%integrate(dt=delta_t, error_code=error_code)

    if(error_code /= ALL_OK) then
      write(std_error, '(a)') 'Something went wrong in the time integration, saving to disk and exiting...'
      write(std_out, '(a)') 'Something went wrong in the time integration, saving to disk and exiting...'
      call contour_writer%write_contour(master, time, iteration)
      error stop 1
    end if

    time = time + delta_t
    iteration = iteration + 1
    call master%set_time(time, delta_t, iteration)
    call timer%log_time(iteration, delta_t)

    ! I/O
    if(time >= next_output_time) then
      next_output_time = next_output_time + contour_interval_dt
      write(std_out, '(a, es10.3, a)') 'Saving Contour, Next Output Time: ', &
        next_output_time * t_0 * io_time_units, ' '//trim(io_time_label)
      call contour_writer%write_contour(master, time, iteration)
      last_io_iteration = iteration
    end if

    if(.not. input%use_constant_delta_t) delta_t = get_timestep(cfl=input%cfl, master=master)
  end do

  call timer%stop()

  ! One final contour output (if necessary)
  if(iteration > last_io_iteration) call contour_writer%write_contour(master, time, iteration)

  deallocate(master)

  call timer%output_stats()
  close(std_error)
end program
