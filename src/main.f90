program cato
#ifdef USE_OPENMP
  use omp_lib
#endif /* USE_OPENMP */

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_out => output_unit, std_error => error_unit
  use mod_contour_writer, only: contour_writer_t
  use mod_globals, only: print_version_stats, open_debug_files
  use mod_units, only: set_output_unit_system, io_time_label, io_time_units
  use mod_input, only: input_t
  use mod_nondimensionalization, only: set_scale_factors, t_0
  use mod_timing, only: timer_t, get_timestep
  use mod_finite_volume_schemes, only: finite_volume_scheme_t, make_fv_scheme
  use mod_fluid, only: fluid_t, new_fluid
  use mod_integrand, only: integrand_t
  use mod_grid, only: grid_t
  use mod_eos, only: set_equation_of_state

  implicit none

  character(150) :: command_line_arg
  character(50) :: input_filename
  integer(ik) :: error_code = 0
  type(input_t) :: input
  class(finite_volume_scheme_t), pointer :: fv
  class(fluid_t), pointer :: U

  type(contour_writer_t) :: contour_writer
  type(timer_t) :: timer
  real(rk) :: time = 0.0_rk
  integer(ik), dimension(3) :: bounds
  real(rk) :: delta_t = 0.0_rk
  real(rk) :: max_cs = 0.0_rk
  real(rk) :: min_dx = 0.0_rk
  real(rk) :: next_output_time = 0.0_rk
  integer(ik) :: iteration = 0
  logical :: file_exists = .false.
  real(rk) :: contour_interval_dt, max_time
  ! real(rk), dimension(:,:), allocatable :: sound_speed

  print *, 'std_error', std_error
  open(std_error, file='cato.error')

  ! ascii art for the heck of it :)
  write(std_out, '(a)')
  write(std_out, '(a)') "  ______      ___   .___________.  ______  "
  write(std_out, '(a)') " /      |    /   \  |           | /  __  \ "
  write(std_out, '(a)') "|  ,----'   /  ^  \ `---|  |----`|  |  |  |"
  write(std_out, '(a)') "|  |       /  /_\  \    |  |     |  |  |  |"
  write(std_out, '(a)') "|  `----. /  _____  \   |  |     |  `--'  |"
  write(std_out, '(a)') " \______|/__/     \__\  |__|      \______/ "
  write(std_out, '(a)')

#ifdef USE_OPENMP
  write(std_out, '(a, i0, a)') "Running with ", omp_get_max_threads(), " OpenMP threads"
#endif /* USE_OPENMP */

  call open_debug_files()

  call print_version_stats()

  call get_command_argument(1, command_line_arg)
  input_filename = trim(command_line_arg)

  inquire(file=trim(input_filename), exist=file_exists)

  if(.not. file_exists) then
    error stop 'Error in main(), CATO input (std_out.ini) file not found, exiting...'
  end if

  call input%read_from_ini(input_filename)
  call set_scale_factors(time_scale=input%reference_time, &
                         length_scale=input%reference_length, &
                         density_scale=input%reference_density)

  ! Non-dimensionalize
  contour_interval_dt = input%contour_interval_dt / t_0
  max_time = input%max_time / t_0

  call set_equation_of_state(input)
  call set_output_unit_system(input%unit_system)

  fv => make_fv_scheme(input)
  U => new_fluid(input, fv)
  contour_writer = contour_writer_t(input=input)

  if(input%restart_from_file) then
    time = fv%time / t_0
    next_output_time = time + contour_interval_dt
    iteration = fv%iteration
  else

    call contour_writer%write_contour(U, fv, time, iteration)
  end if

  print *
  write(std_out, '(a)') '--------------------------------------------'
  write(std_out, '(a)') ' Starting time loop:'
  write(std_out, '(a)') '--------------------------------------------'
  print *

  call timer%start()
  delta_t = 0.1_rk * get_timestep(cfl=input%cfl, fv=fv, fluid=U)
  call fv%set_time(time, delta_t, iteration)

  if(input%restart_from_file) then
    write(std_out, '(a, es10.3)') "Starting time:", time
    write(std_out, '(a, es10.3)') "Starting timestep:", delta_t
  end if

  bounds = lbound(U%conserved_vars)

  do while(time < max_time .and. iteration < input%max_iterations)

    write(std_out, '(2(a, es10.3), a)') 'Time =', time * io_time_units * t_0, ' '//trim(io_time_label)//', Delta t =', delta_t * t_0, ' s'

    call fv%apply_source_terms(conserved_vars=U%conserved_vars, &
                               lbounds=bounds)
    ! Integrate in time
    call U%integrate(fv, delta_t, error_code)

    if(error_code /= 0) then
      write(std_error, '(a)') 'Something went wrong in the time integration, saving to disk and exiting...'
      write(std_out, '(a)') 'Something went wrong in the time integration, saving to disk and exiting...'
      call contour_writer%write_contour(U, fv, time, iteration)
      error stop
    end if

    time = time + delta_t
    iteration = iteration + 1
    call fv%set_time(time, delta_t, iteration)
    call timer%log_time(iteration, delta_t)

    ! I/O
    if(time >= next_output_time) then
      next_output_time = next_output_time + contour_interval_dt
      write(std_out, '(a, es10.3, a)') 'Saving Contour, Next Output Time: ', &
        next_output_time * t_0 * io_time_units, ' '//trim(io_time_label)
      call contour_writer%write_contour(U, fv, time, iteration)
    end if

    delta_t = get_timestep(cfl=input%cfl, fv=fv, fluid=U)
  end do

  call timer%stop()
  deallocate(fv)
  deallocate(U)

  call timer%output_stats()
  close(std_error)
end program
