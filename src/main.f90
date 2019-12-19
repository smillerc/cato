program fvleg

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, output_unit
  use mod_contour_writer, only: contour_writer_t
  use mod_globals, only: print_version_stats
  use mod_input, only: input_t
  use mod_finite_volume_schemes, only: finite_volume_scheme_t
  use mod_fvleg, only: fvleg_t, new_fvleg
  use mod_grid, only: grid_t
  ! use mod_grid_factory, only: grid_factory_t
  implicit none

  character(150) :: command_line_arg
  character(50) :: input_filename
  type(input_t) :: input
  ! class(fvleg_t), pointer  :: fv
  class(finite_volume_scheme_t), pointer :: fv
  ! type(fvleg_t) :: fv
  ! type(fvleg_t), target :: fv

  type(contour_writer_t) :: contour_writer
  real(rk) :: time = 0.0_rk
  real(rk) :: delta_t
  real(rk) :: next_output_time = 0.0_rk
  integer(ik) :: iteration = 0

  ! ascii art for the heck of it :)
  write(output_unit, '(a)')
  write(output_unit, '(a)') "_______________   ____.____     ___________ ________            ________  ________"
  write(output_unit, '(a)') "\_   _____/\   \ /   /|    |    \_   _____//  _____/            \_____  \ \______ \"
  write(output_unit, '(a)') " |    __)   \   Y   / |    |     |    __)_/   \  ___    ______   /  ____/  |    |  \"
  write(output_unit, '(a)') " |     \     \     /  |    |___  |        \    \_\  \  /_____/  /       \  |    `   \"
  write(output_unit, '(a)') " \___  /      \___/   |_______ \/_______  /\______  /           \_______ \/_______  /"
  write(output_unit, '(a)') "     \/                       \/        \/        \/                    \/        \/"
  write(output_unit, '(a)')

  call print_version_stats()

  call get_command_argument(1, command_line_arg)
  input_filename = trim(command_line_arg)

  ! allocate(input_t :: input)
  call input%read_from_ini(input_filename)

  fv => new_fvleg(input)

  contour_writer = contour_writer_t(input=input)

  delta_t = input%initial_delta_t

  print *, 'Starting time loop:'
  do while(time < input%max_time .and. iteration < input%max_iterations)
    print *
    write(*, '(a)') '--------------------------------------------'
    write(*, '(2(a, 1x, es10.3))') 'Time:', time, ' Delta t:', delta_t
    write(*, '(a)') '--------------------------------------------'
    print *
    call fv%integrate(delta_t)

    if(time >= next_output_time) then
      next_output_time = next_output_time + input%contour_interval_dt
      write(*, '(a, es10.3))') 'Saving Contour, Next Output Time: ', next_output_time
      call contour_writer%write_contour(fv, time, iteration)
    end if
    time = time + delta_t
    iteration = iteration + 1
  end do

end program
