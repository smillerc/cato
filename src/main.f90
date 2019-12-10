program fvleg

  use iso_fortran_env, only: ik => int32, rk => real64, output_unit
  use mod_contour_writer, only: contour_writer_t
  use mod_globals, only: print_version_stats
  use mod_input, only: input_t
  use mod_finite_volume_schemes, only: finite_volume_scheme_t
  use mod_fvleg, only: fvleg_t
  use mod_grid, only: grid_t
  ! use mod_grid_factory, only: grid_factory_t
  implicit none

  character(150) :: command_line_arg
  character(50) :: input_filename
  type(input_t) :: input
  ! class(finite_volume_scheme_t), allocatable :: fv
  type(fvleg_t) :: fv
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
  call fv%initialize(input)

  contour_writer = contour_writer_t(input=input)

  delta_t = input%initial_delta_t
  do while(time < input%max_time)
    write(*, '(2(a, 1x, g0.4))') 'Time:', time, ' Delta t:', delta_t

    call fv%integrate(delta_t)

    if(time >= next_output_time) then
      call contour_writer%write_contour(fv, time, iteration)
      next_output_time = next_output_time + input%contour_interval_dt
    end if

    time = time + delta_t
    iteration = iteration + 1
  end do

end program
