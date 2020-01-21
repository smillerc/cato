program cato

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, output_unit
  use mod_contour_writer, only: contour_writer_t
  use mod_globals, only: print_version_stats
  use mod_input, only: input_t
  use mod_timing, only: timer_t
  use mod_finite_volume_schemes, only: finite_volume_scheme_t, make_fv_scheme
  use mod_fluid, only: fluid_t, new_fluid
  use mod_integrand, only: integrand_t
  use mod_grid, only: grid_t
  use mod_eos, only: set_equation_of_state, eos
  ! use mod_grid_factory, only: grid_factory_t
  implicit none

  character(150) :: command_line_arg
  character(50) :: input_filename
  type(input_t) :: input
  class(finite_volume_scheme_t), pointer :: fv
  class(fluid_t), pointer :: U

  type(contour_writer_t) :: contour_writer
  type(timer_t) :: timer
  real(rk) :: time = 0.0_rk
  real(rk) :: delta_t
  real(rk) :: c_sound = 0.0_rk
  real(rk) :: next_output_time = 0.0_rk
  integer(ik) :: iteration = 0
  logical :: file_exists = .false.

  ! ascii art for the heck of it :)
  write(output_unit, '(a)')
  write(output_unit, '(a)') "  ______      ___   .___________.  ______  "
  write(output_unit, '(a)') " /      |    /   \  |           | /  __  \ "
  write(output_unit, '(a)') "|  ,----'   /  ^  \ `---|  |----`|  |  |  |"
  write(output_unit, '(a)') "|  |       /  /_\  \    |  |     |  |  |  |"
  write(output_unit, '(a)') "|  `----. /  _____  \   |  |     |  `--'  |"
  write(output_unit, '(a)') " \______|/__/     \__\  |__|      \______/ "
  write(output_unit, '(a)')

  call print_version_stats()

  call get_command_argument(1, command_line_arg)
  input_filename = trim(command_line_arg)

  inquire(file=trim(input_filename), exist=file_exists)

  if(.not. file_exists) then
    error stop 'Error in main(), CATO input (*.ini) file not found, exiting...'
  end if

  call input%read_from_ini(input_filename)

  fv => make_fv_scheme(input)
  U => new_fluid(input, fv)
  call set_equation_of_state(input)

  contour_writer = contour_writer_t(input=input)

  print *
  write(*, '(a)') '--------------------------------------------'
  print *, 'Starting time loop:'
  write(*, '(a)') '--------------------------------------------'
  print *

  call timer%start()

  do while(time < input%max_time .and. iteration < input%max_iterations)

    ! Set the timestep via the CFL constraint
    c_sound = maxval(U%get_sound_speed())
    delta_t = min(fv%grid%min_dx, fv%grid%min_dx) * input%cfl / c_sound

    write(*, '(2(a, 1x, es10.3))') 'Time =', time, ', delta t = ', delta_t

    call fv%apply_source_terms(U%conserved_vars, lbound(U%conserved_vars))

    ! First put conserved vars in ghost layers
    call fv%apply_conserved_vars_bc(U%conserved_vars, lbound(U%conserved_vars))

    ! Reference state is a neighbor average (which is why edges needed ghost U vars)
    call fv%calculate_reference_state(U%conserved_vars, lbound(U%conserved_vars))

    ! Now we can reconstruct the entire domain
    call fv%reconstruct(U%conserved_vars, lbound(U%conserved_vars))

    ! Apply the reconstructed state to the ghost layers
    call fv%apply_reconstructed_state_bc()
    call fv%apply_cell_gradient_bc()

    ! Evolve U at each edge
    call fv%evolve_domain()

    if(time >= next_output_time) then
      next_output_time = next_output_time + input%contour_interval_dt
      write(*, '(a, es10.3)') 'Saving Contour, Next Output Time: ', next_output_time
      call contour_writer%write_contour(U, fv, time, iteration)
    end if

    ! Integrate in time
    call U%integrate(fv, delta_t)

    time = time + delta_t
    iteration = iteration + 1
    call fv%set_time(time)
  end do

  call timer%stop()
  ! Write the final step
  call contour_writer%write_contour(U, fv, time, iteration)

  deallocate(fv)
  deallocate(U)

  call timer%output_stats()
end program
