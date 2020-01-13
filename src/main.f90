program fvleg

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, output_unit
  use mod_contour_writer, only: contour_writer_t
  use mod_globals, only: print_version_stats
  use mod_input, only: input_t
  use mod_timing, only: timer_t
  use mod_finite_volume_schemes, only: finite_volume_scheme_t, make_fv_scheme
  use mod_fluid, only: fluid_t, new_fluid
  use mod_integrand, only: integrand_t
  use mod_grid, only: grid_t
  use mod_eos, only: eos
  ! use mod_grid_factory, only: grid_factory_t
  implicit none

  character(150) :: command_line_arg
  character(50) :: input_filename
  type(input_t) :: input
  ! class(fvleg_t), pointer  :: fv
  class(finite_volume_scheme_t), pointer :: fv
  class(fluid_t), pointer :: U
  ! type(fvleg_t) :: fv
  ! type(fvleg_t), target :: fv

  type(contour_writer_t) :: contour_writer
  type(timer_t) :: timer
  real(rk) :: time = 0.0_rk
  real(rk) :: delta_t
  real(rk) :: c_sound = 0.0_rk
  real(rk) :: next_output_time = 0.0_rk
  integer(ik) :: iteration = 0
  integer(ik) :: alloc_status

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

  call input%read_from_ini(input_filename)

  fv => make_fv_scheme(input)
  U => new_fluid(input, fv)

  contour_writer = contour_writer_t(input=input)

  print *
  write(*, '(a)') '--------------------------------------------'
  print *, 'Starting time loop:'
  write(*, '(a)') '--------------------------------------------'
  print *

  call timer%start()

  do while(time < input%max_time .and. iteration < input%max_iterations)

    associate(rho=>U%conserved_vars(1, :, :), &
              p=>U%conserved_vars(4, :, :), gamma=>eos%get_gamma())
      c_sound = maxval(sqrt(abs(gamma * p / rho)))
    end associate

    delta_t = min(fv%grid%min_dx, fv%grid%min_dx) * input%cfl / c_sound

    write(*, '(2(a, 1x, es10.3))') 'Time =', time, ' Delta t = ', delta_t

    ! print*, 'minval(U%conserved_vars(1,:,:))', minval(U%conserved_vars(1,:,:)), 'maxval(U%conserved_vars(1,:,:))', maxval(U%conserved_vars(1,:,:))
    call fv%apply_source_terms()

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
  end do

  call timer%stop()
  ! Write the final step
  call contour_writer%write_contour(U, fv, time, iteration)

  deallocate(fv)
  deallocate(U)

  call timer%output_stats()
end program
