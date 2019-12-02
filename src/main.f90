program fvleg

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t
  use mod_finite_volume_schemes, only: finite_volume_scheme_t
  use mod_grid, only: grid_t
  use mod_grid_factory, only: grid_factory_t
  implicit none

  character(150) :: command_line_arg
  character(50) :: input_filename
  class(input_t), allocatable :: input
  class(grid_factory_t), allocatable :: grid_factory
  class(grid_t), pointer :: grid => null()
  class(finite_volume_scheme_t), allocatable :: fv
  real(rk) :: t = 0.0_rk, delta_t

  call get_command_argument(1, command_line_arg)
  input_filename = trim(command_line_arg)

  allocate(input_t :: input)
  call input%read_from_ini(input_filename)

  grid_factory = grid_factory_t()
  grid => grid_factory%create_grid(input)

  delta_t = input%initial_delta_t
  do while(t < input%max_time)
    write(*, '(2(a, 1x, g0.4))') 'Time:', t, ' Delta t:', delta_t
    ! U = U
    ! call fv%integrate(delta_t)

    t = t + delta_t

  end do

end program
