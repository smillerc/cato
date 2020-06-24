module mod_finite_volume_schemes
  !< Summary: Provide the master puppeteer class that directs each separate physics package. In this case,
  !<          it's just the fluid package for now
  !< Date: 06/23/2020
  !< Author: Sam Miller
  !< Notes:

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use mod_globals, only: debug_print
  use mod_units
  use mod_input, only: input_t
  use mod_eos, only: eos
  use mod_grid, only: grid_t
  use mod_grid_factory, only: grid_factory
  use hdf5_interface, only: hdf5_file
  use mod_nondimensionalization, only: set_scale_factors
  use mod_fluid, only: fluid_t, new_fluid

  implicit none
  private
  public :: finite_volume_scheme_t, make_fv_scheme

  type :: finite_volume_scheme_t
    !< This is a puppeteer class that manages the grid, fluid solver, and any other future physics packages

    character(len=32) :: title = ''
    integer(ik) :: iteration = 0 !< iteration coount
    real(rk) :: dt = 0.0_rk      !< time step
    real(rk) :: time = 0.0_rk    !< simulation time
    class(grid_t), allocatable :: grid   !< grid topology
    class(fluid_t), allocatable :: fluid !< fluid physics

  contains
    procedure, public :: initialize
    procedure, public :: integrate
    procedure, public :: set_time
    final :: finalize
  end type

  interface make_fv_scheme
    module procedure :: constructor
  end interface

contains

  function constructor(input) result(fv)
    class(input_t), intent(in) :: input
    type(finite_volume_scheme_t), pointer :: fv

    allocate(fv)
    call fv%initialize(input)

  end function constructor

  subroutine initialize(self, input)
    !< Construct the finite volume local evolution Galerkin (fvleg) scheme
    class(finite_volume_scheme_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    ! Locals
    class(grid_t), pointer :: grid => null()
    type(hdf5_file) :: h5
    integer(ik) :: alloc_status
    integer(ik), dimension(:, :), allocatable :: ghost_layers
    character(32) :: str_buff

    alloc_status = 0

    call debug_print('Initializing finite_volume_scheme_t', __FILE__, __LINE__)

    if(input%restart_from_file) then
      call h5%initialize(filename=trim(input%restart_file), status='old', action='r')
      call h5%get('/time', self%time)
      call h5%get('/iteration', self%iteration)

      call h5%readattr('/time', 'units', str_buff)
      select case(trim(str_buff))
      case('ns')
        self%time = self%time * ns_to_s
      case('s')
        ! Do nothing, since seconds is what the code works in
      case default
        error stop "Unknown time units in .h5 file. Acceptable units are 'ns' or 's'."
      end select
      call h5%finalize()
    end if

    self%title = trim(input%title)

    grid => grid_factory(input)
    allocate(self%grid, source=grid, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%grid"
    deallocate(grid)

    ! Now that the grid is set, we can finish setting up the scale factors
    ! for non-dimensionalization. The grid sets the length scale automatically
    ! so that the smallest cell edge length is 1.
    call set_scale_factors(pressure_scale=input%reference_pressure, &
                           density_scale=input%reference_density)

    self%fluid = new_fluid(input, self%grid)
  end subroutine initialize

  subroutine finalize(self)
    !< Implementation of the class cleanup
    type(finite_volume_scheme_t), intent(inout) :: self
    call debug_print('Running finite_volume_scheme_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%grid)) deallocate(self%grid)
    if(allocated(self%fluid)) deallocate(self%fluid)
  end subroutine finalize

  subroutine integrate(self, dt)
    !< Advance the simulation forward in time
    class(finite_volume_scheme_t), intent(inout) :: self
    real(rk), intent(in) :: dt !< timestep

    self%iteration = self%iteration + 1
    self%time = self%time + dt
    call self%fluid%integrate(dt=dt, grid=self%grid)
  end subroutine integrate

  subroutine set_time(self, time, dt, iteration)
    !< Set the time statistics
    class(finite_volume_scheme_t), intent(inout) :: self
    real(rk), intent(in) :: time          !< simulation time
    real(rk), intent(in) :: dt            !< time-step
    integer(ik), intent(in) :: iteration  !< iteration

    self%time = time
    self%iteration = iteration
    self%dt = dt

    call self%fluid%set_time(time, dt, iteration)
  end subroutine set_time

end module mod_finite_volume_schemes
