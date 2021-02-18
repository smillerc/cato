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

module mod_master_puppeteer
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
  use mod_grid_block, only: grid_block_t
  ! use mod_grid_block_1d, only: grid_block_1d_t
  use mod_grid_block_2d, only: grid_block_2d_t
  ! use mod_grid_block_3d, only: grid_block_3d_t
  use mod_grid_factory, only: grid_factory
  use mod_field, only: field_2d_t
  use hdf5_interface, only: hdf5_file
  use mod_nondimensionalization
  use mod_fluid, only: fluid_t, new_fluid
  use mod_source, only: source_t, new_source

  implicit none
  private
  public :: master_puppeteer_t, make_master

  type :: master_puppeteer_t
    !< This is a puppeteer class that manages the grid, fluid solver, and any other future physics packages

    character(len=32) :: title = ''
    integer(ik) :: iteration = 0 !< iteration coount
    real(rk) :: dt = 0.0_rk      !< time step
    real(rk) :: time = 0.0_rk    !< simulation time
    class(grid_block_t), allocatable :: grid   !< grid topology
    class(fluid_t), allocatable :: fluid !< fluid physics
    logical :: is_1d = .false.
    logical :: is_2d = .false.
    logical :: is_3d = .false.
    class(source_t), allocatable :: source_term !< source_term, i.e. gravity, laser, etc..

    logical :: do_source_terms = .false.

  contains
    procedure, public :: initialize
    procedure, public :: integrate
    procedure, public :: set_time
    procedure, public :: get_timestep
    final :: finalize
  endtype

  interface make_master
    module procedure :: constructor
  endinterface

contains

  function constructor(input) result(master)
    class(input_t), intent(in) :: input
    type(master_puppeteer_t), pointer :: master

    allocate(master)
    call master%initialize(input)

  endfunction constructor

  subroutine initialize(self, input)
    !< Construct the puppeteer class
    class(master_puppeteer_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    ! Locals
    class(grid_block_t), pointer :: grid => null()
    class(fluid_t), pointer :: fluid => null()
    class(source_t), pointer :: source_term => null()
    type(hdf5_file) :: h5
    integer(ik) :: alloc_status
    character(32) :: str_buff

    alloc_status = 0

    call debug_print('Initializing master_puppeteer_t', __FILE__, __LINE__)

    self%do_source_terms = input%enable_source_terms

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
      endselect
      call h5%finalize()
    endif

    self%title = trim(input%title)

    grid => grid_factory(input)
    allocate(self%grid, source=grid, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate master_puppeteer_t%grid"
    deallocate(grid)

    self%time = self%time * t_to_nondim

    fluid => new_fluid(input, self%grid, time=self%time)
    allocate(self%fluid, source=fluid, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate master_puppeteer_t%fluid"
    deallocate(fluid)

    if(input%enable_source_terms) then
      source_term => new_source(input, self%grid, time=self%time)
      allocate(self%source_term, source=source_term, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate master_puppeteer_t%source_term"
      deallocate(source_term)
    endif

    call debug_print('Done initializing master_puppeteer_t', __FILE__, __LINE__)
  endsubroutine initialize

  subroutine finalize(self)
    !< Implementation of the class cleanup
    type(master_puppeteer_t), intent(inout) :: self
    call debug_print('Running master_puppeteer_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%grid)) deallocate(self%grid)
    if(allocated(self%fluid)) deallocate(self%fluid)
  endsubroutine finalize

  subroutine integrate(self, dt, error_code)
    !< Advance the simulation forward in time
    class(master_puppeteer_t), intent(inout) :: self
    integer(ik), intent(out) :: error_code
    real(rk), optional, intent(in) :: dt

    class(field_2d_t), allocatable, save :: d_dt_source

    real(rk) :: min_dt

    if(present(dt)) then
      self%dt = dt
    else
      self%dt = 0.0_rk
    endif

    self%iteration = self%iteration + 1
    min_dt = self%get_timestep()

    if(self%dt > min_dt) then
      if(this_image() == 1) then
        write(*, '(a,i0,a,2(es16.6,a))') "Warning, the input dt (", self%dt, &
          ") is larger than the max allowable dt (", min_dt, &
          ") based on the CFL condition"
      endif
    else
      self%dt = min_dt
    endif

    self%time = self%time + self%dt

    if(self%do_source_terms) then
      d_dt_source = self%source_term%integrate(dt=self%dt, density=self%fluid%rho)
    endif

    select type(grid => self%grid)
    class is(grid_block_2d_t)

      call self%fluid%integrate(dt=self%dt, source_term=d_dt_source, &
                                grid=self%grid, error_code=error_code)
    endselect

  endsubroutine integrate

  real(rk) function get_timestep(self) result(dt)
    !< Get the maximum allowable timestep from each physics package
    class(master_puppeteer_t), intent(inout) :: self

    select type(grid => self%grid)
    class is(grid_block_2d_t)
      dt = self%fluid%get_timestep(self%grid)
    endselect
  endfunction

  subroutine set_time(self, time, iteration)
    !< Set the time statistics
    class(master_puppeteer_t), intent(inout) :: self
    real(rk), intent(in) :: time          !< simulation time
    integer(ik), intent(in) :: iteration  !< iteration

    self%time = time
    self%iteration = iteration
    call self%fluid%set_time(time, iteration)
  endsubroutine set_time

endmodule mod_master_puppeteer
