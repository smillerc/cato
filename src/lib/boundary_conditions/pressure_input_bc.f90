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

module mod_pressure_input_bc
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use mod_globals, only: debug_print
  use mod_grid, only: grid_t
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_nondimensionalization, only: p_0, rho_0, t_0
  use mod_input, only: input_t
  use mod_inlet_outlet, only: subsonic_inlet, subsonic_outlet, supersonic_inlet, supersonic_outlet
  use linear_interpolation_module, only: linear_interp_1d
  use mod_eos, only: eos
  implicit none

  private
  public :: pressure_input_bc_t, pressure_input_bc_constructor

  type, extends(boundary_condition_t) :: pressure_input_bc_t
    logical :: constant_pressure = .false.
    logical :: constant_density = .false.
    real(rk) :: pressure_input = 0.0_rk
    real(rk) :: temperature_input = 273.0_rk ! K
    real(rk) :: density_input = 5e-3_rk !< Boundary density g/cc
    real(rk) :: scale_factor = 1.0_rk !< scaling factor (e.g. 1x, 2x) to scale the input up/down

    ! Temporal inputs
    type(linear_interp_1d) :: temporal_pressure_input
    type(linear_interp_1d) :: temporal_density_input
    character(len=100) :: input_filename

    real(rk), dimension(:), allocatable :: edge_rho !< boundary (ghost/edge/etc) density
    real(rk), dimension(:), allocatable :: edge_u   !< boundary (ghost/edge/etc) x-velocity
    real(rk), dimension(:), allocatable :: edge_v   !< boundary (ghost/edge/etc) y-velocity
    real(rk), dimension(:), allocatable :: edge_p   !< boundary (ghost/edge/etc) pressure

  contains
    procedure, private :: read_pressure_input
    procedure, private :: get_desired_pressure
    procedure, private :: get_desired_density
    procedure, public :: apply => apply_pressure_input_primitive_var_bc
    final :: finalize
  end type
contains

  function pressure_input_bc_constructor(location, input, grid, time) result(bc)
    type(pressure_input_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input
    class(grid_t), intent(in) :: grid
    real(rk), intent(in) :: time

    allocate(bc)
    bc%name = 'pressure_input'
    bc%location = location
    call bc%set_indices(grid)

    bc%time = time
    bc%density_input = input%bc_density / rho_0
    bc%constant_pressure = input%apply_constant_bc_pressure
    bc%scale_factor = input%bc_pressure_scale_factor
    if(.not. bc%constant_pressure) then
      bc%input_filename = trim(input%bc_pressure_input_file)
      call bc%read_pressure_input()
    else
      bc%pressure_input = input%constant_bc_pressure_value / p_0 ! non-dimensionalize
    end if

  end function pressure_input_bc_constructor

  subroutine read_pressure_input(self)
    !< Read the pressure input file and initialize the linear iterpolated object. This
    !< makes it easy to get the input pressure at any given time
    class(pressure_input_bc_t), intent(inout) :: self

    real(rk), dimension(:), allocatable :: time_sec       !< User input for time in seconds
    real(rk), dimension(:), allocatable :: pressure_barye !< User input for pressure in barye
    real(rk), dimension(:), allocatable :: density_gcc    !< User input for density in g/cc
    character(len=300) :: line_buffer
    logical :: has_header_line
    logical :: file_exists
    integer(ik) :: input_unit, io_status, interp_status
    integer(ik) :: line, nlines

    file_exists = .false.
    has_header_line = .false.

    inquire(file=trim(self%input_filename), exist=file_exists)

    if(.not. file_exists) then
      error stop 'Error in pressure_input_bc_t%read_pressure_input(); pressure input file not found, exiting...'
    end if

    ! Open and determine how many lines there are
    open(newunit=input_unit, file=trim(self%input_filename))
    nlines = 0
    do
      read(input_unit, *, iostat=io_status) line_buffer
      if(nlines == 0) then
        line_buffer = adjustl(line_buffer)
        if(line_buffer(1:1) == '#') has_header_line = .true.
      end if

      if(io_status /= 0) exit
      nlines = nlines + 1
    end do
    close(input_unit)

    ! Re-open and now allocate the arrays to the proper size based on the # of lines in the file
    open(newunit=input_unit, file=trim(self%input_filename))
    if(has_header_line) then
      read(input_unit, *, iostat=io_status) ! skip the first line
      nlines = nlines - 1
    end if

    allocate(time_sec(nlines))
    allocate(pressure_barye(nlines))
    allocate(density_gcc(nlines))
    time_sec = 0.0_rk
    pressure_barye = 0.0_rk
    density_gcc = 0.0_rk

    do line = 1, nlines
      read(input_unit, *, iostat=io_status) time_sec(line), pressure_barye(line), density_gcc(line)
      if(io_status /= 0) exit
    end do
    close(input_unit)

    ! Apply scaling and non-dimensionalization
    time_sec = time_sec / t_0
    pressure_barye = (pressure_barye / p_0) * self%scale_factor
    density_gcc = (density_gcc / rho_0) * self%scale_factor

    self%max_time = maxval(time_sec)

    ! Initialize the linear interpolated data object so we can query the pressure at any time
    call self%temporal_pressure_input%initialize(time_sec, pressure_barye, interp_status)
    if(interp_status /= 0) then
      write(*, *) time_sec
      write(*, *) pressure_barye
      write(*, *) density_gcc
      write(*, *) interp_status
      error stop "Error initializing pressure_input_bc_t%temporal_pressure_input"
    end if

    call self%temporal_density_input%initialize(time_sec, density_gcc, interp_status)
    if(interp_status /= 0) then
      write(*, *) time_sec
      write(*, *) pressure_barye
      write(*, *) density_gcc
      write(*, *) interp_status
      error stop "Error initializing pressure_input_bc_t%temporal_density_input"
    end if

    deallocate(time_sec)
    deallocate(pressure_barye)
    deallocate(density_gcc)
  end subroutine read_pressure_input

  real(rk) function get_desired_pressure(self) result(desired_pressure)
    class(pressure_input_bc_t), intent(inout) :: self
    integer(ik) :: interp_stat

    if(self%constant_pressure) then
      desired_pressure = self%pressure_input
    else
      call self%temporal_pressure_input%evaluate(x=self%get_time(), f=desired_pressure, istat=interp_stat)
      if(interp_stat /= 0) then
        error stop "Unable to interpolate pressure within pressure_input_bc_t%get_desired_pressure()"
      end if
    end if

    if(desired_pressure < tiny(1.0_rk)) then
      error stop "Error in pressure_input_bc_t%get_desired_pressure: desired_pressure < tiny(1.0_rk)"
    end if

  end function get_desired_pressure

  real(rk) function get_desired_density(self) result(desired_density)
    class(pressure_input_bc_t), intent(inout) :: self
    integer(ik) :: interp_stat

    if(self%constant_pressure) then
      desired_density = self%density_input
    else
      call self%temporal_density_input%evaluate(x=self%get_time(), f=desired_density, istat=interp_stat)

      if(interp_stat /= 0) then
        error stop "Unable to interpolate density within pressure_input_bc_t%get_desired_density()"
      end if
    end if

    if(desired_density < tiny(1.0_rk)) then
      error stop "Error in pressure_input_bc_t%get_desired_density: desired_density < tiny(1.0_rk)"
    end if
  end function get_desired_density

  subroutine apply_pressure_input_primitive_var_bc(self, rho, u, v, p, lbounds)
    !< Apply pressure_input boundary conditions to the conserved state vector field
    class(pressure_input_bc_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: rho
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: u
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: v
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: p

    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    logical :: inflow
    logical :: outflow

    real(rk) :: desired_boundary_density  !< user-specified boundary value for density
    real(rk) :: desired_boundary_pressure !< user-specified boundary value for pressure
    real(rk) :: boundary_density, gamma, boundary_cs
    real(rk) :: mach_u, mach_v, mach, cs
    real(rk), dimension(:), allocatable :: edge_pressure
    real(rk), dimension(:), allocatable :: domain_rho
    real(rk), dimension(:), allocatable :: domain_u
    real(rk), dimension(:), allocatable :: domain_v
    real(rk), dimension(:), allocatable :: domain_p
    real(rk), dimension(4) :: boundary_prim_vars, domain_prim_vars
    real(rk), dimension(2) :: boundary_norm

    gamma = eos%get_gamma()
    inflow = .false.
    outflow = .false.

    desired_boundary_pressure = self%get_desired_pressure()
    desired_boundary_density = self%get_desired_density()

    if(allocated(self%edge_rho)) deallocate(self%edge_rho)
    if(allocated(self%edge_u)) deallocate(self%edge_u)
    if(allocated(self%edge_v)) deallocate(self%edge_v)
    if(allocated(self%edge_p)) deallocate(self%edge_p)

    associate(left => self%ilo, right => self%ihi, bottom => self%jlo, top => self%jhi, &
              left_ghost => minval(self%ilo_ghost), &
              right_ghost => maxval(self%ihi_ghost), &
              bottom_ghost => minval(self%jlo_ghost), &
              top_ghost => maxval(self%jhi_ghost), &
              n => self%n_ghost_layers)

      select case(self%location)
      case('+x', '-x')
        allocate(self%edge_rho(bottom_ghost:top_ghost))
        allocate(self%edge_u(bottom_ghost:top_ghost))
        allocate(self%edge_v(bottom_ghost:top_ghost))
        allocate(self%edge_p(bottom_ghost:top_ghost))
        allocate(domain_rho(bottom_ghost:top_ghost))
        allocate(domain_u(bottom_ghost:top_ghost))
        allocate(domain_v(bottom_ghost:top_ghost))
        allocate(domain_p(bottom_ghost:top_ghost))
      case('+y', '-y')
        allocate(self%edge_rho(left_ghost:right_ghost))
        allocate(self%edge_u(left_ghost:right_ghost))
        allocate(self%edge_v(left_ghost:right_ghost))
        allocate(self%edge_p(left_ghost:right_ghost))
        allocate(domain_rho(left_ghost:right_ghost))
        allocate(domain_u(left_ghost:right_ghost))
        allocate(domain_v(left_ghost:right_ghost))
        allocate(domain_p(left_ghost:right_ghost))
      end select

      self%edge_rho = 0.0_rk
      self%edge_u = 0.0_rk
      self%edge_v = 0.0_rk
      self%edge_p = 0.0_rk

      domain_rho = 0.0_rk
      domain_u = 0.0_rk
      domain_v = 0.0_rk
      domain_p = 0.0_rk

      select case(self%location)
      case('+x')
        call debug_print('Running pressure_input_bc_t%apply_pressure_input_primitive_var_bc() +x', __FILE__, __LINE__)

        domain_rho = rho(right, bottom_ghost:top_ghost)
        domain_u = u(right, bottom_ghost:top_ghost)
        domain_v = 0.0_rk
        domain_p = p(right, bottom_ghost:top_ghost)

        boundary_norm = [1.0_rk, 0.0_rk] ! outward normal vector at +x

        do j = bottom, top
          associate(rho_d => domain_rho(j), u_d => domain_u(j), &
                    v_d => domain_v(j), p_d => domain_p(j))
            domain_prim_vars = [rho_d, u_d, v_d, p_d]

            call eos%sound_speed(p=p_d, rho=rho_d, cs=cs)
            call eos%sound_speed(p=desired_boundary_pressure, rho=desired_boundary_density, cs=boundary_cs)
            mach_u = u_d / cs
            ! print*, 'domain mach', mach_u, 'ghost mach',  u_d / boundary_cs
            ! print*, "maxval(domain_p) > desired_boundary_pressure ", maxval(domain_p(bottom:top)) > desired_boundary_pressure, maxval(domain_p(bottom:top)), desired_boundary_pressure

            if(mach_u > 0.0_rk) then ! outlet
              if(mach_u <= 1.0_rk .or. maxval(domain_p(bottom:top)) < desired_boundary_pressure) then
                ! print*, 'sub outlet'
                ! error stop
                boundary_prim_vars = subsonic_outlet(domain_prim_vars=domain_prim_vars, &
                                                     exit_pressure=desired_boundary_pressure, &
                                                     boundary_norm=boundary_norm)
              else
                boundary_prim_vars = supersonic_outlet(domain_prim_vars=domain_prim_vars)
              end if
            else ! inlet
              if(abs(mach_u) > 1.0_rk) then
                boundary_prim_vars = [desired_boundary_density, domain_u, domain_v, desired_boundary_pressure]
                ! error stop 'Error in pressure_input_bc_t%apply_pressure_input_primitive_var_bc(): ' // &
                !            'Supersonic inlet not configured yet'
              else
                boundary_prim_vars = subsonic_inlet(domain_prim_vars=domain_prim_vars, &
                                                    boundary_norm=boundary_norm, &
                                                    inlet_density=desired_boundary_density, &
                                                    inlet_total_press=desired_boundary_pressure, &
                                                    inlet_flow_angle=0.0_rk)
              end if
            end if
          end associate
          self%edge_rho(j) = boundary_prim_vars(1)
          self%edge_u(j) = boundary_prim_vars(2)
          self%edge_v(j) = 0.0_rk
          self%edge_p(j) = boundary_prim_vars(4)
        end do

        do i = 1, self%n_ghost_layers
          rho(self%ihi_ghost(i), bottom:top) = self%edge_rho(bottom:top)
          u(self%ihi_ghost(i), bottom:top) = self%edge_u(bottom:top)
          v(self%ihi_ghost(i), bottom:top) = self%edge_v(bottom:top)
          p(self%ihi_ghost(i), bottom:top) = self%edge_p(bottom:top)
        end do
      case default
        error stop "Unsupported location to apply the bc at in "// &
          "pressure_input_bc_t%apply_pressure_input_primitive_var_bc()"
      end select

    end associate
    if(allocated(edge_pressure)) deallocate(edge_pressure)

  end subroutine apply_pressure_input_primitive_var_bc

  subroutine finalize(self)
    type(pressure_input_bc_t), intent(inout) :: self
    call debug_print('Running pressure_input_bc_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%ilo_ghost)) deallocate(self%ilo_ghost)
    if(allocated(self%ihi_ghost)) deallocate(self%ihi_ghost)
    if(allocated(self%jlo_ghost)) deallocate(self%jlo_ghost)
    if(allocated(self%jhi_ghost)) deallocate(self%jhi_ghost)

    if(allocated(self%edge_rho)) deallocate(self%edge_rho)
    if(allocated(self%edge_u)) deallocate(self%edge_u)
    if(allocated(self%edge_v)) deallocate(self%edge_v)
    if(allocated(self%edge_p)) deallocate(self%edge_p)
  end subroutine finalize
end module mod_pressure_input_bc
