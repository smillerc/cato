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
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
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
    real(rk) :: pressure_input = 0.0_rk
    real(rk) :: temperature_input = 273.0_rk ! K
    real(rk) :: density_input = 5e-3_rk !< Boundary density g/cc
    real(rk) :: scale_factor = 1.0_rk !< scaling factor (e.g. 1x, 2x) to scale the input up/down

    ! Temporal inputs
    type(linear_interp_1d) :: temporal_pressure_input
    character(len=100) :: input_filename

    real(rk), dimension(:), allocatable :: edge_rho !< boundary (ghost/edge/etc) density
    real(rk), dimension(:), allocatable :: edge_u   !< boundary (ghost/edge/etc) x-velocity
    real(rk), dimension(:), allocatable :: edge_v   !< boundary (ghost/edge/etc) y-velocity
    real(rk), dimension(:), allocatable :: edge_p   !< boundary (ghost/edge/etc) pressure

  contains
    procedure, private :: read_pressure_input
    procedure, private :: get_desired_pressure
    procedure, public :: apply_primitive_var_bc => apply_pressure_input_primitive_var_bc
    procedure, public :: apply_reconstructed_state_bc => apply_pressure_input_reconstructed_state_bc
    procedure, public :: apply_gradient_bc
    final :: finalize
  end type
contains

  function pressure_input_bc_constructor(location, input, grid) result(bc)
    type(pressure_input_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input
    class(grid_t), intent(in) :: grid

    allocate(bc)
    bc%name = 'pressure_input'
    bc%location = location
    call bc%set_indices(grid)

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

    real(rk), dimension(:), allocatable :: time_sec
    real(rk), dimension(:), allocatable :: pressure_barye
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
    time_sec = 0.0_rk
    pressure_barye = 0.0_rk

    do line = 1, nlines
      read(input_unit, *, iostat=io_status) time_sec(line), pressure_barye(line)
      if(io_status /= 0) exit
    end do
    close(input_unit)

    ! Apply scaling and non-dimensionalization
    time_sec = time_sec / t_0
    pressure_barye = (pressure_barye / p_0) * self%scale_factor

    self%max_time = maxval(time_sec)

    ! Initialize the linear interpolated data object so we can query the pressure at any time
    call self%temporal_pressure_input%initialize(time_sec, pressure_barye, interp_status)
    if(interp_status /= 0) then
      write(*, *) time_sec
      write(*, *) pressure_barye
      write(*, *) interp_status
      error stop "Error initializing pressure_input_bc_t%temporal_pressure_input"
    end if

    deallocate(time_sec)
    deallocate(pressure_barye)
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
  end function get_desired_pressure

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

    real(rk) :: desired_boundary_pressure, boundary_density, gamma
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

    if(allocated(self%edge_rho)) deallocate(self%edge_rho)
    if(allocated(self%edge_u)) deallocate(self%edge_u)
    if(allocated(self%edge_v)) deallocate(self%edge_v)
    if(allocated(self%edge_p)) deallocate(self%edge_p)

    associate(left=>self%ilo, right=>self%ihi, bottom=>self%jlo, top=>self%jhi, &
              left_ghost=>self%ilo_ghost, right_ghost=>self%ihi_ghost, &
              bottom_ghost=>self%jlo_ghost, top_ghost=>self%jhi_ghost, &
              n=>self%n_ghost_layers)

      select case(self%location)
      case('+x', '-x')
        allocate(self%edge_rho(bottom_ghost(1):top_ghost(n)))
        allocate(self%edge_u(bottom_ghost(1):top_ghost(n)))
        allocate(self%edge_v(bottom_ghost(1):top_ghost(n)))
        allocate(self%edge_p(bottom_ghost(1):top_ghost(n)))
        allocate(domain_rho(bottom_ghost(1):top_ghost(n)))
        allocate(domain_u(bottom_ghost(1):top_ghost(n)))
        allocate(domain_v(bottom_ghost(1):top_ghost(n)))
        allocate(domain_p(bottom_ghost(1):top_ghost(n)))
      case('+y', '-y')
        allocate(self%edge_rho(left_ghost(1):right_ghost(n)))
        allocate(self%edge_u(left_ghost(1):right_ghost(n)))
        allocate(self%edge_v(left_ghost(1):right_ghost(n)))
        allocate(self%edge_p(left_ghost(1):right_ghost(n)))
        allocate(domain_rho(left_ghost(1):right_ghost(n)))
        allocate(domain_u(left_ghost(1):right_ghost(n)))
        allocate(domain_v(left_ghost(1):right_ghost(n)))
        allocate(domain_p(left_ghost(1):right_ghost(n)))
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

        domain_rho = rho(right, bottom_ghost(1):top_ghost(n))
        domain_u = u(right, bottom_ghost(1):top_ghost(n))
        domain_v = 0.0_rk
        domain_p = p(right, bottom_ghost(1):top_ghost(n))
        boundary_norm = [1.0_rk, 0.0_rk] ! outward normal vector at +x

        do j = bottom, top
          associate(rho_d=>domain_rho(j), u_d=>domain_u(j), &
                    v_d=>domain_v(j), p_d=>domain_p(j))
            domain_prim_vars = [rho_d, u_d, v_d, p_d]

            call eos%sound_speed(p=p_d, rho=rho_d, cs=cs)
            mach_u = u_d / cs

            if(mach_u > 0.0_rk) then ! outlet
              if(mach_u > 1.0_rk) then

                boundary_prim_vars = supersonic_outlet(domain_prim_vars=domain_prim_vars)
              else
                boundary_prim_vars = subsonic_outlet(domain_prim_vars=domain_prim_vars, &
                                                     exit_pressure=desired_boundary_pressure, &
                                                     boundary_norm=boundary_norm)
              end if
            else ! inlet
              if(abs(mach_u) > 1.0_rk) then
              error stop 'Error in pressure_input_bc_t%apply_pressure_input_primitive_var_bc(): Supersonic inlet not configured yet'
              else
                boundary_prim_vars = subsonic_inlet(domain_prim_vars=domain_prim_vars, &
                                                    boundary_norm=boundary_norm, &
                                                    inlet_density=self%density_input, &
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
          rho(right_ghost(i), bottom:top) = self%edge_rho(bottom:top)
          u(right_ghost(i), bottom:top) = self%edge_u(bottom:top)
          v(right_ghost(i), bottom:top) = self%edge_v(bottom:top)
          p(right_ghost(i), bottom:top) = self%edge_p(bottom:top)
        end do
      case default
        error stop "Unsupported location to apply the bc at in "// &
          "pressure_input_bc_t%apply_pressure_input_primitive_var_bc()"
      end select

    end associate
    if(allocated(edge_pressure)) deallocate(edge_pressure)

  end subroutine apply_pressure_input_primitive_var_bc

  subroutine apply_pressure_input_reconstructed_state_bc(self, recon_rho, recon_u, recon_v, recon_p, lbounds)
    !< Apply pressure boundary conditions to the reconstructed state vector field

    class(pressure_input_bc_t), intent(inout) :: self

    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_rho
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_u
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_v
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_p

    integer(ik) :: i, p

    associate(left=>self%ilo, right=>self%ihi, bottom=>self%jlo, top=>self%jhi, &
              left_ghost=>self%ilo_ghost, right_ghost=>self%ihi_ghost, &
              bottom_ghost=>self%jlo_ghost, top_ghost=>self%jhi_ghost)
      select case(self%location)
      case('+x')
        call debug_print('Running pressure_input_bc_t%apply_pressure_input_reconstructed_state_bc() +x', __FILE__, __LINE__)

        do i = 1, self%n_ghost_layers
          do p = lbound(recon_rho, dim=1), ubound(recon_rho, dim=1)
            recon_rho(p, right_ghost(i), :) = self%edge_rho
            recon_u(p, right_ghost(i), :) = self%edge_u
            recon_v(p, right_ghost(i), :) = self%edge_v
            recon_p(p, right_ghost(i), :) = self%edge_p
          end do
        end do

        ! case('-x')
        !   reconstructed_state(:, :, :, left_ghost, top_ghost) = reconstructed_state(:, :, :, right, bottom)
        !   reconstructed_state(:, :, :, left_ghost, bottom_ghost) = reconstructed_state(:, :, :, right, top)
        !   reconstructed_state(:, :, :, left_ghost, bottom:top) = reconstructed_state(:, :, :, right, bottom:top)
        ! case('+y')
        !   reconstructed_state(:, :, :, left_ghost, top_ghost) = reconstructed_state(:, :, :, right, bottom)
        !   reconstructed_state(:, :, :, right_ghost, top_ghost) = reconstructed_state(:, :, :, left, bottom)
        !   reconstructed_state(:, :, :, left:right, top_ghost) = reconstructed_state(:, :, :, left:right, bottom)
        ! case('-y')
        !   reconstructed_state(:, :, :, left_ghost, bottom_ghost) = reconstructed_state(:, :, :, right, top)
        !   reconstructed_state(:, :, :, right_ghost, bottom_ghost) = reconstructed_state(:, :, :, left, top)
        !   reconstructed_state(:, :, :, left:right, bottom_ghost) = reconstructed_state(:, :, :, left:right, top)
      case default
        error stop "Unsupported location to apply the bc at in pressure_input_bc_t%apply_pressure_input_reconstructed_state_bc()"
      end select
    end associate

    if(allocated(self%edge_rho)) deallocate(self%edge_rho)
    if(allocated(self%edge_u)) deallocate(self%edge_u)
    if(allocated(self%edge_v)) deallocate(self%edge_v)
    if(allocated(self%edge_p)) deallocate(self%edge_p)
  end subroutine apply_pressure_input_reconstructed_state_bc

  subroutine apply_gradient_bc(self, grad_x, grad_y, lbounds)
    class(pressure_input_bc_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: grad_x
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: grad_y
    integer(ik) :: i

    associate(left_ghost=>self%ilo_ghost, right_ghost=>self%ihi_ghost, &
              bottom_ghost=>self%jlo_ghost, top_ghost=>self%jhi_ghost)

      select case(self%location)
      case('+x')
        call debug_print('Running pressure_input_bc_t%apply_gradient_bc() +x', __FILE__, __LINE__)
        do i = 1, self%n_ghost_layers
          grad_x(right_ghost(i), :) = 0.0_rk
          grad_y(right_ghost(i), :) = 0.0_rk
        end do

      case('-x')
        call debug_print('Running pressure_input_bc_t%apply_gradient_bc() -x', __FILE__, __LINE__)
        do i = 1, self%n_ghost_layers
          grad_x(left_ghost(i), :) = 0.0_rk
          grad_y(left_ghost(i), :) = 0.0_rk
        end do

      case('+y')
        call debug_print('Running pressure_input_bc_t%apply_gradient_bc() +y', __FILE__, __LINE__)
        do i = 1, self%n_ghost_layers
          grad_x(:, top_ghost(i)) = 0.0_rk
          grad_y(:, top_ghost(i)) = 0.0_rk
        end do
      case('-y')
        call debug_print('Running pressure_input_bc_t%apply_gradient_bc() -y', __FILE__, __LINE__)
        do i = 1, self%n_ghost_layers
          grad_x(:, bottom_ghost(i)) = 0.0_rk
          grad_y(:, bottom_ghost(i)) = 0.0_rk
        end do
      case default
        error stop "Unsupported location to apply the bc at in pressure_input_bc_t%apply_gradient_bc()"
      end select
    end associate

  end subroutine apply_gradient_bc

  subroutine finalize(self)
    type(pressure_input_bc_t), intent(inout) :: self
    call debug_print('Running pressure_input_bc_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%ilo_ghost)) deallocate(self%ilo_ghost)
    if(allocated(self%ihi_ghost)) deallocate(self%ihi_ghost)
    if(allocated(self%jlo_ghost)) deallocate(self%jlo_ghost)
    if(allocated(self%jhi_ghost)) deallocate(self%jhi_ghost)
  end subroutine finalize
end module mod_pressure_input_bc
