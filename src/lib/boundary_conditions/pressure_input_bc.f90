module mod_pressure_input_bc
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_input, only: input_t
  use linear_interpolation_module, only: linear_interp_1d
  use mod_eos, only: eos
  implicit none

  private
  public :: pressure_input_bc_t, pressure_input_bc_constructor

  type, extends(boundary_condition_t) :: pressure_input_bc_t
    logical :: constant_pressure = .false.
    real(rk) :: pressure_input = 0.0_rk
    real(rk) :: scale_factor = 1.0_rk !< scaling factor (e.g. 1x, 2x) to scale the input up/down

    ! Temporal inputs
    type(linear_interp_1d) :: temporal_pressure_input
    character(len=100) :: input_filename

    real(rk), dimension(:, :), allocatable :: edge_primitive_vars
    !< ((rho, u ,v, p), i); Conserved variables for the ghost cells along the boundary. This is saved b/c
    !< it is reapplied to the reconstructed bc state for the ghost cells (besides for the conserved var bc)

  contains
    procedure, private :: read_pressure_input
    procedure, private :: get_desired_pressure
    procedure, public :: apply_primitive_var_bc => apply_pressure_input_primitive_var_bc
    procedure, public :: apply_reconstructed_state_bc => apply_pressure_input_reconstructed_state_bc
    procedure, public :: apply_cell_gradient_bc => apply_pressure_input_cell_gradient_bc
    procedure, public :: copy => copy_pressure_input_bc
  end type
contains

  function pressure_input_bc_constructor(location, input) result(bc)
    type(pressure_input_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input

    allocate(bc)
    bc%name = 'pressure_input'
    bc%location = location
    bc%constant_pressure = input%apply_constant_bc_pressure
    bc%scale_factor = input%bc_pressure_scale_factor
    if(.not. bc%constant_pressure) then
      bc%input_filename = trim(input%bc_pressure_input_file)
      call bc%read_pressure_input()
    else
      bc%pressure_input = input%constant_bc_pressure_value
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
  end subroutine

  subroutine copy_pressure_input_bc(out_bc, in_bc)
    class(boundary_condition_t), intent(in) :: in_bc
    class(pressure_input_bc_t), intent(inout) :: out_bc

    call debug_print('Running pressure_input_bc_t%copy_pressure_input_bc()', __FILE__, __LINE__)

    out_bc%name = in_bc%name
    out_bc%location = in_bc%location
    select type(in_bc)
    class is(pressure_input_bc_t)
      out_bc%constant_pressure = in_bc%constant_pressure
      out_bc%pressure_input = in_bc%pressure_input
      out_bc%temporal_pressure_input = in_bc%temporal_pressure_input
      out_bc%input_filename = in_bc%input_filename
    class default
      error stop 'pressure_input_bc_t%copy_pressure_input_bc: unsupported in_bc class'
    end select

  end subroutine

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

    desired_pressure = desired_pressure * self%scale_factor
  end function get_desired_pressure

  subroutine apply_pressure_input_primitive_var_bc(self, primitive_vars, lbounds)
    !< Apply pressure_input boundary conditions to the conserved state vector field
    class(pressure_input_bc_t), intent(inout) :: self
    integer(ik), dimension(3), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(inout) :: primitive_vars
    !< ((rho, u ,v, p), i, j); Conserved variables for each cell

    integer(ik) :: left         !< Min i real cell index
    integer(ik) :: right        !< Max i real cell index
    integer(ik) :: bottom       !< Min j real cell index
    integer(ik) :: top          !< Max j real cell index
    integer(ik) :: left_ghost   !< Min i ghost cell index
    integer(ik) :: right_ghost  !< Max i ghost cell index
    integer(ik) :: bottom_ghost !< Min j ghost cell index
    integer(ik) :: top_ghost    !< Max j ghost cell index

    logical :: inflow
    logical :: outflow

    real(rk) :: desired_boundary_pressure, boundary_density, gamma, ave_u, min_u, max_u, u, v, rho_new, ave_rho
    real(rk) :: mach
    real(rk), dimension(:), allocatable :: edge_pressure

    gamma = eos%get_gamma()

    left_ghost = lbound(primitive_vars, dim=2)
    right_ghost = ubound(primitive_vars, dim=2)
    bottom_ghost = lbound(primitive_vars, dim=3)
    top_ghost = ubound(primitive_vars, dim=3)
    left = left_ghost + 1
    right = right_ghost - 1
    bottom = bottom_ghost + 1
    top = top_ghost - 1

    inflow = .false.
    outflow = .false.
    mach = 1.5_rk
    ! Since we know the pressure/density in the fluid edge cell,
    ! and we also know what we want the pressure to be on the boundary (ghost)
    ! We use the isentropic relation rho2/rho2 = (P2/P1)^(gamma-1) to find the density
    ! in the ghost zone
    !           |
    !  Fluid    |   Ghost
    !  P1,      |    P2
    !  rho_1    |  rho2
    !           |

    select case(self%location)
    case('+x')
      call debug_print('Running pressure_input_bc_t%apply_pressure_input_primitive_var_bc() +x', __FILE__, __LINE__)

      if(allocated(self%edge_primitive_vars)) deallocate(self%edge_primitive_vars)
      allocate(self%edge_primitive_vars(4, bottom_ghost:top_ghost))
      self%edge_primitive_vars = 0.0_rk

      allocate(edge_pressure(bottom:top))

      min_u = minval(primitive_vars(2, right, bottom:top))
      max_u = maxval(primitive_vars(2, right, bottom:top))
      ave_u = sum(primitive_vars(2, right, bottom:top)) / size(primitive_vars(2, right, bottom:top))
      ave_rho = sum(primitive_vars(1, right, bottom:top)) / size(primitive_vars(1, right, bottom:top))
      primitive_vars(1, right, bottom:top) = ave_rho

      ! defaults
      u = ave_u
      v = 0.0_rk

      if(ave_u < 0.0_rk) then
        inflow = .true.
        outflow = .false.
        u = min_u
      else
        inflow = .false.
        outflow = .true.
        u = max_u
      end if

      ! Default to zero-gradient
      self%edge_primitive_vars(1, :) = sum(primitive_vars(1, right, bottom:top)) / size(primitive_vars(1, right, bottom:top))
      self%edge_primitive_vars(2, :) = u
      self%edge_primitive_vars(3, :) = v

      edge_pressure(bottom:top) = sum(primitive_vars(4, right, bottom:top)) / size(primitive_vars(4, right, bottom:top))
      desired_boundary_pressure = self%get_desired_pressure()

      self%edge_primitive_vars(4, bottom:top) = edge_pressure(bottom:top)

      if(desired_boundary_pressure <= 0.0_rk) then
        ! Default to zero-gradient if the input pressure goes <= 0
        self%edge_primitive_vars(1, bottom:top) = primitive_vars(1, right, bottom:top)
        self%edge_primitive_vars(2, bottom:top) = primitive_vars(2, right, bottom:top)
        self%edge_primitive_vars(3, bottom:top) = primitive_vars(3, right, bottom:top)
        self%edge_primitive_vars(4, bottom:top) = primitive_vars(4, right, bottom:top)
      else
        ! Pressure on the right-most column of fluid cells

        ! Find the desired density from isentropic relations rho2 = f(rho1, P1, P2)
        associate(rho_1=>primitive_vars(1, right, bottom:top), &
                  P_1=>edge_pressure(bottom:top), &
                  P_2=>desired_boundary_pressure)

          if(self%get_time() > self%max_time) then  ! switch to a zero gradient
            self%edge_primitive_vars(1, bottom:top) = primitive_vars(1, right, bottom:top)
            self%edge_primitive_vars(4, bottom:top) = edge_pressure
          else

            self%edge_primitive_vars(4, bottom:top) = P_2

            ! Outflow
            if(outflow) then
              ! print*, 'outflow!'
              rho_new = abs((gamma * P_2) / (u / mach)) ! Set density to make the flow supersonic
              if(rho_new > ave_rho) then
                self%edge_primitive_vars(1, bottom:top) = ave_rho
              else
                self%edge_primitive_vars(1, bottom:top) = rho_new
              end if

            else ! Inflow
              ! print*, 'inflow!'
              self%edge_primitive_vars(1, bottom:top) = rho_1 * (P_2 / P_1)**(1.0_rk / gamma)
            end if
            ! write(*,'(a, 12(es10.3,1x))') 'edge p  ', P_1
            ! write(*,'(a, 12(es10.3,1x))') 'input p ', P_2
            ! write(*,'(a, 12(es10.3,1x))') 'edge rho', rho_1
            ! write(*,'(a, 12(es10.3,1x))') 'new rho ', self%edge_primitive_vars(1, bottom:top)
            ! write(*,'(a, 12(es10.3,1x))') 'new u   ', self%edge_primitive_vars(2, bottom:top)
            ! write(*,'(a, 12(es10.3,1x))') 'new p   ', self%edge_primitive_vars(4, bottom:top)
            ! print*

          end if
        end associate

      end if

      primitive_vars(:, right_ghost, bottom:top) = self%edge_primitive_vars(:, bottom:top)
    case default
      error stop "Unsupported location to apply the bc at in "// &
        "pressure_input_bc_t%apply_pressure_input_primitive_var_bc()"
    end select

    if(allocated(edge_pressure)) deallocate(edge_pressure)

  end subroutine apply_pressure_input_primitive_var_bc

  subroutine apply_pressure_input_reconstructed_state_bc(self, reconstructed_state, lbounds)
    !< Apply pressure_input boundary conditions to the reconstructed state vector field

    class(pressure_input_bc_t), intent(inout) :: self
    integer(ik), dimension(5), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):, &
                        lbounds(4):, lbounds(5):), intent(inout) :: reconstructed_state
    !< ((rho, u ,v, p), point, node/midpoint, i, j); Reconstructed state for each cell

    integer(ik) :: left         !< Min i real cell index
    integer(ik) :: right        !< Max i real cell index
    integer(ik) :: bottom       !< Min j real cell index
    integer(ik) :: top          !< Max j real cell index
    integer(ik) :: left_ghost   !< Min i ghost cell index
    integer(ik) :: right_ghost  !< Max i ghost cell index
    integer(ik) :: bottom_ghost !< Min j ghost cell index
    integer(ik) :: top_ghost    !< Max j ghost cell index
    integer(ik) :: n, p, n_points

    n_points = ubound(reconstructed_state, dim=2)
    left_ghost = lbound(reconstructed_state, dim=4)
    right_ghost = ubound(reconstructed_state, dim=4)
    bottom_ghost = lbound(reconstructed_state, dim=5)
    top_ghost = ubound(reconstructed_state, dim=5)
    left = left_ghost + 1
    right = right_ghost - 1
    bottom = bottom_ghost + 1
    top = top_ghost - 1

    select case(self%location)
    case('+x')
      call debug_print('Running pressure_input_bc_t%apply_pressure_input_reconstructed_state_bc() +x', __FILE__, __LINE__)
      do n = 1, 2
        do p = 1, n_points
          reconstructed_state(:, p, n, right_ghost, :) = self%edge_primitive_vars
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

    if(allocated(self%edge_primitive_vars)) deallocate(self%edge_primitive_vars)
  end subroutine apply_pressure_input_reconstructed_state_bc

  subroutine apply_pressure_input_cell_gradient_bc(self, cell_gradient, lbounds)
    !< Apply pressure_input boundary conditions to the reconstructed state vector field

    class(pressure_input_bc_t), intent(in) :: self
    integer(ik), dimension(4), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):, &
                        lbounds(4):), intent(inout) :: cell_gradient
    !< ((rho, u ,v, p), (d/dx, d/dy), i, j); Gradient of each cell's primitive variables

    integer(ik) :: left         !< Min i real cell index
    integer(ik) :: right        !< Max i real cell index
    integer(ik) :: bottom       !< Min j real cell index
    integer(ik) :: top          !< Max j real cell index
    integer(ik) :: left_ghost   !< Min i ghost cell index
    integer(ik) :: right_ghost  !< Max i ghost cell index
    integer(ik) :: bottom_ghost !< Min j ghost cell index
    integer(ik) :: top_ghost    !< Max j ghost cell index

    left_ghost = lbound(cell_gradient, dim=3)
    right_ghost = ubound(cell_gradient, dim=3)
    bottom_ghost = lbound(cell_gradient, dim=4)
    top_ghost = ubound(cell_gradient, dim=4)
    left = left_ghost + 1
    right = right_ghost - 1
    bottom = bottom_ghost + 1
    top = top_ghost - 1

    select case(self%location)
    case('+x')
      call debug_print('Running pressure_input_bc_t%apply_periodic_cell_gradient_bc() +x', __FILE__, __LINE__)
      cell_gradient(:, :, right_ghost, :) = 0.0_rk
      ! cell_gradient(1, :, right_ghost, :) = -cell_gradient(1, :, right, :)
      ! cell_gradient(2:3, :, right_ghost, :) = 0.0_rk !cell_gradient(2:3, :, right, :)
      ! cell_gradient(4, :, right_ghost, :) = -cell_gradient(4, :, right, :)
    case('-x')
      call debug_print('Running pressure_input_bc_t%apply_periodic_cell_gradient_bc() -x', __FILE__, __LINE__)
      cell_gradient(:, :, right_ghost, :) = 0.0_rk
    case('+y')
      call debug_print('Running pressure_input_bc_t%apply_periodic_cell_gradient_bc() +y', __FILE__, __LINE__)
      cell_gradient(:, :, :, top_ghost) = 0.0_rk
    case('-y')
      call debug_print('Running pressure_input_bc_t%apply_periodic_cell_gradient_bc() -y', __FILE__, __LINE__)
      cell_gradient(:, :, :, bottom_ghost) = 0.0_rk
    case default
      error stop "Unsupported location to apply the bc at in pressure_input_bc_t%apply_pressure_input_cell_gradient_bc()"
    end select

  end subroutine apply_pressure_input_cell_gradient_bc
end module mod_pressure_input_bc
