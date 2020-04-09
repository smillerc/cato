module mod_vacuum_bc
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use mod_floating_point_utils, only: near_zero

  implicit none

  private
  public :: vacuum_bc_t, vacuum_bc_constructor

  type, extends(boundary_condition_t) :: vacuum_bc_t
    real(rk) :: vacuum_density = 1e-5_rk
    real(rk) :: vacuum_pressure = 1e3 ! approx 1 atmosphere

    real(rk), dimension(:, :), allocatable :: edge_primitive_vars
  contains
    procedure, public :: apply_primitive_var_bc => apply_vacuum_primitive_var_bc
    procedure, public :: apply_reconstructed_state_bc => apply_vacuum_reconstructed_state_bc
    procedure, public :: apply_cell_gradient_bc => apply_vacuum_cell_gradient_bc
    procedure, public :: copy => copy_vacuum_bc
  end type
contains

  function vacuum_bc_constructor(location, input) result(bc)
    type(vacuum_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input

    allocate(bc)
    bc%name = 'vacuum'
    bc%location = location
  end function vacuum_bc_constructor

  subroutine copy_vacuum_bc(out_bc, in_bc)
    class(boundary_condition_t), intent(in) :: in_bc
    class(vacuum_bc_t), intent(inout) :: out_bc
  end subroutine

  subroutine apply_vacuum_primitive_var_bc(self, primitive_vars, lbounds)
    !< Apply vacuum boundary conditions to the conserved state vector field
    class(vacuum_bc_t), intent(inout) :: self
    integer(ik), dimension(3), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(inout) :: primitive_vars
    !< ((rho, u ,v, p), i, j); Conserved variables for each cell

    real(rk) :: min_density
    real(rk) :: min_pressure

    integer(ik) :: left         !< Min i real cell index
    integer(ik) :: right        !< Max i real cell index
    integer(ik) :: bottom       !< Min j real cell index
    integer(ik) :: top          !< Max j real cell index
    integer(ik) :: left_ghost   !< Min i ghost cell index
    integer(ik) :: right_ghost  !< Max i ghost cell index
    integer(ik) :: bottom_ghost !< Min j ghost cell index
    integer(ik) :: top_ghost    !< Max j ghost cell index

    integer(ik) :: i
    real(rk) :: gamma, cs, mach
    real(rk), dimension(:), allocatable :: edge_velocity
    real(rk), dimension(:), allocatable :: edge_density

    gamma = eos%get_gamma()

    left_ghost = lbound(primitive_vars, dim=2)
    right_ghost = ubound(primitive_vars, dim=2)
    bottom_ghost = lbound(primitive_vars, dim=3)
    top_ghost = ubound(primitive_vars, dim=3)
    left = left_ghost + 1
    right = right_ghost - 1
    bottom = bottom_ghost + 1
    top = top_ghost - 1

    select case(self%location)
    case('+x')
      call debug_print('Running vacuum_bc_t%apply_vacuum_primitive_var_bc() +x', __FILE__, __LINE__)

      ! Get the density (will be overwitten below)
      allocate(self%edge_primitive_vars(4, bottom_ghost:top_ghost))
      allocate(edge_velocity(bottom:top))
      allocate(edge_density(bottom:top))

      self%edge_primitive_vars(1, bottom:top) = self%vacuum_density ! primitive_vars(1, right, bottom:top)

      ! Zero-gradient in velocity
      self%edge_primitive_vars(2, bottom:top) = primitive_vars(2, right, bottom:top)
      self%edge_primitive_vars(3, bottom:top) = primitive_vars(3, right, bottom:top)

      ! Set the edge pressure to "vacuum"
      self%edge_primitive_vars(4, bottom:top) = self%vacuum_pressure ! primitive_vars(4, right, bottom:top)

      ! associate(u=>primitive_vars(2, right, bottom:top), &
      !           v=>primitive_vars(3, right, bottom:top), &
      !           desired_mach=>2.5_rk, P => self%vacuum_pressure)

      ! Set the edge density (if the density is greater than the imposed "vacuum" value)
      ! do i = bottom, top

      ! edge_velocity(i) = sqrt(u(i)**2 + v(i)**2)
      ! cs = eos%sound_speed(pressure=primitive_vars(4, right, i), density=primitive_vars(1, right, i))
      ! mach = abs(edge_velocity(i)) / cs
      ! ! print*, i, bottom, top, edge_velocity(i), cs, mach

      ! ! Make it supersonic
      ! if (mach < desired_mach .and. .not. near_zero(edge_velocity(i))) then
      !   ! print*, 'Attempting to make the bc supersonic'
      !   ! Set density to make the flow supersonic
      !   if (abs(edge_velocity(i)) > 0.0_rk) then
      !     edge_density(i) = abs((gamma * P) / (edge_velocity(i) / desired_mach))
      !   end if

      !   if (self%edge_primitive_vars(1, i) > edge_density(i)) then
      !     self%edge_primitive_vars(1, i) = edge_density(i)
      !   end if

      !   self%edge_primitive_vars(4, i) = self%vacuum_pressure

      ! ! else ! the edge is already supersonic, so make it zero gradient
      ! !   print*, 'bc is already supersonic, Mach=', mach
      ! end if
      ! end do
      ! end associate

      ! Now apply the final values to the ghost region
      primitive_vars(:, right_ghost, bottom:top) = self%edge_primitive_vars(:, bottom:top)
      ! print*, 'rho', primitive_vars(1, right_ghost, bottom:top)
      ! print*, 'u',   primitive_vars(2, right_ghost, bottom:top)
      ! print*, 'v',   primitive_vars(3, right_ghost, bottom:top)
      ! print*, 'p',   primitive_vars(4, right_ghost, bottom:top)
    case default
      error stop "Unsupported location to apply the bc at in vacuum_bc_t%apply_vacuum_cell_gradient_bc()"
    end select

    if(allocated(edge_velocity)) deallocate(edge_velocity)
    if(allocated(edge_density)) deallocate(edge_density)
  end subroutine apply_vacuum_primitive_var_bc

  subroutine apply_vacuum_reconstructed_state_bc(self, reconstructed_state, lbounds)
    !< Apply vacuum boundary conditions to the reconstructed state vector field

    class(vacuum_bc_t), intent(inout) :: self
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
      call debug_print('Running vacuum_bc_t%apply_vacuum_reconstructed_state_bc() +x', __FILE__, __LINE__)
      do n = 1, 2
        do p = 1, n_points
          reconstructed_state(:, p, n, right_ghost, :) = self%edge_primitive_vars
        end do
      end do

      ! case('-x')
      !   call debug_print('Running vacuum_bc_t%apply_vacuum_reconstructed_state_bc() -x', __FILE__, __LINE__)
      !   reconstructed_state(:, :, :, left_ghost, :) = reconstructed_state(:, :, :, left, :)
      ! case('+y')
      !   call debug_print('Running vacuum_bc_t%apply_vacuum_reconstructed_state_bc() +y', __FILE__, __LINE__)
      !   reconstructed_state(:, :, :, :, top_ghost) = reconstructed_state(:, :, :, :, top)
      ! case('-y')
      !   call debug_print('Running vacuum_bc_t%apply_vacuum_reconstructed_state_bc() -y', __FILE__, __LINE__)
      !   reconstructed_state(:, :, :, :, bottom_ghost) = reconstructed_state(:, :, :, :, bottom)
    case default
      error stop "Unsupported location to apply the bc at in vacuum_bc_t%apply_vacuum_reconstructed_state_bc()"
    end select

    if(allocated(self%edge_primitive_vars)) deallocate(self%edge_primitive_vars)
  end subroutine apply_vacuum_reconstructed_state_bc

  subroutine apply_vacuum_cell_gradient_bc(self, cell_gradient, lbounds)
    !< Apply vacuum boundary conditions to the reconstructed state vector field

    class(vacuum_bc_t), intent(in) :: self
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
      cell_gradient(:, :, right_ghost, :) = 0.0_rk
    case('-x')
      cell_gradient(:, :, left_ghost, :) = 0.0_rk
    case('+y')
      cell_gradient(:, :, :, top_ghost) = 0.0_rk
    case('-y')
      cell_gradient(:, :, :, bottom_ghost) = 0.0_rk
    case default
      error stop "Unsupported location to apply the bc at in vacuum_bc_t%apply_vacuum_cell_gradient_bc()"
    end select

  end subroutine apply_vacuum_cell_gradient_bc
end module mod_vacuum_bc
