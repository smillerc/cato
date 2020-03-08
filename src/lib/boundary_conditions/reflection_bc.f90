module mod_reflection_bc
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_input, only: input_t

  implicit none

  private
  public :: reflection_bc_t, reflection_bc_constructor

  type, extends(boundary_condition_t) :: reflection_bc_t
  contains
    procedure, public :: apply_primitive_var_bc => apply_reflection_primitive_var_bc
    procedure, public :: apply_reconstructed_state_bc => apply_reflection_reconstructed_state_bc
    procedure, public :: apply_cell_gradient_bc => apply_reflection_cell_gradient_bc
    procedure, public :: copy => copy_reflection_bc
  end type
contains
  function reflection_bc_constructor(location, input) result(bc)
    type(reflection_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input

    allocate(bc)
    bc%name = 'reflection'
    bc%location = location
  end function reflection_bc_constructor

  subroutine copy_reflection_bc(out_bc, in_bc)
    class(boundary_condition_t), intent(in) :: in_bc
    class(reflection_bc_t), intent(inout) :: out_bc
  end subroutine

  subroutine apply_reflection_primitive_var_bc(self, primitive_vars, lbounds)
    !< Apply reflection boundary conditions to the conserved state vector field
    class(reflection_bc_t), intent(inout) :: self
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

    left_ghost = lbound(primitive_vars, dim=2)
    right_ghost = ubound(primitive_vars, dim=2)
    bottom_ghost = lbound(primitive_vars, dim=3)
    top_ghost = ubound(primitive_vars, dim=3)
    left = left_ghost + 1
    right = right_ghost - 1
    bottom = bottom_ghost + 1
    top = top_ghost - 1

    if(size(primitive_vars, dim=1) /= 4) then
      error stop "Error in reflection_bc_t%apply_reflection_primitive_var_bc(), dimension 1 /= 4 (rho,u,v,p)"
    end if

    select case(self%location)
    case('+x')
      call debug_print('Running reflection_bc_t%apply_reflection_primitive_var_bc() +x', __FILE__, __LINE__)
      primitive_vars(1, right_ghost, :) = primitive_vars(1, right, :)      ! density
      primitive_vars(2, right_ghost, :) = -primitive_vars(2, right, :)     ! x velocity
      primitive_vars(3, right_ghost, :) = primitive_vars(3, right, :)      ! y velocity
      primitive_vars(4, right_ghost, :) = primitive_vars(4, right, :)      ! pressure
    case('-x')
      call debug_print('Running reflection_bc_t%apply_reflection_primitive_var_bc() -x', __FILE__, __LINE__)
      primitive_vars(1, left_ghost, :) = primitive_vars(1, left, :)      ! density
      primitive_vars(2, left_ghost, :) = -primitive_vars(2, left, :)     ! x velocity
      primitive_vars(3, left_ghost, :) = primitive_vars(3, left, :)      ! y velocity
      primitive_vars(4, left_ghost, :) = primitive_vars(4, left, :)      ! pressure
    case('+y')
      call debug_print('Running reflection_bc_t%apply_reflection_primitive_var_bc() +y', __FILE__, __LINE__)
      primitive_vars(1, :, top_ghost) = primitive_vars(1, :, top)      ! density
      primitive_vars(2, :, top_ghost) = primitive_vars(2, :, top)      ! x velocity
      primitive_vars(3, :, top_ghost) = -primitive_vars(3, :, top)     ! y velocity
      primitive_vars(4, :, top_ghost) = primitive_vars(4, :, top)      ! pressure
    case('-y')
      call debug_print('Running reflection_bc_t%apply_reflection_primitive_var_bc() -y', __FILE__, __LINE__)
      primitive_vars(1, :, bottom_ghost) = primitive_vars(1, :, bottom)      ! density
      primitive_vars(2, :, bottom_ghost) = primitive_vars(2, :, bottom)      ! x velocity
      primitive_vars(3, :, bottom_ghost) = -primitive_vars(3, :, bottom)     ! y velocity
      primitive_vars(4, :, bottom_ghost) = primitive_vars(4, :, bottom)      ! pressure
    case default
      error stop "Unsupported location to apply the bc at in symmetry_bc_t%apply_symmetry_cell_gradient_bc()"
    end select

  end subroutine apply_reflection_primitive_var_bc

  subroutine apply_reflection_reconstructed_state_bc(self, reconstructed_state, lbounds)
    !< Apply reflection boundary conditions to the reconstructed state vector field

    class(reflection_bc_t), intent(inout) :: self
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

    left_ghost = lbound(reconstructed_state, dim=4)
    right_ghost = ubound(reconstructed_state, dim=4)
    bottom_ghost = lbound(reconstructed_state, dim=5)
    top_ghost = ubound(reconstructed_state, dim=5)
    left = left_ghost + 1
    right = right_ghost - 1
    bottom = bottom_ghost + 1
    top = top_ghost - 1

    if(size(reconstructed_state, dim=1) /= 4) then
      error stop "Error in reflection_bc_t%apply_reflection_cell_gradient_bc(), dimension 1 /= 4 (rho,u,v,p)"
    end if

    if(size(reconstructed_state, dim=2) /= 4) then
      error stop "Error in reflection_bc_t%apply_reflection_cell_gradient_bc(), dimension 2 /= 2 (node 1 - 4)"
    end if

    if(size(reconstructed_state, dim=3) /= 2) then
      error stop "Error in reflection_bc_t%apply_reflection_cell_gradient_bc(), dimension 3 /= 2 (corner, midpoint)"
    end if

    select case(self%location)
    case('+x')
      call debug_print('Running reflection_bc_t%apply_reflection_reconstructed_state_bc() +x', __FILE__, __LINE__)
      reconstructed_state(1, :, :, right_ghost, :) = reconstructed_state(1, :, :, right, :)  ! density
      reconstructed_state(2, :, :, right_ghost, :) = -reconstructed_state(2, :, :, right, :) ! x velocity
      reconstructed_state(3, :, :, right_ghost, :) = reconstructed_state(3, :, :, right, :)  ! y velocity
      reconstructed_state(4, :, :, right_ghost, :) = reconstructed_state(4, :, :, right, :)  ! pressure
    case('-x')
      call debug_print('Running reflection_bc_t%apply_reflection_reconstructed_state_bc() -x', __FILE__, __LINE__)
      reconstructed_state(1, :, :, left_ghost, :) = reconstructed_state(1, :, :, left, :)  ! density
      reconstructed_state(2, :, :, left_ghost, :) = -reconstructed_state(2, :, :, left, :) ! x velocity
      reconstructed_state(3, :, :, left_ghost, :) = reconstructed_state(3, :, :, left, :)  ! y velocity
      reconstructed_state(4, :, :, left_ghost, :) = reconstructed_state(4, :, :, left, :)  ! pressure
    case('+y')
      call debug_print('Running reflection_bc_t%apply_reflection_reconstructed_state_bc() +y', __FILE__, __LINE__)
      reconstructed_state(1, :, :, :, top_ghost) = reconstructed_state(1, :, :, :, top)  ! density
      reconstructed_state(2, :, :, :, top_ghost) = reconstructed_state(2, :, :, :, top)  ! x velocity
      reconstructed_state(3, :, :, :, top_ghost) = -reconstructed_state(3, :, :, :, top) ! y velocity
      reconstructed_state(4, :, :, :, top_ghost) = reconstructed_state(4, :, :, :, top)  ! pressure
    case('-y')
      call debug_print('Running reflection_bc_t%apply_reflection_reconstructed_state_bc() -y', __FILE__, __LINE__)
      reconstructed_state(1, :, :, :, bottom_ghost) = reconstructed_state(1, :, :, :, bottom)  ! density
      reconstructed_state(2, :, :, :, bottom_ghost) = reconstructed_state(2, :, :, :, bottom)  ! x velocity
      reconstructed_state(3, :, :, :, bottom_ghost) = -reconstructed_state(3, :, :, :, bottom) ! y velocity
      reconstructed_state(4, :, :, :, bottom_ghost) = reconstructed_state(4, :, :, :, bottom)  ! pressure
    case default
      error stop "Unsupported location to apply the bc at in reflection_bc_t%apply_symmetry_reconstructed_state_bc()"
    end select

  end subroutine apply_reflection_reconstructed_state_bc

  subroutine apply_reflection_cell_gradient_bc(self, cell_gradient, lbounds)
    !< Apply reflection boundary conditions to the reconstructed state vector field

    class(reflection_bc_t), intent(in) :: self
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

    if(size(cell_gradient, dim=1) /= 4) then
      error stop "Error in reflection_bc_t%apply_reflection_cell_gradient_bc(), dimension 1 /= 4 (rho,u,v,p)"
    end if

    if(size(cell_gradient, dim=2) /= 2) then
      error stop "Error in reflection_bc_t%apply_reflection_cell_gradient_bc(), dimension 2 /= 2 (d/dx, d/dy)"
    end if

    select case(self%location)
    case('+x')
      call debug_print('Running reflection_bc_t%apply_periodic_cell_gradient_bc() +x', __FILE__, __LINE__)
      cell_gradient(:, 1, right_ghost, :) = -cell_gradient(:, 1, right, :)
      cell_gradient(:, 2, right_ghost, :) = cell_gradient(:, 2, right, :)
    case('-x')
      call debug_print('Running reflection_bc_t%apply_periodic_cell_gradient_bc() -x', __FILE__, __LINE__)
      cell_gradient(:, 1, left_ghost, :) = -cell_gradient(:, 1, left, :)
      cell_gradient(:, 2, left_ghost, :) = cell_gradient(:, 2, left, :)
    case('+y')
      call debug_print('Running reflection_bc_t%apply_periodic_cell_gradient_bc() +y', __FILE__, __LINE__)
      cell_gradient(:, 1, :, top_ghost) = cell_gradient(:, 1, :, top)
      cell_gradient(:, 2, :, top_ghost) = -cell_gradient(:, 2, :, top)
    case('-y')
      call debug_print('Running reflection_bc_t%apply_periodic_cell_gradient_bc() -y', __FILE__, __LINE__)
      cell_gradient(:, 1, :, bottom_ghost) = cell_gradient(:, 1, :, bottom)
      cell_gradient(:, 2, :, bottom_ghost) = -cell_gradient(:, 2, :, bottom)
    case default
      error stop "Unsupported location to apply the bc at in reflection_bc_t%apply_reflection_cell_gradient_bc()"
    end select

  end subroutine apply_reflection_cell_gradient_bc

end module mod_reflection_bc
