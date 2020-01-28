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

    ! integer(ik) :: left         !< Min i real cell index
    ! integer(ik) :: right        !< Max i real cell index
    ! integer(ik) :: bottom       !< Min j real cell index
    ! integer(ik) :: top          !< Max j real cell index
    ! integer(ik) :: left_ghost   !< Min i ghost cell index
    ! integer(ik) :: right_ghost  !< Max i ghost cell index
    ! integer(ik) :: bottom_ghost !< Min j ghost cell index
    ! integer(ik) :: top_ghost    !< Max j ghost cell index

    ! left_ghost = lbound(primitive_vars, dim=2)
    ! right_ghost = ubound(primitive_vars, dim=2)
    ! bottom_ghost = lbound(primitive_vars, dim=3)
    ! top_ghost = ubound(primitive_vars, dim=3)
    ! left = left_ghost + 1
    ! right = right_ghost - 1
    ! bottom = bottom_ghost + 1
    ! top = top_ghost - 1

    ! select case(self%location)
    ! case('+x')
    !   primitive_vars(:, right_ghost, top_ghost) = primitive_vars(:, left, bottom)
    !   primitive_vars(:, right_ghost, bottom_ghost) = primitive_vars(:, left, top)
    !   primitive_vars(:, right_ghost, bottom:top) = primitive_vars(:, left, bottom:top)
    ! case('-x')
    !   primitive_vars(:, left_ghost, top_ghost) = primitive_vars(:, right, bottom)
    !   primitive_vars(:, left_ghost, bottom_ghost) = primitive_vars(:, right, top)
    !   primitive_vars(:, left_ghost, bottom:top) = primitive_vars(:, right, bottom:top)
    ! case('+y')
    !   primitive_vars(:, left_ghost, top_ghost) = primitive_vars(:, right, bottom)
    !   primitive_vars(:, right_ghost, top_ghost) = primitive_vars(:, left, bottom)
    !   primitive_vars(:, left:right, top_ghost) = primitive_vars(:, left:right, bottom)
    ! case('-y')
    !   primitive_vars(:, left_ghost, bottom_ghost) = primitive_vars(:, right, top)
    !   primitive_vars(:, right_ghost, bottom_ghost) = primitive_vars(:, left, top)
    !   primitive_vars(:, left:right, bottom_ghost) = primitive_vars(:, left:right, top)
    ! case default
    !   error stop "Unsupported location to apply the bc at in reflection_bc_t%apply_reflection_cell_gradient_bc()"
    ! end select

  end subroutine apply_reflection_primitive_var_bc

  subroutine apply_reflection_reconstructed_state_bc(self, reconstructed_state, lbounds)
    !< Apply reflection boundary conditions to the reconstructed state vector field

    class(reflection_bc_t), intent(inout) :: self
    integer(ik), dimension(5), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):, &
                        lbounds(4):, lbounds(5):), intent(inout) :: reconstructed_state
    !< ((rho, u ,v, p), point, node/midpoint, i, j); Reconstructed state for each cell

    ! integer(ik) :: left         !< Min i real cell index
    ! integer(ik) :: right        !< Max i real cell index
    ! integer(ik) :: bottom       !< Min j real cell index
    ! integer(ik) :: top          !< Max j real cell index
    ! integer(ik) :: left_ghost   !< Min i ghost cell index
    ! integer(ik) :: right_ghost  !< Max i ghost cell index
    ! integer(ik) :: bottom_ghost !< Min j ghost cell index
    ! integer(ik) :: top_ghost    !< Max j ghost cell index

    ! left_ghost = lbound(reconstructed_state, dim=4)
    ! right_ghost = ubound(reconstructed_state, dim=4)
    ! bottom_ghost = lbound(reconstructed_state, dim=5)
    ! top_ghost = ubound(reconstructed_state, dim=5)
    ! left = left_ghost + 1
    ! right = right_ghost - 1
    ! bottom = bottom_ghost + 1
    ! top = top_ghost - 1

    ! select case(self%location)
    ! case('+x')
    !   reconstructed_state(:, :, :, right_ghost, top_ghost) = reconstructed_state(:, :, :, left, bottom)
    !   reconstructed_state(:, :, :, right_ghost, bottom_ghost) = reconstructed_state(:, :, :, left, top)
    !   reconstructed_state(:, :, :, right_ghost, bottom:top) = reconstructed_state(:, :, :, left, bottom:top)
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
    ! case default
    !   error stop "Unsupported location to apply the bc at in reflection_bc_t%apply_reflection_reconstructed_state_bc()"
    ! end select

  end subroutine apply_reflection_reconstructed_state_bc

  subroutine apply_reflection_cell_gradient_bc(self, cell_gradient, lbounds)
    !< Apply reflection boundary conditions to the reconstructed state vector field

    class(reflection_bc_t), intent(in) :: self
    integer(ik), dimension(4), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):, &
                        lbounds(4):), intent(inout) :: cell_gradient
    !< ((rho, u ,v, p), point, node/midpoint, i, j); Reconstructed state for each cell

    ! integer(ik) :: left         !< Min i real cell index
    ! integer(ik) :: right        !< Max i real cell index
    ! integer(ik) :: bottom       !< Min j real cell index
    ! integer(ik) :: top          !< Max j real cell index
    ! integer(ik) :: left_ghost   !< Min i ghost cell index
    ! integer(ik) :: right_ghost  !< Max i ghost cell index
    ! integer(ik) :: bottom_ghost !< Min j ghost cell index
    ! integer(ik) :: top_ghost    !< Max j ghost cell index

    ! left_ghost = lbound(cell_gradient, dim=3)
    ! right_ghost = ubound(cell_gradient, dim=3)
    ! bottom_ghost = lbound(cell_gradient, dim=4)
    ! top_ghost = ubound(cell_gradient, dim=4)
    ! left = left_ghost + 1
    ! right = right_ghost - 1
    ! bottom = bottom_ghost + 1
    ! top = top_ghost - 1

    ! select case(self%location)
    ! case('+x')
    !   cell_gradient(:, :, right_ghost, top_ghost) = cell_gradient(:, :, left, bottom)
    !   cell_gradient(:, :, right_ghost, bottom_ghost) = cell_gradient(:, :, left, top)
    !   cell_gradient(:, :, right_ghost, bottom:top) = cell_gradient(:, :, left, bottom:top)
    ! case('-x')
    !   cell_gradient(:, :, left_ghost, top_ghost) = cell_gradient(:, :, right, bottom)
    !   cell_gradient(:, :, left_ghost, bottom_ghost) = cell_gradient(:, :, right, top)
    !   cell_gradient(:, :, left_ghost, bottom:top) = cell_gradient(:, :, right, bottom:top)
    ! case('+y')
    !   cell_gradient(:, :, left_ghost, top_ghost) = cell_gradient(:, :, right, bottom)
    !   cell_gradient(:, :, right_ghost, top_ghost) = cell_gradient(:, :, left, bottom)
    !   cell_gradient(:, :, left:right, top_ghost) = cell_gradient(:, :, left:right, bottom)
    ! case('-y')
    !   cell_gradient(:, :, left_ghost, bottom_ghost) = cell_gradient(:, :, right, top)
    !   cell_gradient(:, :, right_ghost, bottom_ghost) = cell_gradient(:, :, left, top)
    !   cell_gradient(:, :, left:right, bottom_ghost) = cell_gradient(:, :, left:right, top)
    ! case default
    !   error stop "Unsupported location to apply the bc at in reflection_bc_t%apply_reflection_cell_gradient_bc()"
    ! end select

  end subroutine apply_reflection_cell_gradient_bc

end module mod_reflection_bc
