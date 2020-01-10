module mod_symmetry_bc
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_input, only: input_t

  implicit none

  private
  public :: symmetry_bc_t, symmetry_bc_constructor

  type, extends(boundary_condition_t) :: symmetry_bc_t
  contains
    procedure, public :: apply_conserved_var_bc => apply_symmetry_conserved_var_bc
    procedure, public :: apply_reconstructed_state_bc => apply_symmetry_reconstructed_state_bc
    procedure, public :: apply_cell_gradient_bc => apply_symmetry_cell_gradient_bc
    procedure, public :: copy => copy_symmetry_bc
  end type

contains

  function symmetry_bc_constructor(location, input) result(bc)
    type(symmetry_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input

    allocate(bc)
    bc%name = 'symmetry'
    bc%location = location
  end function symmetry_bc_constructor

  subroutine copy_symmetry_bc(out_bc, in_bc)
    class(boundary_condition_t), intent(in) :: in_bc
    class(symmetry_bc_t), intent(inout) :: out_bc

    call debug_print('Calling symmetry_bc_t%copy_symmetry_bc()', __FILE__, __LINE__)
    out_bc%name = in_bc%name
    out_bc%location = in_bc%location

  end subroutine

  subroutine apply_symmetry_conserved_var_bc(self, conserved_vars)
    !< Apply symmetry boundary conditions to the conserved state vector field
    class(symmetry_bc_t), intent(inout) :: self
    real(rk), dimension(:, 0:, 0:), intent(inout) :: conserved_vars
    !< ((rho, u ,v, p), i, j); Conserved variables for each cell

    integer(ik) :: left         !< Min i real cell index
    integer(ik) :: right        !< Max i real cell index
    integer(ik) :: bottom       !< Min j real cell index
    integer(ik) :: top          !< Max j real cell index
    integer(ik) :: left_ghost   !< Min i ghost cell index
    integer(ik) :: right_ghost  !< Max i ghost cell index
    integer(ik) :: bottom_ghost !< Min j ghost cell index
    integer(ik) :: top_ghost    !< Max j ghost cell index

    left_ghost = lbound(conserved_vars, dim=2)
    right_ghost = ubound(conserved_vars, dim=2)
    bottom_ghost = lbound(conserved_vars, dim=3)
    top_ghost = ubound(conserved_vars, dim=3)
    left = left_ghost + 1
    right = right_ghost - 1
    bottom = bottom_ghost + 1
    top = top_ghost - 1

    select case(self%location)
    case('+x')
      conserved_vars(1, right_ghost, :) = conserved_vars(1, right, :)      ! density
      conserved_vars(2, right_ghost, :) = -1 * conserved_vars(2, right, :) ! x velocity
      conserved_vars(3, right_ghost, :) = conserved_vars(3, right, :)      ! y velocity
      conserved_vars(4, right_ghost, :) = conserved_vars(4, right, :)      ! pressure
    case('-x')
      conserved_vars(1, left_ghost, :) = conserved_vars(1, left, :)      ! density
      conserved_vars(2, left_ghost, :) = -1 * conserved_vars(2, left, :) ! x velocity
      conserved_vars(3, left_ghost, :) = conserved_vars(3, left, :)      ! y velocity
      conserved_vars(4, left_ghost, :) = conserved_vars(4, left, :)      ! pressure
    case('+y')
      conserved_vars(1, :, top_ghost) = conserved_vars(1, :, top)      ! density
      conserved_vars(2, :, top_ghost) = conserved_vars(2, :, top)      ! x velocity
      conserved_vars(3, :, top_ghost) = -1 * conserved_vars(3, :, top) ! y velocity
      conserved_vars(4, :, top_ghost) = conserved_vars(4, :, top)      ! pressure
    case('-y')
      conserved_vars(1, :, bottom_ghost) = conserved_vars(1, :, bottom)      ! density
      conserved_vars(2, :, bottom_ghost) = conserved_vars(2, :, bottom)      ! x velocity
      conserved_vars(3, :, bottom_ghost) = -1 * conserved_vars(3, :, bottom) ! y velocity
      conserved_vars(4, :, bottom_ghost) = conserved_vars(4, :, bottom)      ! pressure
    case default
      error stop "Unsupported location to apply the bc at in symmetry_bc_t%apply_symmetry_cell_gradient_bc()"
    end select

  end subroutine apply_symmetry_conserved_var_bc

  subroutine apply_symmetry_reconstructed_state_bc(self, reconstructed_state)
    !< Apply symmetry boundary conditions to the reconstructed state vector field

    class(symmetry_bc_t), intent(in) :: self
    real(rk), dimension(:, :, :, 0:, 0:), intent(inout) :: reconstructed_state
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

    select case(self%location)
    case('+x')
      reconstructed_state(1, :, :, right_ghost, :) = reconstructed_state(1, :, :, right, :)
      reconstructed_state(2, :, :, right_ghost, :) = -1 * reconstructed_state(2, :, :, right, :)
      reconstructed_state(3, :, :, right_ghost, :) = reconstructed_state(3, :, :, right, :)
      reconstructed_state(4, :, :, right_ghost, :) = reconstructed_state(4, :, :, right, :)
    case('-x')
      reconstructed_state(1, :, :, left_ghost, :) = reconstructed_state(1, :, :, left, :)
      reconstructed_state(2, :, :, left_ghost, :) = -1 * reconstructed_state(2, :, :, left, :)
      reconstructed_state(3, :, :, left_ghost, :) = reconstructed_state(3, :, :, left, :)
      reconstructed_state(4, :, :, left_ghost, :) = reconstructed_state(4, :, :, left, :)
    case('+y')
      reconstructed_state(1, :, :, :, top_ghost) = reconstructed_state(1, :, :, :, top)
      reconstructed_state(2, :, :, :, top_ghost) = reconstructed_state(2, :, :, :, top)
      reconstructed_state(3, :, :, :, top_ghost) = -1 * reconstructed_state(3, :, :, :, top)
      reconstructed_state(4, :, :, :, top_ghost) = reconstructed_state(4, :, :, :, top)
    case('-y')
      reconstructed_state(1, :, :, :, bottom_ghost) = reconstructed_state(1, :, :, :, bottom)
      reconstructed_state(2, :, :, :, bottom_ghost) = reconstructed_state(2, :, :, :, bottom)
      reconstructed_state(3, :, :, :, bottom_ghost) = -1 * reconstructed_state(3, :, :, :, bottom)
      reconstructed_state(4, :, :, :, bottom_ghost) = reconstructed_state(4, :, :, :, bottom)
    case default
      error stop "Unsupported location to apply the bc at in symmetry_bc_t%apply_symmetry_reconstructed_state_bc()"
    end select

  end subroutine apply_symmetry_reconstructed_state_bc

  subroutine apply_symmetry_cell_gradient_bc(self, cell_gradient)
    !< Apply symmetry boundary conditions to the cell gradient state vector field

    class(symmetry_bc_t), intent(in) :: self
    real(rk), dimension(:, :, 0:, 0:), intent(inout) :: cell_gradient
    !< ((d/dx, d/dy), (rho, u ,v, p), i, j); Gradient of each cell's conserved quantities

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
      cell_gradient(1, :, right_ghost, :) = -1 * cell_gradient(1, :, right, :)
      cell_gradient(2, :, right_ghost, :) = cell_gradient(2, :, right, :)
    case('-x')
      cell_gradient(1, :, left_ghost, :) = -1 * cell_gradient(1, :, left, :)
      cell_gradient(2, :, left_ghost, :) = cell_gradient(2, :, left, :)
    case('+y')
      cell_gradient(1, :, :, top_ghost) = cell_gradient(1, :, :, top)
      cell_gradient(2, :, :, top_ghost) = -1 * cell_gradient(2, :, :, top)
    case('-y')
      cell_gradient(1, :, :, bottom_ghost) = cell_gradient(1, :, :, bottom)
      cell_gradient(2, :, :, bottom_ghost) = -1 * cell_gradient(2, :, :, bottom)
    case default
      error stop "Unsupported location to apply the bc at in symmetry_bc_t%apply_symmetry_cell_gradient_bc()"
    end select

  end subroutine apply_symmetry_cell_gradient_bc

end module mod_symmetry_bc
