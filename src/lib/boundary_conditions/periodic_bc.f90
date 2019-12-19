module mod_periodic_bc

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_input, only: input_t

  implicit none

  private
  public :: periodic_bc_t, periodic_bc_constructor

  type, extends(boundary_condition_t) :: periodic_bc_t
  contains
    ! procedure, public :: initialize => init_periodic_bc
    procedure, public :: apply_conserved_var_bc => apply_periodic_conserved_var_bc
    procedure, public :: apply_reconstructed_state_bc => apply_periodic_reconstructed_state_bc
    procedure, public :: apply_cell_gradient_bc => apply_periodic_cell_gradient_bc
    procedure, private, nopass :: apply_plus_x_conserved_vars
    procedure, private, nopass :: apply_plus_y_conserved_vars
    procedure, private, nopass :: apply_minus_x_conserved_vars
    procedure, private, nopass :: apply_minus_y_conserved_vars
    procedure, private, nopass :: apply_plus_x_reconstructed_state
    procedure, private, nopass :: apply_plus_y_reconstructed_state
    procedure, private, nopass :: apply_minus_x_reconstructed_state
    procedure, private, nopass :: apply_minus_y_reconstructed_state
    procedure, private, nopass :: apply_plus_x_cell_gradient
    procedure, private, nopass :: apply_plus_y_cell_gradient
    procedure, private, nopass :: apply_minus_x_cell_gradient
    procedure, private, nopass :: apply_minus_y_cell_gradient
    procedure, public :: copy => copy_periodic_bc

  end type periodic_bc_t

contains

  function periodic_bc_constructor(location) result(bc)
    type(periodic_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    allocate(bc)
    bc%name = 'periodic'
    bc%location = location
  end function periodic_bc_constructor

  subroutine copy_periodic_bc(out_bc, in_bc)
    class(boundary_condition_t), intent(in) :: in_bc
    class(periodic_bc_t), intent(inout) :: out_bc

    call debug_print('Calling boundary_condition_t%copy_periodic_bc()', __FILE__, __LINE__)
    ! if (allocated(out_bc%name)) deallocate(out_bc%name)
    ! allocate(out_bc%name, source=in_bc%name)
    out_bc%name = in_bc%name
    out_bc%location = in_bc%location
  end subroutine

  subroutine apply_periodic_conserved_var_bc(self, conserved_vars)
    !< Apply periodic boundary conditions to the conserved state vector field

    class(periodic_bc_t), intent(in) :: self
    real(rk), dimension(:, 0:, 0:), intent(inout) :: conserved_vars
    !< ((rho, u ,v, p), i, j); Conserved variables for each cell

    select case(self%location)
    case('+x')
      call self%apply_plus_x_conserved_vars(conserved_vars)
    case('-x')
      call self%apply_minus_x_conserved_vars(conserved_vars)
    case('+y')
      call self%apply_plus_y_conserved_vars(conserved_vars)
    case('-y')
      call self%apply_minus_y_conserved_vars(conserved_vars)
    case default
      error stop "Unsupported location to apply the bc at in periodic_bc_t%apply_periodic_conserved_var_bc()"
    end select

  end subroutine apply_periodic_conserved_var_bc

  subroutine apply_periodic_reconstructed_state_bc(self, reconstructed_state)
    !< Apply periodic boundary conditions to the reconstructed state vector field

    class(periodic_bc_t), intent(in) :: self
    real(rk), dimension(:, :, :, 0:, 0:), intent(inout) :: reconstructed_state
    !< ((rho, u ,v, p), point, node/midpoint, i, j); Reconstructed state for each cell

    select case(self%location)
    case('+x')
      call self%apply_plus_x_reconstructed_state(reconstructed_state)
    case('-x')
      call self%apply_minus_x_reconstructed_state(reconstructed_state)
    case('+y')
      call self%apply_plus_y_reconstructed_state(reconstructed_state)
    case('-y')
      call self%apply_minus_y_reconstructed_state(reconstructed_state)
    case default
      error stop "Unsupported location to apply the bc at in periodic_bc_t%apply_periodic_reconstructed_state_bc()"
    end select

  end subroutine apply_periodic_reconstructed_state_bc

  subroutine apply_periodic_cell_gradient_bc(self, cell_gradient)
    !< Apply periodic boundary conditions to the reconstructed state vector field

    class(periodic_bc_t), intent(in) :: self
    real(rk), dimension(:, :, 0:, 0:), intent(inout) :: cell_gradient
    !< ((rho, u ,v, p), point, node/midpoint, i, j); Reconstructed state for each cell

    select case(self%location)
    case('+x')
      call self%apply_plus_x_cell_gradient(cell_gradient)
    case('-x')
      call self%apply_minus_x_cell_gradient(cell_gradient)
    case('+y')
      call self%apply_plus_y_cell_gradient(cell_gradient)
    case('-y')
      call self%apply_minus_y_cell_gradient(cell_gradient)
    case default
      error stop "Unsupported location to apply the bc at in periodic_bc_t%apply_periodic_cell_gradient_bc()"
    end select

  end subroutine apply_periodic_cell_gradient_bc

  pure subroutine apply_plus_x_conserved_vars(conserved_vars)
    !< Implementation of applying the +x boundary condition for the conserved variable state

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

    conserved_vars(:, right_ghost, top_ghost) = conserved_vars(:, left, bottom)
    conserved_vars(:, right_ghost, bottom_ghost) = conserved_vars(:, left, top)
    conserved_vars(:, right_ghost, bottom:top) = conserved_vars(:, left, bottom:top)

  end subroutine apply_plus_x_conserved_vars

  pure subroutine apply_plus_y_conserved_vars(conserved_vars)
    !< Implementation of applying the +y boundary condition for the conserved variable state

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

    conserved_vars(:, left_ghost, top_ghost) = conserved_vars(:, right, bottom)
    conserved_vars(:, right_ghost, top_ghost) = conserved_vars(:, left, bottom)
    conserved_vars(:, left:right, top_ghost) = conserved_vars(:, left:right, bottom)

  end subroutine apply_plus_y_conserved_vars

  pure subroutine apply_minus_x_conserved_vars(conserved_vars)
    !< Implementation of applying the -x boundary condition for the conserved variable state

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

    conserved_vars(:, left_ghost, top_ghost) = conserved_vars(:, right, bottom)
    conserved_vars(:, left_ghost, bottom_ghost) = conserved_vars(:, right, top)
    conserved_vars(:, left_ghost, bottom:top) = conserved_vars(:, right, bottom:top)

  end subroutine apply_minus_x_conserved_vars

  pure subroutine apply_minus_y_conserved_vars(conserved_vars)
    !< Implementation of applying the -y boundary condition for the conserved variable state

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

    conserved_vars(:, left_ghost, bottom_ghost) = conserved_vars(:, right, top)
    conserved_vars(:, right_ghost, bottom_ghost) = conserved_vars(:, left, top)
    conserved_vars(:, left:right, bottom_ghost) = conserved_vars(:, left:right, top)

  end subroutine apply_minus_y_conserved_vars

  pure subroutine apply_plus_x_reconstructed_state(reconstructed_state)
    !< Implementation of applying the +x boundary condition for the cell gradient

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

    reconstructed_state(:, :, :, right_ghost, top_ghost) = reconstructed_state(:, :, :, left, bottom)
    reconstructed_state(:, :, :, right_ghost, bottom_ghost) = reconstructed_state(:, :, :, left, top)
    reconstructed_state(:, :, :, right_ghost, bottom:top) = reconstructed_state(:, :, :, left, bottom:top)

  end subroutine apply_plus_x_reconstructed_state

  pure subroutine apply_plus_y_reconstructed_state(reconstructed_state)
    !< Implementation of applying the +y boundary condition for the cell gradient

    real(rk), dimension(:, :, :, 0:, 0:), intent(inout) :: reconstructed_state

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

    reconstructed_state(:, :, :, left_ghost, top_ghost) = reconstructed_state(:, :, :, right, bottom)
    reconstructed_state(:, :, :, right_ghost, top_ghost) = reconstructed_state(:, :, :, left, bottom)
    reconstructed_state(:, :, :, left:right, top_ghost) = reconstructed_state(:, :, :, left:right, bottom)

  end subroutine apply_plus_y_reconstructed_state

  pure subroutine apply_minus_x_reconstructed_state(reconstructed_state)
    !< Implementation of applying the -x boundary condition for the cell gradient

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

    reconstructed_state(:, :, :, left_ghost, top_ghost) = reconstructed_state(:, :, :, right, bottom)
    reconstructed_state(:, :, :, left_ghost, bottom_ghost) = reconstructed_state(:, :, :, right, top)
    reconstructed_state(:, :, :, left_ghost, bottom:top) = reconstructed_state(:, :, :, right, bottom:top)

  end subroutine apply_minus_x_reconstructed_state

  pure subroutine apply_minus_y_reconstructed_state(reconstructed_state)
    !< Implementation of applying the -y boundary condition for the cell gradient

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

    reconstructed_state(:, :, :, left_ghost, bottom_ghost) = reconstructed_state(:, :, :, right, top)
    reconstructed_state(:, :, :, right_ghost, bottom_ghost) = reconstructed_state(:, :, :, left, top)
    reconstructed_state(:, :, :, left:right, bottom_ghost) = reconstructed_state(:, :, :, left:right, top)

  end subroutine apply_minus_y_reconstructed_state

  pure subroutine apply_plus_x_cell_gradient(cell_gradient)
    !< Implementation of applying the +x boundary condition for the cell gradient

    real(rk), dimension(:, :, 0:, 0:), intent(inout) :: cell_gradient
    !< ((rho, u ,v, p), point, node/midpoint, i, j); Reconstructed state for each cell

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

    cell_gradient(:, :, right_ghost, top_ghost) = cell_gradient(:, :, left, bottom)
    cell_gradient(:, :, right_ghost, bottom_ghost) = cell_gradient(:, :, left, top)
    cell_gradient(:, :, right_ghost, bottom:top) = cell_gradient(:, :, left, bottom:top)

  end subroutine apply_plus_x_cell_gradient

  pure subroutine apply_plus_y_cell_gradient(cell_gradient)
    !< Implementation of applying the +y boundary condition for the cell gradient

    real(rk), dimension(:, :, 0:, 0:), intent(inout) :: cell_gradient

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

    cell_gradient(:, :, left_ghost, top_ghost) = cell_gradient(:, :, right, bottom)
    cell_gradient(:, :, right_ghost, top_ghost) = cell_gradient(:, :, left, bottom)
    cell_gradient(:, :, left:right, top_ghost) = cell_gradient(:, :, left:right, bottom)

  end subroutine apply_plus_y_cell_gradient

  pure subroutine apply_minus_x_cell_gradient(cell_gradient)
    !< Implementation of applying the -x boundary condition for the cell gradient

    real(rk), dimension(:, :, 0:, 0:), intent(inout) :: cell_gradient
    !< ((rho, u ,v, p), point, node/midpoint, i, j); Reconstructed state for each cell

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

    cell_gradient(:, :, left_ghost, top_ghost) = cell_gradient(:, :, right, bottom)
    cell_gradient(:, :, left_ghost, bottom_ghost) = cell_gradient(:, :, right, top)
    cell_gradient(:, :, left_ghost, bottom:top) = cell_gradient(:, :, right, bottom:top)

  end subroutine apply_minus_x_cell_gradient

  pure subroutine apply_minus_y_cell_gradient(cell_gradient)
    !< Implementation of applying the -y boundary condition for the cell gradient

    real(rk), dimension(:, :, 0:, 0:), intent(inout) :: cell_gradient
    !< ((rho, u ,v, p), point, node/midpoint, i, j); Reconstructed state for each cell

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

    cell_gradient(:, :, left_ghost, bottom_ghost) = cell_gradient(:, :, right, top)
    cell_gradient(:, :, right_ghost, bottom_ghost) = cell_gradient(:, :, left, top)
    cell_gradient(:, :, left:right, bottom_ghost) = cell_gradient(:, :, left:right, top)

  end subroutine apply_minus_y_cell_gradient

end module mod_periodic_bc
