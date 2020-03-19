module mod_vacuum_bc
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_input, only: input_t

  implicit none

  private
  public :: vacuum_bc_t, vacuum_bc_constructor

  type, extends(boundary_condition_t) :: vacuum_bc_t
    real(rk) :: vacuum_density = 1e-3_rk
    real(rk) :: vacuum_pressure = 1e8 ! approx 1 atmosphere
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

      ! min_density = min(minval(primitive_vars(1, right-5:right, bottom:top)), self%vacuum_density)
      ! self%vacuum_density = min_density
      ! if (self%vacuum_density < 1e-4_rk) self%vacuum_density = 1e-4_rk

      min_pressure = min(minval(primitive_vars(4, right - 5:right, bottom:top)), self%vacuum_pressure)
      self%vacuum_pressure = min_pressure
      ! write(*,'(a, 1(es12.3))') "Vacuum pressure: ", self%vacuum_pressure

      if(self%vacuum_pressure < 1e8_rk) then
        self%vacuum_pressure = 1e8_rk
        primitive_vars(1, right:right_ghost, :) = self%vacuum_density
        primitive_vars(2:3, right_ghost, :) = primitive_vars(2:3, right, :)
        primitive_vars(4, right:right_ghost, :) = self%vacuum_pressure
      else
        primitive_vars(1, right_ghost, :) = primitive_vars(1, right, :)
        primitive_vars(2:3, right_ghost, :) = primitive_vars(2:3, right, :)
        primitive_vars(4, right_ghost, :) = self%vacuum_pressure
      end if

      ! write(*,'(a, 2(es12.3))') "Vacuum density, pressure: ", self%vacuum_density, self%vacuum_pressure
      ! primitive_vars(1, right_ghost, :) = self%vacuum_density
      primitive_vars(1, right_ghost, :) = primitive_vars(1, right, :)
      primitive_vars(2:3, right_ghost, :) = primitive_vars(2:3, right, :)
      primitive_vars(4, right_ghost, :) = self%vacuum_pressure
      ! case('-x')
      !   call debug_print('Running vacuum_bc_t%apply_vacuum_primitive_var_bc() -x', __FILE__, __LINE__)
      !   primitive_vars(:, left_ghost, :) = primitive_vars(:, left, :)
      ! case('+y')
      !   call debug_print('Running vacuum_bc_t%apply_vacuum_primitive_var_bc() +y', __FILE__, __LINE__)
      !   primitive_vars(:, :, top_ghost) = primitive_vars(:, :, top)
      ! case('-y')
      !   call debug_print('Running vacuum_bc_t%apply_vacuum_primitive_var_bc() -y', __FILE__, __LINE__)
      !   primitive_vars(:, :, bottom_ghost) = primitive_vars(:, :, bottom)
    case default
      error stop "Unsupported location to apply the bc at in vacuum_bc_t%apply_vacuum_cell_gradient_bc()"
    end select

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
      ! reconstructed_state(1, :, :, right_ghost, :) = self%vacuum_density
      if(self%vacuum_pressure < 1e-8_rk) then
        reconstructed_state(:, :, :, right_ghost, :) = reconstructed_state(:, :, :, right, :)
      else
        reconstructed_state(1, :, :, right_ghost, :) = reconstructed_state(1, :, :, right, :)
        reconstructed_state(2:3, :, :, right_ghost, :) = reconstructed_state(2:3, :, :, right, :)
        reconstructed_state(4, :, :, right_ghost, :) = self%vacuum_pressure
      end if

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
