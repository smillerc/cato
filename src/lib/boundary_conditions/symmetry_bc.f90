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
    procedure, public :: apply_primitive_var_bc => apply_symmetry_primitive_var_bc
    procedure, public :: apply_reconstructed_state_bc => apply_symmetry_reconstructed_state_bc
    procedure, public :: copy => copy_symmetry_bc
  end type
contains
  function symmetry_bc_constructor(location, input, ghost_layers) result(bc)
    type(symmetry_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input
    integer(ik), dimension(:, :), intent(in) :: ghost_layers
    !< (ilo_layers(n), ihi_layers(n), jlo_layers(n), jhi_layers(n)); indices to the ghost layers.
    !< The ilo_layers type var can be scalar or a vector, e.g. ilo_layers = [-1,0] or ilo_layers = 0

    allocate(bc)
    bc%name = 'symmetry'
    bc%location = location
    call bc%set_indices(ghost_layers)
  end function symmetry_bc_constructor

  subroutine copy_symmetry_bc(out_bc, in_bc)
    class(boundary_condition_t), intent(in) :: in_bc
    class(symmetry_bc_t), intent(inout) :: out_bc
  end subroutine

  subroutine apply_symmetry_primitive_var_bc(self, rho, u, v, p, lbounds)

    class(symmetry_bc_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: rho
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: u
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: v
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: p

    integer(ik) :: left         !< Min i real cell index
    integer(ik) :: right        !< Max i real cell index
    integer(ik) :: bottom       !< Min j real cell index
    integer(ik) :: top          !< Max j real cell index
    integer(ik) :: left_ghost   !< Min i ghost cell index
    integer(ik) :: right_ghost  !< Max i ghost cell index
    integer(ik) :: bottom_ghost !< Min j ghost cell index
    integer(ik) :: top_ghost    !< Max j ghost cell index

    left_ghost = lbound(rho, dim=1)
    right_ghost = ubound(rho, dim=1)
    bottom_ghost = lbound(rho, dim=2)
    top_ghost = ubound(rho, dim=2)
    left = left_ghost + 1
    right = right_ghost - 1
    bottom = bottom_ghost + 1
    top = top_ghost - 1

    select case(self%location)
    case('+x')
      call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() +x', __FILE__, __LINE__)
      rho(right_ghost, :) = rho(right, :)
      u(right_ghost, :) = -u(right, :)
      v(right_ghost, :) = v(right, :)
      p(right_ghost, :) = p(right, :)
    case('-x')
      call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() -x', __FILE__, __LINE__)
      rho(left_ghost, :) = rho(left, :)
      u(left_ghost, :) = -u(left, :)
      v(left_ghost, :) = v(left, :)
      p(left_ghost, :) = p(left, :)
    case('+y')
      call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() +y', __FILE__, __LINE__)
      rho(:, top_ghost) = rho(:, top)
      u(:, top_ghost) = u(:, top)
      v(:, top_ghost) = -v(:, top)
      p(:, top_ghost) = p(:, top)
    case('-y')
      call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() -y', __FILE__, __LINE__)
      rho(:, bottom_ghost) = rho(:, bottom)
      u(:, bottom_ghost) = u(:, bottom)
      v(:, bottom_ghost) = -v(:, bottom)
      p(:, bottom_ghost) = p(:, bottom)
    case default
      error stop "Unsupported location to apply the bc at in symmetry_bc_t%apply_symmetry_cell_gradient_bc()"
    end select

  end subroutine apply_symmetry_primitive_var_bc

  subroutine apply_symmetry_reconstructed_state_bc(self, recon_rho, recon_u, recon_v, recon_p, lbounds)
    !< Apply zero gradient boundary conditions to the reconstructed state vector field

    class(symmetry_bc_t), intent(inout) :: self

    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_rho
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_u
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_v
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_p

    integer(ik) :: left         !< Min i real cell index
    integer(ik) :: right        !< Max i real cell index
    integer(ik) :: bottom       !< Min j real cell index
    integer(ik) :: top          !< Max j real cell index
    integer(ik) :: left_ghost   !< Min i ghost cell index
    integer(ik) :: right_ghost  !< Max i ghost cell index
    integer(ik) :: bottom_ghost !< Min j ghost cell index
    integer(ik) :: top_ghost    !< Max j ghost cell index

    left_ghost = lbound(recon_rho, dim=2)
    right_ghost = ubound(recon_rho, dim=2)
    bottom_ghost = lbound(recon_rho, dim=3)
    top_ghost = ubound(recon_rho, dim=3)
    left = left_ghost + 1
    right = right_ghost - 1
    bottom = bottom_ghost + 1
    top = top_ghost - 1

    select case(self%location)
    case('+x')
      call debug_print('Running symmetry_bc_t%apply_symmetry_reconstructed_state_bc() +x', __FILE__, __LINE__)
      recon_rho(:, right_ghost, :) = recon_rho(:, right, :)
      recon_u(:, right_ghost, :) = -recon_u(:, right, :)
      recon_v(:, right_ghost, :) = recon_v(:, right, :)
      recon_p(:, right_ghost, :) = recon_p(:, right, :)
    case('-x')
      call debug_print('Running symmetry_bc_t%apply_symmetry_reconstructed_state_bc() -x', __FILE__, __LINE__)
      recon_rho(:, left_ghost, :) = recon_rho(:, left, :)
      recon_u(:, left_ghost, :) = -recon_u(:, left, :)
      recon_v(:, left_ghost, :) = recon_v(:, left, :)
      recon_p(:, left_ghost, :) = recon_p(:, left, :)
    case('+y')
      call debug_print('Running symmetry_bc_t%apply_symmetry_reconstructed_state_bc() +y', __FILE__, __LINE__)
      recon_rho(:, :, top_ghost) = recon_rho(:, :, top)
      recon_u(:, :, top_ghost) = recon_u(:, :, top)
      recon_v(:, :, top_ghost) = -recon_v(:, :, top)
      recon_p(:, :, top_ghost) = recon_p(:, :, top)
    case('-y')
      call debug_print('Running symmetry_bc_t%apply_symmetry_reconstructed_state_bc() -y', __FILE__, __LINE__)
      recon_rho(:, :, bottom_ghost) = recon_rho(:, :, bottom)
      recon_u(:, :, bottom_ghost) = recon_u(:, :, bottom)
      recon_v(:, :, bottom_ghost) = -recon_v(:, :, bottom)
      recon_p(:, :, bottom_ghost) = recon_p(:, :, bottom)
    case default
      error stop "Unsupported location to apply the bc at in symmetry_bc_t%apply_symmetry_reconstructed_state_bc()"
    end select

  end subroutine apply_symmetry_reconstructed_state_bc

end module mod_symmetry_bc
