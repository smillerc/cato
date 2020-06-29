module mod_symmetry_bc
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_grid, only: grid_t
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_input, only: input_t

  implicit none

  private
  public :: symmetry_bc_t, symmetry_bc_constructor

  type, extends(boundary_condition_t) :: symmetry_bc_t
  contains
    procedure, public :: apply_primitive_var_bc => apply_symmetry_primitive_var_bc
    procedure, public :: apply_reconstructed_state_bc => apply_symmetry_reconstructed_state_bc
    procedure, public :: apply_gradient_bc
    final :: finalize
  end type
contains
  function symmetry_bc_constructor(location, input, grid) result(bc)
    type(symmetry_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input
    class(grid_t), intent(in) :: grid

    allocate(bc)
    bc%name = 'symmetry'
    bc%location = location
    call bc%set_indices(grid)
  end function symmetry_bc_constructor

  subroutine apply_symmetry_primitive_var_bc(self, rho, u, v, p, lbounds)

    class(symmetry_bc_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: rho
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: u
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: v
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: p

    integer(ik) :: i

    associate(left=>self%ilo, right=>self%ihi, bottom=>self%jlo, top=>self%jhi, &
              left_ghost=>self%ilo_ghost, right_ghost=>self%ihi_ghost, &
              bottom_ghost=>self%jlo_ghost, top_ghost=>self%jhi_ghost)

      select case(self%location)
      case('+x')
        call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() +x', __FILE__, __LINE__)

        do i = 1, self%n_ghost_layers
          rho(right_ghost(i), :) = rho(right + (i - 1), :)
          u(right_ghost(i), :) = -u(right + (i - 1), :)
          v(right_ghost(i), :) = v(right + (i - 1), :)
          p(right_ghost(i), :) = p(right + (i - 1), :)
        end do

      case('-x')
        call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() -x', __FILE__, __LINE__)

        do i = 1, self%n_ghost_layers
          rho(left_ghost(i), :) = rho(left + (i - 1), :)
          u(left_ghost(i), :) = -u(left + (i - 1), :)
          v(left_ghost(i), :) = v(left + (i - 1), :)
          p(left_ghost(i), :) = p(left + (i - 1), :)
        end do

      case('+y')
        call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() +y', __FILE__, __LINE__)

        do i = 1, self%n_ghost_layers
          rho(:, top_ghost(i)) = rho(:, top + (i - 1))
          u(:, top_ghost(i)) = u(:, top + (i - 1))
          v(:, top_ghost(i)) = -v(:, top + (i - 1))
          p(:, top_ghost(i)) = p(:, top + (i - 1))
        end do

      case('-y')
        call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() -y', __FILE__, __LINE__)
        do i = 1, self%n_ghost_layers
          rho(:, bottom_ghost(i)) = rho(:, bottom + (i - 1))
          u(:, bottom_ghost(i)) = u(:, bottom + (i - 1))
          v(:, bottom_ghost(i)) = -v(:, bottom + (i - 1))
          p(:, bottom_ghost(i)) = p(:, bottom + (i - 1))
        end do

      case default
        error stop "Unsupported location to apply the bc at in symmetry_bc_t%apply_symmetry_primitive_var_bc()"
      end select
    end associate

  end subroutine apply_symmetry_primitive_var_bc

  subroutine apply_symmetry_reconstructed_state_bc(self, recon_rho, recon_u, recon_v, recon_p, lbounds)
    !< Apply zero gradient boundary conditions to the reconstructed state vector field

    class(symmetry_bc_t), intent(inout) :: self

    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_rho
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_u
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_v
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_p

    integer(ik) :: i

    associate(left=>self%ilo, right=>self%ihi, bottom=>self%jlo, top=>self%jhi, &
              left_ghost=>self%ilo_ghost, right_ghost=>self%ihi_ghost, &
              bottom_ghost=>self%jlo_ghost, top_ghost=>self%jhi_ghost)

      select case(self%location)
      case('+x')
        call debug_print('Running periodic_bc_t%apply_periodic_reconstructed_state_bc() +x', __FILE__, __LINE__)

        do i = 1, self%n_ghost_layers
          recon_rho(:, right_ghost(i), :) = recon_rho(:, right + (i - 1), :)
          recon_u(:, right_ghost(i), :) = -recon_u(:, right + (i - 1), :)
          recon_v(:, right_ghost(i), :) = recon_v(:, right + (i - 1), :)
          recon_p(:, right_ghost(i), :) = recon_p(:, right + (i - 1), :)
        end do

      case('-x')
        call debug_print('Running periodic_bc_t%apply_periodic_reconstructed_state_bc() -x', __FILE__, __LINE__)

        do i = 1, self%n_ghost_layers
          recon_rho(:, left_ghost(i), :) = recon_rho(:, left + (i - 1), :)
          recon_u(:, left_ghost(i), :) = -recon_u(:, left + (i - 1), :)
          recon_v(:, left_ghost(i), :) = recon_v(:, left + (i - 1), :)
          recon_p(:, left_ghost(i), :) = recon_p(:, left + (i - 1), :)
        end do
      case('+y')
        call debug_print('Running periodic_bc_t%apply_periodic_reconstructed_state_bc() +y', __FILE__, __LINE__)

        do i = 1, self%n_ghost_layers
          recon_rho(:, :, top_ghost(i)) = recon_rho(:, :, top + (i - 1))
          recon_u(:, :, top_ghost(i)) = recon_u(:, :, top + (i - 1))
          recon_v(:, :, top_ghost(i)) = -recon_v(:, :, top + (i - 1))
          recon_p(:, :, top_ghost(i)) = recon_p(:, :, top + (i - 1))
        end do
      case('-y')
        call debug_print('Running periodic_bc_t%apply_periodic_reconstructed_state_bc() -y', __FILE__, __LINE__)

        do i = 1, self%n_ghost_layers
          recon_rho(:, :, bottom_ghost(i)) = recon_rho(:, :, bottom + (i - 1))
          recon_u(:, :, bottom_ghost(i)) = recon_u(:, :, bottom + (i - 1))
          recon_v(:, :, bottom_ghost(i)) = -recon_v(:, :, bottom + (i - 1))
          recon_p(:, :, bottom_ghost(i)) = recon_p(:, :, bottom + (i - 1))
        end do
      case default
        error stop "Unsupported location to apply the bc at in periodic_bc_t%apply_periodic_reconstructed_state_bc()"
      end select
    end associate

  end subroutine apply_symmetry_reconstructed_state_bc

  subroutine apply_gradient_bc(self, grad_x, grad_y, lbounds)
    class(symmetry_bc_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: grad_x
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: grad_y
    integer(ik) :: i

    associate(left=>self%ilo, right=>self%ihi, bottom=>self%jlo, top=>self%jhi, &
              left_ghost=>self%ilo_ghost, right_ghost=>self%ihi_ghost, &
              bottom_ghost=>self%jlo_ghost, top_ghost=>self%jhi_ghost)

      select case(self%location)
      case('+x')
        call debug_print('Running periodic_bc_t%apply_gradient_bc() +x', __FILE__, __LINE__)
        do i = 1, self%n_ghost_layers
          grad_x(right_ghost(i), :) = -grad_x(left + (i - 1), :)
          grad_y(right_ghost(i), :) = grad_y(left + (i - 1), :)
        end do
      case('-x')
        call debug_print('Running periodic_bc_t%apply_gradient_bc() -x', __FILE__, __LINE__)
        do i = 1, self%n_ghost_layers
          grad_x(left_ghost(i), :) = -grad_x(right + (i - 1), :)
          grad_y(left_ghost(i), :) = grad_y(right + (i - 1), :)
        end do
      case('+y')
        call debug_print('Running periodic_bc_t%apply_gradient_bc() +y', __FILE__, __LINE__)
        do i = 1, self%n_ghost_layers
          grad_x(:, top_ghost(i)) = grad_x(:, bottom + (i - 1))
          grad_y(:, top_ghost(i)) = -grad_y(:, bottom + (i - 1))
        end do
      case('-y')
        call debug_print('Running periodic_bc_t%apply_gradient_bc() -y', __FILE__, __LINE__)
        do i = 1, self%n_ghost_layers
          grad_x(:, bottom_ghost(i)) = grad_x(:, top + (i - 1))
          grad_y(:, bottom_ghost(i)) = -grad_y(:, top + (i - 1))
        end do
      case default
        error stop "Unsupported location to apply the bc at in periodic_bc_t%apply_gradient_bc()"
      end select
    end associate
  end subroutine apply_gradient_bc

  subroutine finalize(self)
    type(symmetry_bc_t), intent(inout) :: self
    call debug_print('Running symmetry_bc_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%ilo_ghost)) deallocate(self%ilo_ghost)
    if(allocated(self%ihi_ghost)) deallocate(self%ihi_ghost)
    if(allocated(self%jlo_ghost)) deallocate(self%jlo_ghost)
    if(allocated(self%jhi_ghost)) deallocate(self%jhi_ghost)
  end subroutine
end module mod_symmetry_bc
