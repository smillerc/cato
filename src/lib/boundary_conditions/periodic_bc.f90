module mod_periodic_bc

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_input, only: input_t

  implicit none

  private
  public :: periodic_bc_t, periodic_bc_constructor

  type, extends(boundary_condition_t) :: periodic_bc_t
    logical :: do_corners = .false.
  contains
    ! procedure, public :: initialize => init_periodic_bc
    procedure, public :: apply_primitive_var_bc => apply_periodic_primitive_var_bc
    procedure, public :: apply_reconstructed_state_bc => apply_periodic_reconstructed_state_bc
    procedure, public :: apply_gradient_bc
  end type periodic_bc_t

contains

  function periodic_bc_constructor(location, input, ghost_layers) result(bc)
    type(periodic_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input
    integer(ik), dimension(:, :), intent(in) :: ghost_layers
    !< (ilo_layers(n), ihi_layers(n), jlo_layers(n), jhi_layers(n)); indices to the ghost layers.
    !< The ilo_layers type var can be scalar or a vector, e.g. ilo_layers = [-1,0] or ilo_layers = 0

    allocate(bc)
    bc%name = 'periodic'
    bc%location = location
    call bc%set_indices(ghost_layers)

    if(trim(input%plus_x_bc) == 'periodic' .and. &
       trim(input%minus_x_bc) == 'periodic' .and. &
       trim(input%plus_y_bc) == 'periodic' .and. &
       trim(input%minus_y_bc) == 'periodic') then

      bc%do_corners = .true.
    end if

  end function periodic_bc_constructor

  subroutine apply_periodic_primitive_var_bc(self, rho, u, v, p, lbounds)

    class(periodic_bc_t), intent(inout) :: self
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
        call debug_print('Running periodic_bc_t%apply_periodic_primitive_var_bc() +x', __FILE__, __LINE__)
        if(self%do_corners) then
          do i = 1, self%n_ghost_layers
            rho(right_ghost(i), top_ghost(i)) = rho(left + (i - 1), bottom + (i - 1))
            u(right_ghost(i), top_ghost(i)) = u(left + (i - 1), bottom + (i - 1))
            v(right_ghost(i), top_ghost(i)) = v(left + (i - 1), bottom + (i - 1))
            p(right_ghost(i), top_ghost(i)) = p(left + (i - 1), bottom + (i - 1))
            rho(right_ghost(i), bottom_ghost(i)) = rho(left + (i - 1), top + (i - 1))
            u(right_ghost(i), bottom_ghost(i)) = u(left + (i - 1), top + (i - 1))
            v(right_ghost(i), bottom_ghost(i)) = v(left + (i - 1), top + (i - 1))
            p(right_ghost(i), bottom_ghost(i)) = p(left + (i - 1), top + (i - 1))

            rho(right_ghost(i), bottom:top) = rho(left + (i - 1), bottom:top)
            u(right_ghost(i), bottom:top) = u(left + (i - 1), bottom:top)
            v(right_ghost(i), bottom:top) = v(left + (i - 1), bottom:top)
            p(right_ghost(i), bottom:top) = p(left + (i - 1), bottom:top)
          end do

        else
          do i = 1, self%n_ghost_layers
            rho(right_ghost(i), :) = rho(left + (i - 1), :)
            u(right_ghost(i), :) = u(left + (i - 1), :)
            v(right_ghost(i), :) = v(left + (i - 1), :)
            p(right_ghost(i), :) = p(left + (i - 1), :)
          end do
        end if

      case('-x')
        call debug_print('Running periodic_bc_t%apply_periodic_primitive_var_bc() -x', __FILE__, __LINE__)
        if(self%do_corners) then
          do i = 1, self%n_ghost_layers
            rho(left_ghost(i), top_ghost(i)) = rho(right + (i - 1), bottom + (i - 1))
            u(left_ghost(i), top_ghost(i)) = u(right + (i - 1), bottom + (i - 1))
            v(left_ghost(i), top_ghost(i)) = v(right + (i - 1), bottom + (i - 1))
            p(left_ghost(i), top_ghost(i)) = p(right + (i - 1), bottom + (i - 1))
            rho(left_ghost(i), bottom_ghost(i)) = rho(right + (i - 1), top + (i - 1))
            u(left_ghost(i), bottom_ghost(i)) = u(right + (i - 1), top + (i - 1))
            v(left_ghost(i), bottom_ghost(i)) = v(right + (i - 1), top + (i - 1))
            p(left_ghost(i), bottom_ghost(i)) = p(right + (i - 1), top + (i - 1))

            rho(left_ghost(i), bottom:top) = rho(right + (i - 1), bottom:top)
            u(left_ghost(i), bottom:top) = u(right + (i - 1), bottom:top)
            v(left_ghost(i), bottom:top) = v(right + (i - 1), bottom:top)
            p(left_ghost(i), bottom:top) = p(right + (i - 1), bottom:top)
          end do

        else
          do i = 1, self%n_ghost_layers
            rho(left_ghost(i), :) = rho(right + (i - 1), :)
            u(left_ghost(i), :) = u(right + (i - 1), :)
            v(left_ghost(i), :) = v(right + (i - 1), :)
            p(left_ghost(i), :) = p(right + (i - 1), :)
          end do
        end if

      case('+y')
        call debug_print('Running periodic_bc_t%apply_periodic_primitive_var_bc() +y', __FILE__, __LINE__)
        if(self%do_corners) then
          do i = 1, self%n_ghost_layers
            rho(left_ghost(i), top_ghost(i)) = rho(right + (i - 1), bottom + (i - 1))
            u(left_ghost(i), top_ghost(i)) = u(right + (i - 1), bottom + (i - 1))
            v(left_ghost(i), top_ghost(i)) = v(right + (i - 1), bottom + (i - 1))
            p(left_ghost(i), top_ghost(i)) = p(right + (i - 1), bottom + (i - 1))
            rho(right_ghost(i), top_ghost(i)) = rho(left + (i - 1), bottom + (i - 1))
            u(right_ghost(i), top_ghost(i)) = u(left + (i - 1), bottom + (i - 1))
            v(right_ghost(i), top_ghost(i)) = v(left + (i - 1), bottom + (i - 1))
            p(right_ghost(i), top_ghost(i)) = p(left + (i - 1), bottom + (i - 1))

            rho(left:right, top_ghost(i)) = rho(left:right, bottom + (i - 1))
            u(left:right, top_ghost(i)) = u(left:right, bottom + (i - 1))
            v(left:right, top_ghost(i)) = v(left:right, bottom + (i - 1))
            p(left:right, top_ghost(i)) = p(left:right, bottom + (i - 1))
          end do

        else
          do i = 1, self%n_ghost_layers
            rho(:, top_ghost(i)) = rho(:, bottom + (i - 1))
            u(:, top_ghost(i)) = u(:, bottom + (i - 1))
            v(:, top_ghost(i)) = v(:, bottom + (i - 1))
            p(:, top_ghost(i)) = p(:, bottom + (i - 1))
          end do
        end if

      case('-y')
        call debug_print('Running periodic_bc_t%apply_periodic_primitive_var_bc() -y', __FILE__, __LINE__)
        if(self%do_corners) then
          do i = 1, self%n_ghost_layers
            rho(left_ghost(i), bottom_ghost(i)) = rho(right + (i - 1), top + (i - 1))
            u(left_ghost(i), bottom_ghost(i)) = u(right + (i - 1), top + (i - 1))
            v(left_ghost(i), bottom_ghost(i)) = v(right + (i - 1), top + (i - 1))
            p(left_ghost(i), bottom_ghost(i)) = p(right + (i - 1), top + (i - 1))
            rho(right_ghost(i), bottom_ghost(i)) = rho(left + (i - 1), top + (i - 1))
            u(right_ghost(i), bottom_ghost(i)) = u(left + (i - 1), top + (i - 1))
            v(right_ghost(i), bottom_ghost(i)) = v(left + (i - 1), top + (i - 1))
            p(right_ghost(i), bottom_ghost(i)) = p(left + (i - 1), top + (i - 1))
            rho(left:right, bottom_ghost(i)) = rho(left:right, top + (i - 1))
            u(left:right, bottom_ghost(i)) = u(left:right, top + (i - 1))
            v(left:right, bottom_ghost(i)) = v(left:right, top + (i - 1))
            p(left:right, bottom_ghost(i)) = p(left:right, top + (i - 1))
          end do
        else
          do i = 1, self%n_ghost_layers
            rho(:, bottom_ghost(i)) = rho(:, top + (i - 1))
            u(:, bottom_ghost(i)) = u(:, top + (i - 1))
            v(:, bottom_ghost(i)) = v(:, top + (i - 1))
            p(:, bottom_ghost(i)) = p(:, top + (i - 1))
          end do
        end if

      case default
        error stop "Unsupported location to apply the bc at in periodic_bc_t%apply_periodic_cell_gradient_bc()"
      end select
    end associate

  end subroutine apply_periodic_primitive_var_bc

  subroutine apply_periodic_reconstructed_state_bc(self, recon_rho, recon_u, recon_v, recon_p, lbounds)
    !< Apply periodic boundary conditions to the reconstructed state vector field

    class(periodic_bc_t), intent(inout) :: self

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
        if(self%do_corners) then
          do i = 1, self%n_ghost_layers
            recon_rho(:, right_ghost(i), top_ghost(i)) = recon_rho(:, left + (i - 1), bottom + (i - 1))
            recon_u(:, right_ghost(i), top_ghost(i)) = recon_u(:, left + (i - 1), bottom + (i - 1))
            recon_v(:, right_ghost(i), top_ghost(i)) = recon_v(:, left + (i - 1), bottom + (i - 1))
            recon_p(:, right_ghost(i), top_ghost(i)) = recon_p(:, left + (i - 1), bottom + (i - 1))
            recon_rho(:, right_ghost(i), bottom_ghost(i)) = recon_rho(:, left + (i - 1), top + (i - 1))
            recon_u(:, right_ghost(i), bottom_ghost(i)) = recon_u(:, left + (i - 1), top + (i - 1))
            recon_v(:, right_ghost(i), bottom_ghost(i)) = recon_v(:, left + (i - 1), top + (i - 1))
            recon_p(:, right_ghost(i), bottom_ghost(i)) = recon_p(:, left + (i - 1), top + (i - 1))
            recon_rho(:, right_ghost(i), bottom:top) = recon_rho(:, left + (i - 1), bottom:top)
            recon_u(:, right_ghost(i), bottom:top) = recon_u(:, left + (i - 1), bottom:top)
            recon_v(:, right_ghost(i), bottom:top) = recon_v(:, left + (i - 1), bottom:top)
            recon_p(:, right_ghost(i), bottom:top) = recon_p(:, left + (i - 1), bottom:top)
          end do
        else
          do i = 1, self%n_ghost_layers
            recon_rho(:, right_ghost(i), :) = recon_rho(:, left + (i - 1), :)
            recon_u(:, right_ghost(i), :) = recon_u(:, left + (i - 1), :)
            recon_v(:, right_ghost(i), :) = recon_v(:, left + (i - 1), :)
            recon_p(:, right_ghost(i), :) = recon_p(:, left + (i - 1), :)
          end do
        end if

      case('-x')
        call debug_print('Running periodic_bc_t%apply_periodic_reconstructed_state_bc() -x', __FILE__, __LINE__)
        if(self%do_corners) then
          do i = 1, self%n_ghost_layers
            recon_rho(:, left_ghost(i), top_ghost(i)) = recon_rho(:, right + (i - 1), bottom + (i - 1))
            recon_u(:, left_ghost(i), top_ghost(i)) = recon_u(:, right + (i - 1), bottom + (i - 1))
            recon_v(:, left_ghost(i), top_ghost(i)) = recon_v(:, right + (i - 1), bottom + (i - 1))
            recon_p(:, left_ghost(i), top_ghost(i)) = recon_p(:, right + (i - 1), bottom + (i - 1))
            recon_rho(:, left_ghost(i), bottom_ghost(i)) = recon_rho(:, right + (i - 1), top + (i - 1))
            recon_u(:, left_ghost(i), bottom_ghost(i)) = recon_u(:, right + (i - 1), top + (i - 1))
            recon_v(:, left_ghost(i), bottom_ghost(i)) = recon_v(:, right + (i - 1), top + (i - 1))
            recon_p(:, left_ghost(i), bottom_ghost(i)) = recon_p(:, right + (i - 1), top + (i - 1))
            recon_rho(:, left_ghost(i), bottom:top) = recon_rho(:, right + (i - 1), bottom:top)
            recon_u(:, left_ghost(i), bottom:top) = recon_u(:, right + (i - 1), bottom:top)
            recon_v(:, left_ghost(i), bottom:top) = recon_v(:, right + (i - 1), bottom:top)
            recon_p(:, left_ghost(i), bottom:top) = recon_p(:, right + (i - 1), bottom:top)
          end do
        else
          do i = 1, self%n_ghost_layers
            recon_rho(:, left_ghost(i), :) = recon_rho(:, right + (i - 1), :)
            recon_u(:, left_ghost(i), :) = recon_u(:, right + (i - 1), :)
            recon_v(:, left_ghost(i), :) = recon_v(:, right + (i - 1), :)
            recon_p(:, left_ghost(i), :) = recon_p(:, right + (i - 1), :)
          end do
        end if
      case('+y')
        call debug_print('Running periodic_bc_t%apply_periodic_reconstructed_state_bc() +y', __FILE__, __LINE__)
        if(self%do_corners) then
          do i = 1, self%n_ghost_layers
            recon_rho(:, left_ghost(i), top_ghost(i)) = recon_rho(:, right + (i - 1), bottom + (i - 1))
            recon_u(:, left_ghost(i), top_ghost(i)) = recon_u(:, right + (i - 1), bottom + (i - 1))
            recon_v(:, left_ghost(i), top_ghost(i)) = recon_v(:, right + (i - 1), bottom + (i - 1))
            recon_p(:, left_ghost(i), top_ghost(i)) = recon_p(:, right + (i - 1), bottom + (i - 1))
            recon_rho(:, right_ghost(i), top_ghost(i)) = recon_rho(:, left + (i - 1), bottom + (i - 1))
            recon_u(:, right_ghost(i), top_ghost(i)) = recon_u(:, left + (i - 1), bottom + (i - 1))
            recon_v(:, right_ghost(i), top_ghost(i)) = recon_v(:, left + (i - 1), bottom + (i - 1))
            recon_p(:, right_ghost(i), top_ghost(i)) = recon_p(:, left + (i - 1), bottom + (i - 1))
            recon_rho(:, left:right, top_ghost(i)) = recon_rho(:, left:right, bottom + (i - 1))
            recon_u(:, left:right, top_ghost(i)) = recon_u(:, left:right, bottom + (i - 1))
            recon_v(:, left:right, top_ghost(i)) = recon_v(:, left:right, bottom + (i - 1))
            recon_p(:, left:right, top_ghost(i)) = recon_p(:, left:right, bottom + (i - 1))
          end do
        else
          do i = 1, self%n_ghost_layers
            recon_rho(:, :, top_ghost(i)) = recon_rho(:, :, bottom + (i - 1))
            recon_u(:, :, top_ghost(i)) = recon_u(:, :, bottom + (i - 1))
            recon_v(:, :, top_ghost(i)) = recon_v(:, :, bottom + (i - 1))
            recon_p(:, :, top_ghost(i)) = recon_p(:, :, bottom + (i - 1))
          end do
        end if
      case('-y')
        call debug_print('Running periodic_bc_t%apply_periodic_reconstructed_state_bc() -y', __FILE__, __LINE__)
        if(self%do_corners) then
          do i = 1, self%n_ghost_layers
            recon_rho(:, left_ghost(i), bottom_ghost(i)) = recon_rho(:, right + (i - 1), top + (i - 1))
            recon_u(:, left_ghost(i), bottom_ghost(i)) = recon_u(:, right + (i - 1), top + (i - 1))
            recon_v(:, left_ghost(i), bottom_ghost(i)) = recon_v(:, right + (i - 1), top + (i - 1))
            recon_p(:, left_ghost(i), bottom_ghost(i)) = recon_p(:, right + (i - 1), top + (i - 1))
            recon_rho(:, right_ghost(i), bottom_ghost(i)) = recon_rho(:, left + (i - 1), top + (i - 1))
            recon_u(:, right_ghost(i), bottom_ghost(i)) = recon_u(:, left + (i - 1), top + (i - 1))
            recon_v(:, right_ghost(i), bottom_ghost(i)) = recon_v(:, left + (i - 1), top + (i - 1))
            recon_p(:, right_ghost(i), bottom_ghost(i)) = recon_p(:, left + (i - 1), top + (i - 1))
            recon_rho(:, left:right, bottom_ghost(i)) = recon_rho(:, left:right, top + (i - 1))
            recon_u(:, left:right, bottom_ghost(i)) = recon_u(:, left:right, top + (i - 1))
            recon_v(:, left:right, bottom_ghost(i)) = recon_v(:, left:right, top + (i - 1))
            recon_p(:, left:right, bottom_ghost(i)) = recon_p(:, left:right, top + (i - 1))
          end do
        else
          do i = 1, self%n_ghost_layers
            recon_rho(:, :, bottom_ghost(i)) = recon_rho(:, :, top + (i - 1))
            recon_u(:, :, bottom_ghost(i)) = recon_u(:, :, top + (i - 1))
            recon_v(:, :, bottom_ghost(i)) = recon_v(:, :, top + (i - 1))
            recon_p(:, :, bottom_ghost(i)) = recon_p(:, :, top + (i - 1))
          end do
        end if
      case default
        error stop "Unsupported location to apply the bc at in periodic_bc_t%apply_periodic_reconstructed_state_bc()"
      end select
    end associate

  end subroutine apply_periodic_reconstructed_state_bc

  subroutine apply_gradient_bc(self, grad_x, grad_y, lbounds)
    class(periodic_bc_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: grad_x
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: grad_y
    integer(ik)  :: i

    associate(left=>self%ilo, right=>self%ihi, bottom=>self%jlo, top=>self%jhi, &
              left_ghost=>self%ilo_ghost, right_ghost=>self%ihi_ghost, &
              bottom_ghost=>self%jlo_ghost, top_ghost=>self%jhi_ghost)

      select case(self%location)
      case('+x')
        call debug_print('Running periodic_bc_t%apply_gradient_bc() +x', __FILE__, __LINE__)
        if(self%do_corners) then
          do i = 1, self%n_ghost_layers
            grad_x(right_ghost(i), top_ghost(i)) = grad_x(left + (i - 1), bottom + (i - 1))
            grad_x(right_ghost(i), bottom_ghost(i)) = grad_x(left + (i - 1), top + (i - 1))
            grad_x(right_ghost(i), bottom:top) = grad_x(left + (i - 1), bottom:top)

            grad_y(right_ghost(i), top_ghost(i)) = grad_y(left + (i - 1), bottom + (i - 1))
            grad_y(right_ghost(i), bottom_ghost(i)) = grad_y(left + (i - 1), top + (i - 1))
            grad_y(right_ghost(i), bottom:top) = grad_y(left + (i - 1), bottom:top)
          end do
        else
          do i = 1, self%n_ghost_layers
            grad_x(right_ghost(i), :) = grad_x(left + (i - 1), :)
            grad_y(right_ghost(i), :) = grad_y(left + (i - 1), :)
          end do
        end if

      case('-x')
        call debug_print('Running periodic_bc_t%apply_gradient_bc() -x', __FILE__, __LINE__)
        if(self%do_corners) then
          do i = 1, self%n_ghost_layers
            grad_x(left_ghost(i), top_ghost(i)) = grad_x(right + (i - 1), bottom + (i - 1))
            grad_x(left_ghost(i), bottom_ghost(i)) = grad_x(right + (i - 1), top + (i - 1))
            grad_x(left_ghost(i), bottom:top) = grad_x(right + (i - 1), bottom:top)

            grad_y(left_ghost(i), top_ghost(i)) = grad_y(right + (i - 1), bottom + (i - 1))
            grad_y(left_ghost(i), bottom_ghost(i)) = grad_y(right + (i - 1), top + (i - 1))
            grad_y(left_ghost(i), bottom:top) = grad_y(right + (i - 1), bottom:top)
          end do
        else
          do i = 1, self%n_ghost_layers
            grad_x(left_ghost(i), :) = grad_x(right + (i - 1), :)
            grad_y(left_ghost(i), :) = grad_y(right + (i - 1), :)
          end do
        end if
      case('+y')
        call debug_print('Running periodic_bc_t%apply_gradient_bc() +y', __FILE__, __LINE__)
        if(self%do_corners) then
          do i = 1, self%n_ghost_layers
            grad_x(left_ghost(i), top_ghost(i)) = grad_x(right + (i - 1), bottom + (i - 1))
            grad_x(right_ghost(i), top_ghost(i)) = grad_x(left + (i - 1), bottom + (i - 1))
            grad_x(left:right, top_ghost(i)) = grad_x(left:right, bottom + (i - 1))

            grad_y(left_ghost(i), top_ghost(i)) = grad_y(right + (i - 1), bottom + (i - 1))
            grad_y(right_ghost(i), top_ghost(i)) = grad_y(left + (i - 1), bottom + (i - 1))
            grad_y(left:right, top_ghost(i)) = grad_y(left:right, bottom + (i - 1))
          end do
        else
          do i = 1, self%n_ghost_layers
            grad_x(:, top_ghost(i)) = grad_x(:, bottom + (i - 1))
            grad_y(:, top_ghost(i)) = grad_y(:, bottom + (i - 1))
          end do
        end if
      case('-y')
        call debug_print('Running periodic_bc_t%apply_gradient_bc() -y', __FILE__, __LINE__)
        if(self%do_corners) then
          do i = 1, self%n_ghost_layers
            grad_x(left_ghost(i), bottom_ghost(i)) = grad_x(right + (i - 1), top + (i - 1))
            grad_x(right_ghost(i), bottom_ghost(i)) = grad_x(left + (i - 1), top + (i - 1))
            grad_x(left:right, bottom_ghost(i)) = grad_x(left:right, top + (i - 1))

            grad_y(left_ghost(i), bottom_ghost(i)) = grad_y(right + (i - 1), top + (i - 1))
            grad_y(right_ghost(i), bottom_ghost(i)) = grad_y(left + (i - 1), top + (i - 1))
            grad_y(left:right, bottom_ghost(i)) = grad_y(left:right, top + (i - 1))
          end do
        else
          do i = 1, self%n_ghost_layers
            grad_x(:, bottom_ghost(i)) = grad_x(:, top + (i - 1))
            grad_y(:, bottom_ghost(i)) = grad_y(:, top + (i - 1))
          end do
        end if
      case default
        error stop "Unsupported location to apply the bc at in periodic_bc_t%apply_gradient_bc()"
      end select
    end associate

  end subroutine apply_gradient_bc

  subroutine finalize(self)
    type(periodic_bc_t), intent(inout) :: self
    call debug_print('Running periodic_bc_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%ilo_ghost)) deallocate(self%ilo_ghost)
    if(allocated(self%ihi_ghost)) deallocate(self%ihi_ghost)
    if(allocated(self%jlo_ghost)) deallocate(self%jlo_ghost)
    if(allocated(self%jhi_ghost)) deallocate(self%jhi_ghost)
  end subroutine

end module mod_periodic_bc
