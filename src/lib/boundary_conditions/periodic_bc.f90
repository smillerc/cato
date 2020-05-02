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
    procedure, public :: copy => copy_periodic_bc

  end type periodic_bc_t

contains

  function periodic_bc_constructor(location, input) result(bc)
    type(periodic_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input

    allocate(bc)
    bc%name = 'periodic'
    bc%location = location

    if(trim(input%plus_x_bc) == 'periodic' .and. &
       trim(input%minus_x_bc) == 'periodic' .and. &
       trim(input%plus_y_bc) == 'periodic' .and. &
       trim(input%minus_y_bc) == 'periodic') then

      bc%do_corners = .true.
    end if

  end function periodic_bc_constructor

  subroutine copy_periodic_bc(out_bc, in_bc)
    class(boundary_condition_t), intent(in) :: in_bc
    class(periodic_bc_t), intent(inout) :: out_bc

    call debug_print('Running periodic_bc_t%copy_periodic_bc()', __FILE__, __LINE__)
    ! if (allocated(out_bc%name)) deallocate(out_bc%name)
    ! allocate(out_bc%name, source=in_bc%name)
    out_bc%name = in_bc%name
    out_bc%location = in_bc%location
  end subroutine copy_periodic_bc

  subroutine apply_periodic_primitive_var_bc(self, rho, u, v, p, lbounds)

    class(periodic_bc_t), intent(inout) :: self
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
      call debug_print('Running periodic_bc_t%apply_periodic_primitive_var_bc() +x', __FILE__, __LINE__)
      if(self%do_corners) then
        rho(right_ghost, top_ghost) = rho(left, bottom)
        u(right_ghost, top_ghost) = u(left, bottom)
        v(right_ghost, top_ghost) = v(left, bottom)
        p(right_ghost, top_ghost) = p(left, bottom)
        rho(right_ghost, bottom_ghost) = rho(left, top)
        u(right_ghost, bottom_ghost) = u(left, top)
        v(right_ghost, bottom_ghost) = v(left, top)
        p(right_ghost, bottom_ghost) = p(left, top)
        rho(right_ghost, bottom:top) = rho(left, bottom:top)
        u(right_ghost, bottom:top) = u(left, bottom:top)
        v(right_ghost, bottom:top) = v(left, bottom:top)
        p(right_ghost, bottom:top) = p(left, bottom:top)

      else
        rho(right_ghost, :) = rho(left, :)
        u(right_ghost, :) = u(left, :)
        v(right_ghost, :) = v(left, :)
        p(right_ghost, :) = p(left, :)
      end if

    case('-x')
      call debug_print('Running periodic_bc_t%apply_periodic_primitive_var_bc() -x', __FILE__, __LINE__)
      if(self%do_corners) then
        rho(left_ghost, top_ghost) = rho(right, bottom)
        u(left_ghost, top_ghost) = u(right, bottom)
        v(left_ghost, top_ghost) = v(right, bottom)
        p(left_ghost, top_ghost) = p(right, bottom)
        rho(left_ghost, bottom_ghost) = rho(right, top)
        u(left_ghost, bottom_ghost) = u(right, top)
        v(left_ghost, bottom_ghost) = v(right, top)
        p(left_ghost, bottom_ghost) = p(right, top)
        rho(left_ghost, bottom:top) = rho(right, bottom:top)
        u(left_ghost, bottom:top) = u(right, bottom:top)
        v(left_ghost, bottom:top) = v(right, bottom:top)
        p(left_ghost, bottom:top) = p(right, bottom:top)
      else
        rho(left_ghost, :) = rho(right, :)
        u(left_ghost, :) = u(right, :)
        v(left_ghost, :) = v(right, :)
        p(left_ghost, :) = p(right, :)
      end if

    case('+y')
      call debug_print('Running periodic_bc_t%apply_periodic_primitive_var_bc() +y', __FILE__, __LINE__)
      if(self%do_corners) then
        rho(left_ghost, top_ghost) = rho(right, bottom)
        u(left_ghost, top_ghost) = u(right, bottom)
        v(left_ghost, top_ghost) = v(right, bottom)
        p(left_ghost, top_ghost) = p(right, bottom)
        rho(right_ghost, top_ghost) = rho(left, bottom)
        u(right_ghost, top_ghost) = u(left, bottom)
        v(right_ghost, top_ghost) = v(left, bottom)
        p(right_ghost, top_ghost) = p(left, bottom)
        rho(left:right, top_ghost) = rho(left:right, bottom)
        u(left:right, top_ghost) = u(left:right, bottom)
        v(left:right, top_ghost) = v(left:right, bottom)
        p(left:right, top_ghost) = p(left:right, bottom)
      else
        rho(:, top_ghost) = rho(:, bottom)
        u(:, top_ghost) = u(:, bottom)
        v(:, top_ghost) = v(:, bottom)
        p(:, top_ghost) = p(:, bottom)
      end if

    case('-y')
      call debug_print('Running periodic_bc_t%apply_periodic_primitive_var_bc() -y', __FILE__, __LINE__)
      if(self%do_corners) then
        rho(left_ghost, bottom_ghost) = rho(right, top)
        u(left_ghost, bottom_ghost) = u(right, top)
        v(left_ghost, bottom_ghost) = v(right, top)
        p(left_ghost, bottom_ghost) = p(right, top)
        rho(right_ghost, bottom_ghost) = rho(left, top)
        u(right_ghost, bottom_ghost) = u(left, top)
        v(right_ghost, bottom_ghost) = v(left, top)
        p(right_ghost, bottom_ghost) = p(left, top)
        rho(left:right, bottom_ghost) = rho(left:right, top)
        u(left:right, bottom_ghost) = u(left:right, top)
        v(left:right, bottom_ghost) = v(left:right, top)
        p(left:right, bottom_ghost) = p(left:right, top)
      else
        rho(:, bottom_ghost) = rho(:, top)
        u(:, bottom_ghost) = u(:, top)
        v(:, bottom_ghost) = v(:, top)
        p(:, bottom_ghost) = p(:, top)
      end if

    case default
      error stop "Unsupported location to apply the bc at in periodic_bc_t%apply_periodic_cell_gradient_bc()"
    end select

  end subroutine apply_periodic_primitive_var_bc

  subroutine apply_periodic_reconstructed_state_bc(self, recon_rho, recon_u, recon_v, recon_p, lbounds)
    !< Apply periodic boundary conditions to the reconstructed state vector field

    class(periodic_bc_t), intent(inout) :: self

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
      call debug_print('Running periodic_bc_t%apply_periodic_reconstructed_state_bc() +x', __FILE__, __LINE__)
      if(self%do_corners) then
        recon_rho(:, right_ghost, top_ghost) = recon_rho(:, left, bottom)
        recon_u(:, right_ghost, top_ghost) = recon_u(:, left, bottom)
        recon_v(:, right_ghost, top_ghost) = recon_v(:, left, bottom)
        recon_p(:, right_ghost, top_ghost) = recon_p(:, left, bottom)
        recon_rho(:, right_ghost, bottom_ghost) = recon_rho(:, left, top)
        recon_u(:, right_ghost, bottom_ghost) = recon_u(:, left, top)
        recon_v(:, right_ghost, bottom_ghost) = recon_v(:, left, top)
        recon_p(:, right_ghost, bottom_ghost) = recon_p(:, left, top)
        recon_rho(:, right_ghost, bottom:top) = recon_rho(:, left, bottom:top)
        recon_u(:, right_ghost, bottom:top) = recon_u(:, left, bottom:top)
        recon_v(:, right_ghost, bottom:top) = recon_v(:, left, bottom:top)
        recon_p(:, right_ghost, bottom:top) = recon_p(:, left, bottom:top)
      else
        recon_rho(:, right_ghost, :) = recon_rho(:, left, :)
        recon_u(:, right_ghost, :) = recon_u(:, left, :)
        recon_v(:, right_ghost, :) = recon_v(:, left, :)
        recon_p(:, right_ghost, :) = recon_p(:, left, :)
      end if

    case('-x')
      call debug_print('Running periodic_bc_t%apply_periodic_reconstructed_state_bc() -x', __FILE__, __LINE__)
      if(self%do_corners) then
        recon_rho(:, left_ghost, top_ghost) = recon_rho(:, right, bottom)
        recon_u(:, left_ghost, top_ghost) = recon_u(:, right, bottom)
        recon_v(:, left_ghost, top_ghost) = recon_v(:, right, bottom)
        recon_p(:, left_ghost, top_ghost) = recon_p(:, right, bottom)
        recon_rho(:, left_ghost, bottom_ghost) = recon_rho(:, right, top)
        recon_u(:, left_ghost, bottom_ghost) = recon_u(:, right, top)
        recon_v(:, left_ghost, bottom_ghost) = recon_v(:, right, top)
        recon_p(:, left_ghost, bottom_ghost) = recon_p(:, right, top)
        recon_rho(:, left_ghost, bottom:top) = recon_rho(:, right, bottom:top)
        recon_u(:, left_ghost, bottom:top) = recon_u(:, right, bottom:top)
        recon_v(:, left_ghost, bottom:top) = recon_v(:, right, bottom:top)
        recon_p(:, left_ghost, bottom:top) = recon_p(:, right, bottom:top)
      else
        recon_rho(:, left_ghost, :) = recon_rho(:, right, :)
        recon_u(:, left_ghost, :) = recon_u(:, right, :)
        recon_v(:, left_ghost, :) = recon_v(:, right, :)
        recon_p(:, left_ghost, :) = recon_p(:, right, :)
      end if
    case('+y')
      call debug_print('Running periodic_bc_t%apply_periodic_reconstructed_state_bc() +y', __FILE__, __LINE__)
      if(self%do_corners) then
        recon_rho(:, left_ghost, top_ghost) = recon_rho(:, right, bottom)
        recon_u(:, left_ghost, top_ghost) = recon_u(:, right, bottom)
        recon_v(:, left_ghost, top_ghost) = recon_v(:, right, bottom)
        recon_p(:, left_ghost, top_ghost) = recon_p(:, right, bottom)
        recon_rho(:, right_ghost, top_ghost) = recon_rho(:, left, bottom)
        recon_u(:, right_ghost, top_ghost) = recon_u(:, left, bottom)
        recon_v(:, right_ghost, top_ghost) = recon_v(:, left, bottom)
        recon_p(:, right_ghost, top_ghost) = recon_p(:, left, bottom)
        recon_rho(:, left:right, top_ghost) = recon_rho(:, left:right, bottom)
        recon_u(:, left:right, top_ghost) = recon_u(:, left:right, bottom)
        recon_v(:, left:right, top_ghost) = recon_v(:, left:right, bottom)
        recon_p(:, left:right, top_ghost) = recon_p(:, left:right, bottom)
      else
        recon_rho(:, :, top_ghost) = recon_rho(:, :, bottom)
        recon_u(:, :, top_ghost) = recon_u(:, :, bottom)
        recon_v(:, :, top_ghost) = recon_v(:, :, bottom)
        recon_p(:, :, top_ghost) = recon_p(:, :, bottom)
      end if
    case('-y')
      call debug_print('Running periodic_bc_t%apply_periodic_reconstructed_state_bc() -y', __FILE__, __LINE__)
      if(self%do_corners) then
        recon_rho(:, left_ghost, bottom_ghost) = recon_rho(:, right, top)
        recon_u(:, left_ghost, bottom_ghost) = recon_u(:, right, top)
        recon_v(:, left_ghost, bottom_ghost) = recon_v(:, right, top)
        recon_p(:, left_ghost, bottom_ghost) = recon_p(:, right, top)
        recon_rho(:, right_ghost, bottom_ghost) = recon_rho(:, left, top)
        recon_u(:, right_ghost, bottom_ghost) = recon_u(:, left, top)
        recon_v(:, right_ghost, bottom_ghost) = recon_v(:, left, top)
        recon_p(:, right_ghost, bottom_ghost) = recon_p(:, left, top)
        recon_rho(:, left:right, bottom_ghost) = recon_rho(:, left:right, top)
        recon_u(:, left:right, bottom_ghost) = recon_u(:, left:right, top)
        recon_v(:, left:right, bottom_ghost) = recon_v(:, left:right, top)
        recon_p(:, left:right, bottom_ghost) = recon_p(:, left:right, top)
      else
        recon_rho(:, :, bottom_ghost) = recon_rho(:, :, top)
        recon_u(:, :, bottom_ghost) = recon_u(:, :, top)
        recon_v(:, :, bottom_ghost) = recon_v(:, :, top)
        recon_p(:, :, bottom_ghost) = recon_p(:, :, top)
      end if
    case default
      error stop "Unsupported location to apply the bc at in periodic_bc_t%apply_periodic_reconstructed_state_bc()"
    end select

  end subroutine apply_periodic_reconstructed_state_bc

end module mod_periodic_bc
