module mod_periodic_bc

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_input, only: input_t

  implicit none

  private
  public :: periodic_bc_t

  type, extends(boundary_condition_t) :: periodic_bc_t
  contains
    procedure, public :: initialize => init_periodic_bc
    procedure, public :: apply_bc => apply_periodic_bc
    procedure, private, nopass :: apply_plus_x
    procedure, private, nopass :: apply_plus_y
    procedure, private, nopass :: apply_minus_x
    procedure, private, nopass :: apply_minus_y

  end type periodic_bc_t

contains

  subroutine init_periodic_bc(self, location)
    class(periodic_bc_t), intent(inout) :: self
    ! class(input_t), intent(in) :: input
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)

    self%name = 'periodic'
    self%location = location
  end subroutine init_periodic_bc

  subroutine apply_periodic_bc(self, conserved_vars, reconstructed_state)
    class(periodic_bc_t), intent(in) :: self
    real(rk), dimension(:, 0:, 0:), intent(inout) :: conserved_vars
    real(rk), dimension(:, :, :, 0:, 0:), intent(inout) :: reconstructed_state

    select case(self%location)
    case('+x')
      call self%apply_plus_x(conserved_vars, reconstructed_state)
    case('-x')
      call self%apply_minus_x(conserved_vars, reconstructed_state)
    case('+y')
      call self%apply_plus_y(conserved_vars, reconstructed_state)
    case('-y')
      call self%apply_minus_y(conserved_vars, reconstructed_state)
    case default
      error stop "Unsupported location to apply the bc at in periodic_bc_t%apply_periodic_bc()"
    end select

  end subroutine apply_periodic_bc

  pure subroutine apply_plus_x(conserved_vars, reconstructed_state)
    ! class(periodic_bc_t), intent(in) :: self
    real(rk), dimension(:, 0:, 0:), intent(inout) :: conserved_vars
    real(rk), dimension(:, :, :, 0:, 0:), intent(inout) :: reconstructed_state

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

    ! top right corner (ghost cell) = bottom left (real cell)
    conserved_vars(:, right_ghost, top_ghost) = conserved_vars(:, left, bottom)
    reconstructed_state(:, :, :, right_ghost, top_ghost) = reconstructed_state(:, :, :, left, bottom)

    ! bottom right corner (ghost cell) = top left (real cell)
    conserved_vars(:, right_ghost, bottom_ghost) = conserved_vars(:, left, top)
    reconstructed_state(:, :, :, right_ghost, bottom_ghost) = reconstructed_state(:, :, :, left, top)

    ! right middle section (ghost cells) = left middle (real cells)
    conserved_vars(:, right_ghost, bottom:top) = conserved_vars(:, left, bottom:top)
    reconstructed_state(:, :, :, right_ghost, bottom:top) = reconstructed_state(:, :, :, left, bottom:top)

  end subroutine

  pure subroutine apply_plus_y(conserved_vars, reconstructed_state)
    ! class(periodic_bc_t), intent(in) :: self
    real(rk), dimension(:, 0:, 0:), intent(inout) :: conserved_vars
    real(rk), dimension(:, :, :, 0:, 0:), intent(inout) :: reconstructed_state

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

    ! top left corner (ghost cell) = bottom right (real cell)
    conserved_vars(:, left_ghost, top_ghost) = conserved_vars(:, right, bottom)
    reconstructed_state(:, :, :, left_ghost, top_ghost) = reconstructed_state(:, :, :, right, bottom)

    ! top right corner (ghost cell) = bottom left (real cell)
    conserved_vars(:, right_ghost, top_ghost) = conserved_vars(:, left, bottom)
    reconstructed_state(:, :, :, right_ghost, top_ghost) = reconstructed_state(:, :, :, left, bottom)

    ! top middle section (ghost cells) = bottom middle section (real cells)
    conserved_vars(:, left:right, top_ghost) = conserved_vars(:, left:right, bottom)
    reconstructed_state(:, :, :, left:right, top_ghost) = reconstructed_state(:, :, :, left:right, bottom)

  end subroutine

  pure subroutine apply_minus_x(conserved_vars, reconstructed_state)
    ! class(periodic_bc_t), intent(in) :: self
    real(rk), dimension(:, 0:, 0:), intent(inout) :: conserved_vars
    real(rk), dimension(:, :, :, 0:, 0:), intent(inout) :: reconstructed_state

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
    reconstructed_state(:, :, :, left_ghost, top_ghost) = reconstructed_state(:, :, :, right, bottom)

    conserved_vars(:, left_ghost, bottom_ghost) = conserved_vars(:, right, top)
    reconstructed_state(:, :, :, left_ghost, bottom_ghost) = reconstructed_state(:, :, :, right, top)

    conserved_vars(:, left_ghost, bottom:top) = conserved_vars(:, right, bottom:top)
    reconstructed_state(:, :, :, left_ghost, bottom:top) = reconstructed_state(:, :, :, right, bottom:top)
  end subroutine

  pure subroutine apply_minus_y(conserved_vars, reconstructed_state)
    ! class(periodic_bc_t), intent(in) :: self
    real(rk), dimension(:, 0:, 0:), intent(inout) :: conserved_vars
    real(rk), dimension(:, :, :, 0:, 0:), intent(inout) :: reconstructed_state

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
    reconstructed_state(:, :, :, left_ghost, bottom_ghost) = reconstructed_state(:, :, :, right, top)

    conserved_vars(:, right_ghost, bottom_ghost) = conserved_vars(:, left, top)
    reconstructed_state(:, :, :, right_ghost, bottom_ghost) = reconstructed_state(:, :, :, left, top)

    conserved_vars(:, left:right, bottom_ghost) = conserved_vars(:, left:right, top)
    reconstructed_state(:, :, :, left:right, bottom_ghost) = reconstructed_state(:, :, :, left:right, top)
  end subroutine

end module mod_periodic_bc
