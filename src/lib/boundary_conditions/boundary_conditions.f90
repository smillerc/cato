module mod_boundary_conditions
  !< Define the based boundary condition class

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t

  implicit none

  private
  public :: boundary_condition_t

  type, abstract :: boundary_condition_t
    character(len=32) :: name = ''
    character(len=2) :: location = '' !< Location (+x, -x, +y, or -y)
    real(rk), private :: time = 0.0_rk ! Solution time (for time dependent bc's)
    real(rk) :: max_time = 0.0_rk !< Max time in source (e.g. stop after this)
    integer(ik) :: priority = 0 !< Certain b.c.'s take priority over others in terms of when they are applied (periodic is last)
    integer(ik) :: io_unit = 0

    integer(ik) :: n_ghost_layers = 0
    integer(ik), dimension(:), allocatable :: ilo_ghost !< ilo ghost boundary cell indices
    integer(ik), dimension(:), allocatable :: ihi_ghost !< ihi ghost boundary cell indices
    integer(ik), dimension(:), allocatable :: jlo_ghost !< jlo ghost boundary cell indices
    integer(ik), dimension(:), allocatable :: jhi_ghost !< jhi ghost boundary cell indices
    integer(ik) :: ilo = 0 !< ilo real domain cell index
    integer(ik) :: ihi = 0 !< ihi real domain cell index
    integer(ik) :: jlo = 0 !< jlo real domain cell index
    integer(ik) :: jhi = 0 !< jhi real domain cell index
  contains
    procedure, public :: set_time
    procedure, public :: get_time
    procedure, public :: set_indices
    procedure(apply_primitive_var_bc), public, deferred :: apply_primitive_var_bc
    procedure(apply_gradient_bc), public, deferred :: apply_gradient_bc
    procedure(apply_reconstructed_state_bc), public, deferred :: apply_reconstructed_state_bc
  end type boundary_condition_t

  abstract interface
    subroutine apply_primitive_var_bc(self, rho, u, v, p, lbounds)
      import :: boundary_condition_t
      import :: ik, rk
      class(boundary_condition_t), intent(inout) :: self
      integer(ik), dimension(2), intent(in) :: lbounds
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: rho
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: u
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: v
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: p
    end subroutine apply_primitive_var_bc

    subroutine apply_gradient_bc(self, grad_x, grad_y, lbounds)
      import :: boundary_condition_t
      import :: ik, rk
      class(boundary_condition_t), intent(inout) :: self
      integer(ik), dimension(2), intent(in) :: lbounds
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: grad_x
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: grad_y
    end subroutine apply_gradient_bc

    subroutine apply_reconstructed_state_bc(self, recon_rho, recon_u, recon_v, recon_p, lbounds)
      import :: boundary_condition_t
      import :: ik, rk
      class(boundary_condition_t), intent(inout) :: self
      integer(ik), dimension(2), intent(in) :: lbounds
      real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_rho
      real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_u
      real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_v
      real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_p
    end subroutine apply_reconstructed_state_bc

  end interface
contains
  pure subroutine set_time(self, time)
    class(boundary_condition_t), intent(inout) :: self
    real(rk), intent(in) :: time
    self%time = time
  end subroutine

  pure real(rk) function get_time(self) result(time)
    class(boundary_condition_t), intent(in) :: self
    time = self%time
  end function

  subroutine set_indices(self, ghost_layers)
    !< Save the indices for the cells that are tagged as ghost cells. These indices will be used
    !< by the other procedures that apply the boundary conditions

    class(boundary_condition_t), intent(inout) :: self

    integer(ik), dimension(:, :), intent(in) :: ghost_layers
    !< (ilo_layers(n), ihi_layers(n), jlo_layers(n), jhi_layers(n)); indices to the ghost layers.
    !< The ilo_layers type var can be scalar or a vector, e.g. ilo_layers = [-1,0] or ilo_layers = 0

    ! Example:
    ! the cells are numbered as follows: [-1 0 | 1 2 3 4 | 5 6]
    ! the real cells are [1,2,3,4] and the ghost cells are [-1 0] and [5 6]
    self%n_ghost_layers = size(ghost_layers(1, :))

    allocate(self%ilo_ghost(self%n_ghost_layers))
    self%ilo_ghost = ghost_layers(1, :)

    allocate(self%ihi_ghost(self%n_ghost_layers))
    self%ihi_ghost = ghost_layers(2, :)

    allocate(self%jlo_ghost(self%n_ghost_layers))
    self%jlo_ghost = ghost_layers(3, :)

    allocate(self%jhi_ghost(self%n_ghost_layers))
    self%jhi_ghost = ghost_layers(4, :)

    self%ilo = maxval(self%ilo_ghost) + 1
    self%ihi = minval(self%ihi_ghost) - 1

    self%jlo = maxval(self%jlo_ghost) + 1
    self%jhi = minval(self%jhi_ghost) - 1
  end subroutine set_indices
end module mod_boundary_conditions
