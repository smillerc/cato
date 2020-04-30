module mod_abstract_reconstruction

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_regular_2d_grid, only: regular_2d_grid_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t
  use mod_grid, only: grid_t
  use mod_globals, only: debug_print, PRESSURE_FLOOR, DENSITY_FLOOR

  implicit none

  private
  public :: abstract_reconstruction_t

  type, abstract :: abstract_reconstruction_t
    !< Base class for reconstruction operators

    class(grid_t), pointer :: grid => null()
    !< Pointer to the grid object, which should be managed by the finite_volume_scheme_t puppeteer class

    real(rk), dimension(:, :), pointer :: rho => null() !< (i,j); pointer to primitive density data
    real(rk), dimension(:, :), pointer :: u => null()   !< (i,j); pointer to primitive x-velocity data
    real(rk), dimension(:, :), pointer :: v => null()   !< (i,j); pointer to primitive y-velocity data
    real(rk), dimension(:, :), pointer :: p => null()   !< (i,j); pointer to primitive pressure data
    !< Pointer to the primitive variables for each cell (rho, u, v, p)

    integer(ik), public :: order = 0  !< Reconstruction order
    character(:), allocatable, public :: name  !< Name of the reconstruction scheme

    real(rk), dimension(:, :, :, :), allocatable :: cell_gradient
    !< ((rho, u ,v, p), (d/dx, d/dy), i, j); Gradient of each cell's primitive variables

    real(rk), dimension(:, :, :), allocatable :: cell_average
    !< ((rho, u ,v, p), i, j); Cell average (based on nodes and neighbor cells)

    type(slope_limiter_t), public :: limiter  !< Slope limiter (if any)

    logical :: use_post_limiter = .false. !< Use the 'a posteriori' limiter (Kitamura et al.)

    logical :: domain_has_been_reconstructed = .false.
    logical :: cell_averages_found = .false.
  contains
    procedure, public, non_overridable :: set_slope_limiter
    procedure, public, non_overridable :: set_grid_pointer
    procedure, public, non_overridable :: set_primitive_vars_pointer
    procedure, public, non_overridable :: nullify_pointer_members
    procedure, public, non_overridable :: find_extrema
    procedure, public, non_overridable :: interpolate
    procedure(initialize), public, deferred :: initialize
    procedure(reconstruct_point), public, deferred :: reconstruct_point
    procedure(reconstruct_domain), public, deferred :: reconstruct_domain
    ! procedure(copy_recon), public, deferred :: copy
    ! generic :: assignment(=) => copy
  end type abstract_reconstruction_t

  abstract interface
    subroutine initialize(self, input, grid_target)
      import :: abstract_reconstruction_t
      import :: input_t
      import :: grid_t
      import :: ik, rk
      class(abstract_reconstruction_t), intent(inout) :: self
      class(input_t), intent(in) :: input
      class(grid_t), intent(in), target :: grid_target
    end subroutine

    function reconstruct_point(self, xy, cell_ij) result(V_bar)
      !< Reconstruct the value of the primitive variables (V) at location (x,y) based on the
      !> cell average and gradient (if higher order)
      import :: abstract_reconstruction_t
      import :: ik, rk

      class(abstract_reconstruction_t), intent(in) :: self
      real(rk), dimension(2), intent(in) :: xy !< (x,y) position to reconstruct
      integer(ik), dimension(2), intent(in) :: cell_ij !< cell (i,j) indices to reconstruct within
      real(rk), dimension(4) :: V_bar  !< V_bar = reconstructed primitive variables [rho, u, v, p]
    end function reconstruct_point

    subroutine reconstruct_domain(self, reconstructed_domain, lbounds)
      !< Reconstruct each corner/midpoint. This converts the cell centered conserved
      !< quantities [rho, rho u, rho v, e] to reconstructed primitive variables [rho, u, v, p]
      !< based on the chosen reconstruction order, e.g. using a piecewise-linear function based on the
      !< selected cell and it's neighbors.
      import :: abstract_reconstruction_t
      import :: ik, rk
      class(abstract_reconstruction_t), intent(inout) :: self
      integer(ik), dimension(5), intent(in) :: lbounds
      real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):, &
                          lbounds(4):, lbounds(5):), intent(out) :: reconstructed_domain
      !< ((rho, u ,v, p), point, node/midpoint, i, j);
      !< The node/midpoint dimension just selects which set of points,
      !< e.g. 1 - all corners, 2 - all midpoints
    end subroutine reconstruct_domain

    subroutine copy_recon(out_recon, in_recon)
      import :: abstract_reconstruction_t
      class(abstract_reconstruction_t), intent(in) :: in_recon
      class(abstract_reconstruction_t), intent(inout) :: out_recon
    end subroutine
  end interface

contains
  subroutine set_slope_limiter(self, name)
    !< Create the class's slope limiter
    class(abstract_reconstruction_t), intent(inout) :: self
    character(len=*) :: name
    self%limiter = slope_limiter_t(name)
  end subroutine set_slope_limiter

  subroutine set_grid_pointer(self, grid)
    !< Associate the grid with data
    class(abstract_reconstruction_t), intent(inout) :: self
    class(grid_t), intent(in), target :: grid

    if(.not. associated(self%grid)) self%grid => grid
  end subroutine set_grid_pointer

  subroutine set_primitive_vars_pointer(self, rho, u, v, p, lbounds)
    !< Associate the primitive variables with data. The lbounds argument
    !< is due to the way in which the conserved vars array is indexed (due to ghost cells).
    !< This is normaly indexed starting at 0 for the i (2nd) and j (3rd) indices.
    class(abstract_reconstruction_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in), target :: rho
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in), target :: u
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in), target :: v
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in), target :: p

    if(.not. associated(self%rho)) self%rho => rho
    if(.not. associated(self%u)) self%u => u
    if(.not. associated(self%v)) self%v => v
    if(.not. associated(self%p)) self%p => p
  end subroutine set_primitive_vars_pointer

  subroutine nullify_pointer_members(self)
    class(abstract_reconstruction_t), intent(inout) :: self
    if(associated(self%grid)) nullify(self%grid)
    if(associated(self%rho)) nullify(self%rho)
    if(associated(self%u)) nullify(self%u)
    if(associated(self%v)) nullify(self%v)
    if(associated(self%p)) nullify(self%p)
  end subroutine nullify_pointer_members

  subroutine find_extrema(self, i, j, U_max, U_min)
    !< At each corner and midpoint, find the min/max
    class(abstract_reconstruction_t), intent(in) :: self
    real(rk), dimension(4, 4, 2), intent(out) :: U_max !< ((rho, u, v, p), (point 1 - 4), (corner=1/midpoint=2))
    real(rk), dimension(4, 4, 2), intent(out) :: U_min !< ((rho, u, v, p), (point 1 - 4), (corner=1/midpoint=2))
    integer(ik), intent(in) :: i, j
    integer(ik) :: l
    integer(ik), parameter :: c = 1 !< corner index
    integer(ik), parameter :: m = 2 !< midpoint index

    ! call debug_print('Running abstract_reconstruction_t%find_extrema()', __FILE__, __LINE__)

    ! Find the extrema at each node point
    ! associate(U=>self%primitive_vars)

    !   ! C1
    !   do l = 1, 4
    !     U_max(l, 1, c) = max(U(l, i, j), U(l, i - 1, j), U(l, i - 1, j - 1), U(l, i, j - 1))
    !     U_min(l, 1, c) = min(U(l, i, j), U(l, i - 1, j), U(l, i - 1, j - 1), U(l, i, j - 1))
    !   end do

    !   ! C2
    !   do l = 1, 4
    !     U_max(l, 2, c) = max(U(l, i, j), U(l, i + 1, j), U(l, i + 1, j - 1), U(l, i, j - 1))
    !     U_min(l, 2, c) = min(U(l, i, j), U(l, i + 1, j), U(l, i + 1, j - 1), U(l, i, j - 1))
    !   end do

    !   ! C3
    !   do l = 1, 4
    !     U_max(l, 3, c) = max(U(l, i, j), U(l, i + 1, j), U(l, i + 1, j + 1), U(l, i, j + 1))
    !     U_min(l, 3, c) = min(U(l, i, j), U(l, i + 1, j), U(l, i + 1, j + 1), U(l, i, j + 1))
    !   end do

    !   ! C4
    !   do l = 1, 4
    !     U_max(l, 4, c) = max(U(l, i, j), U(l, i, j + 1), U(l, i - 1, j + 1), U(l, i - 1, j))
    !     U_min(l, 4, c) = min(U(l, i, j), U(l, i, j + 1), U(l, i - 1, j + 1), U(l, i - 1, j))
    !   end do

    !   ! M1
    !   do l = 1, 4
    !     U_max(l, 1, m) = max(U(l, i, j), U(l, i, j - 1))
    !     U_min(l, 1, m) = min(U(l, i, j), U(l, i, j - 1))
    !   end do

    !   ! M2
    !   do l = 1, 4
    !     U_max(l, 2, m) = max(U(l, i, j), U(l, i + 1, j))
    !     U_min(l, 2, m) = min(U(l, i, j), U(l, i + 1, j))
    !   end do

    !   ! M3
    !   do l = 1, 4
    !     U_max(l, 3, m) = max(U(l, i, j), U(l, i, j + 1))
    !     U_min(l, 3, m) = min(U(l, i, j), U(l, i, j + 1))
    !   end do

    !   ! M4
    !   do l = 1, 4
    !     U_max(l, 4, m) = max(U(l, i, j), U(l, i - 1, j))
    !     U_min(l, 4, m) = min(U(l, i, j), U(l, i - 1, j))
    !   end do

    ! end associate

  end subroutine find_extrema

  function interpolate(self, i, j, x, y, cell_gradient) result(u_tilde)
    !< Given the cell gradient and location, interpolate the value
    real(rk), dimension(4) :: u_tilde !< (rho, u, v, p); interpolated primitive variables
    class(abstract_reconstruction_t), intent(in) :: self
    integer(ik), intent(in) :: i, j !< cell indices
    real(rk), intent(in) :: x, y !< position to interpolate at
    real(rk), dimension(4, 2), intent(in), optional :: cell_gradient

    real(rk), dimension(4, 2) :: grad_u
    real(rk), dimension(2) :: centroid_xy !< (x,y) location of the cell centroid

    ! centroid_xy = self%grid%get_cell_centroid_xy(i=i, j=j)

    ! ! The gradient can be supplied in the case when we wish to use
    ! ! the limited or unlimited version
    ! if(present(cell_gradient)) then
    !   grad_u = cell_gradient ! the provided (typically limited) gradient
    ! else
    !   grad_u(:, 1) = self%cell_gradient(:, 1, i, j) ! the unlimited gradient
    !   grad_u(:, 2) = self%cell_gradient(:, 2, i, j) ! the unlimited gradient
    ! end if

    ! associate(cell_ave=>self%primitive_vars(:, i, j), &
    !           dU_dx=>grad_u(:, 1), dU_dy=>grad_u(:, 2), &
    !           x_ij=>centroid_xy(1), y_ij=>centroid_xy(2))
    !   u_tilde = cell_ave + dU_dx * (x - x_ij) + dU_dy * (y - y_ij)
    ! end associate
  end function interpolate
end module mod_abstract_reconstruction
