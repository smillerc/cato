module mod_abstract_reconstruction

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_regular_2d_grid, only: regular_2d_grid_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t
  use mod_grid, only: grid_t
  ! use slope_limiter, only: limit

  implicit none

  private
  public :: abstract_reconstruction_t

  type, abstract :: abstract_reconstruction_t
    !< Base class for reconstruction operators

    class(grid_t), pointer :: grid => null()
    !< Pointer to the grid object, which should be managed by the finite_volume_scheme_t puppeteer class

    real(rk), dimension(:, :, :), pointer :: conserved_vars => null()
    !< Pointer to the conserved variables for each cell (rho, rho*u, rho*v, e)

    integer(ik), public :: order = 0  !< Reconstruction order
    character(:), allocatable, public :: name  !< Name of the reconstruction scheme

    real(rk), dimension(:, :, :, :), allocatable :: cell_gradient
    !< ((d/dx, d/dy), (rho, u ,v, p), i, j); Gradient of each cell's primitive variables

    type(slope_limiter_t), public :: limiter  !< Slope limiter (if any)

    logical :: domain_has_been_reconstructed = .false.
  contains
    procedure, public, non_overridable :: set_slope_limiter
    procedure, public, non_overridable :: set_grid_pointer
    procedure, public, non_overridable :: set_conserved_vars_pointer
    procedure, public, non_overridable :: nullify_pointer_members
    procedure(initialize), public, deferred :: initialize
    procedure(reconstruct_point), public, deferred :: reconstruct_point
    procedure(reconstruct_domain), public, deferred :: reconstruct_domain
    procedure(copy_recon), public, deferred :: copy
    generic :: assignment(=) => copy
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

    pure function reconstruct_point(self, xy, cell_ij) result(V_bar)
      !< Reconstruct the value of the primitive variables (V) at location (x,y) based on the
      !> cell average and gradient (if higher order)
      import :: abstract_reconstruction_t
      import :: ik, rk

      class(abstract_reconstruction_t), intent(in) :: self
      real(rk), dimension(2), intent(in) :: xy !< (x,y) position to reconstruct
      integer(ik), dimension(2), intent(in) :: cell_ij !< cell (i,j) indices to reconstruct within
      real(rk), dimension(4) :: V_bar  !< V_bar = reconstructed primitive variables [rho, u, v, p]
    end function reconstruct_point

    pure subroutine reconstruct_domain(self, reconstructed_domain, lbounds)
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

  subroutine set_conserved_vars_pointer(self, conserved_vars, lbounds)
    !< Associate the conserved variables with data. The lbounds argument
    !< is due to the way in which the conserved vars array is indexed (due to ghost cells).
    !< This is normaly indexed starting at 0 for the i (2nd) and j (3rd) indices.
    class(abstract_reconstruction_t), intent(inout) :: self
    integer(ik), dimension(3), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), &
      intent(in), target :: conserved_vars

    if(.not. associated(self%conserved_vars)) self%conserved_vars => conserved_vars
  end subroutine set_conserved_vars_pointer

  subroutine nullify_pointer_members(self)
    class(abstract_reconstruction_t), intent(inout) :: self
    nullify(self%grid)
    nullify(self%conserved_vars)
  end subroutine nullify_pointer_members

end module mod_abstract_reconstruction
