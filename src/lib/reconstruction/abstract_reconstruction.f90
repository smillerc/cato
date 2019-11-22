module mod_abstract_reconstruction

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_regular_2d_grid, only: regular_2d_grid_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t
  use mod_grid, only: grid_t
  ! use slope_limiter, only: limit

  implicit none

  private
  public :: abstract_reconstruction_t

  type, abstract :: abstract_reconstruction_t
    class(grid_t), pointer :: grid
    real(rk), dimension(:, :, :), pointer :: conserved_vars

    integer(ik), public :: order  !< Reconstruction order
    character(:), allocatable, public :: name  !< Name of the reconstruction scheme
    integer(ik), dimension(2), public :: current_cell_ij  !< Current selected i,j indices
    logical, public :: cell_is_selected = .false.  !< Has the cell been selected yet
    real(rk), dimension(4), public :: cell_average = 0.0_rk !< average values of the cell conserved variables
    type(slope_limiter_t), public :: limiter  !< Slope limiter (if any)
  contains
    procedure, public, non_overridable :: set_slope_limiter
    procedure(initialize), public, deferred :: initialize
    procedure(reconstruct_point), public, deferred :: reconstruct_point
    procedure(reconstruct_domain), public, deferred :: reconstruct_domain
    ! procedure(select_and_find_gradient), public, deferred :: select_and_find_gradient
    ! procedure(finalize), private, deferred :: finalize
  end type abstract_reconstruction_t

  abstract interface
    subroutine initialize(self, input)
      import :: abstract_reconstruction_t
      import :: input_t
      class(abstract_reconstruction_t), intent(inout) :: self
      class(input_t), intent(in) :: input
    end subroutine

    pure function reconstruct_point(self, xy, cell_ij) result(U_bar)
      !< Reconstruct the value of the conserved variables (U) at location (x,y) based on the
      !> cell average and gradient (if higher order)
      import :: abstract_reconstruction_t
      import :: ik, rk

      class(abstract_reconstruction_t), intent(in) :: self
      real(rk), dimension(2), intent(in) :: xy !< (x,y) position to reconstruct
      integer(ik), dimension(2), intent(in) :: cell_ij !< cell (i,j) indices to reconstruct within
      real(rk), dimension(4) :: U_bar  !< U_bar = reconstructed [rho, u, v, p]
    end function reconstruct_point

    pure subroutine reconstruct_domain(self, reconstructed_domain)
      import :: abstract_reconstruction_t
      import :: rk

      class(abstract_reconstruction_t), intent(in) :: self
      real(rk), dimension(:, :, :, :, :), intent(inout) :: reconstructed_domain
      !< ((rho, u ,v, p), point, node/midpoint, i, j);
      !< The node/midpoint dimension just selects which set of points,
      !< e.g. 1 - all corners, 2 - all midpoints
    end subroutine reconstruct_domain

  end interface

contains
  subroutine set_slope_limiter(self, name)
    class(abstract_reconstruction_t), intent(inout) :: self
    character(len=*) :: name
    self%limiter = slope_limiter_t(name)
  end subroutine set_slope_limiter

end module mod_abstract_reconstruction
