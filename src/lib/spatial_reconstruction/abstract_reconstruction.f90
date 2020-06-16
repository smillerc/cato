module mod_abstract_reconstruction

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_regular_2d_grid, only: regular_2d_grid_t
  use mod_slope_limiter, only: slope_limiter_t
  ! use mod_flux_limiter, only: flux_limiter_t
  use mod_input, only: input_t
  use mod_grid, only: grid_t
  use mod_globals, only: debug_print

  implicit none

  private
  public :: abstract_reconstruction_t

  type, abstract :: abstract_reconstruction_t
    !< Base class for reconstruction operators

    class(grid_t), pointer :: grid => null()
    !< Pointer to the grid object, which should be managed by the finite_volume_scheme_t puppeteer class

    real(rk), dimension(:, :), pointer :: rho !< cell average density (used to find rho(P') at any (x,y))
    real(rk), dimension(:, :), pointer :: p   !< cell average pressure (used to find p(P') at any (x,y))

    real(rk), dimension(:, :), allocatable :: grad_x_rho !< x-gradient of density (used to find rho(P') at any (x,y))
    real(rk), dimension(:, :), allocatable :: grad_x_p   !< x-gradient of pressure (used to find p(P') at any (x,y))
    real(rk), dimension(:, :), allocatable :: grad_y_rho !< y-gradient of density (used to find rho(P') at any (x,y))
    real(rk), dimension(:, :), allocatable :: grad_y_p   !< y-gradient of pressure (used to find p(P') at any (x,y))

    integer(ik), public :: order = 0  !< Reconstruction order
    character(:), allocatable, public :: name  !< Name of the reconstruction scheme
    logical :: use_post_limiter = .false. !< Use the 'a posteriori' limiter (Kitamura et al.)
  contains
    procedure, public, non_overridable :: set_grid_pointer
    procedure, public, non_overridable :: set_cell_average_pointers
    procedure(initialize), public, deferred :: initialize
    procedure(reconstruct), public, deferred :: reconstruct
    procedure(reconstruct_at_point), public, deferred :: reconstruct_at_point
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

    subroutine reconstruct(self, primitive_var, reconstructed_var, lbounds, name, stage_name)
      import :: abstract_reconstruction_t
      import :: ik, rk
      class(abstract_reconstruction_t), intent(inout) :: self
      integer(ik), dimension(2), intent(in) :: lbounds
      character(len=*), intent(in) :: stage_name
      character(len=*), intent(in) :: name
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(in), contiguous :: primitive_var !< (i,j); cell primitive variable to reconstruct
      real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(out), contiguous :: reconstructed_var
      !< ((corner1:midpoint4), i, j); reconstructed variable, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint

    end subroutine reconstruct

    real(rk) function reconstruct_at_point(self, i, j, x, y, var) result(q)
      import :: abstract_reconstruction_t
      import :: ik, rk
      class(abstract_reconstruction_t), intent(in) :: self
      real(rk), intent(in) :: x, y  !< location within cell
      integer(ik), intent(in) :: i, j !< cell indices
      character(len=*), intent(in) :: var !< variable to reconstruct ('rho', or 'p')
    end function

    subroutine copy_recon(out_recon, in_recon)
      import :: abstract_reconstruction_t
      class(abstract_reconstruction_t), intent(in) :: in_recon
      class(abstract_reconstruction_t), intent(inout) :: out_recon
    end subroutine
  end interface

contains

  subroutine set_grid_pointer(self, grid)
    !< Associate the grid with data
    class(abstract_reconstruction_t), intent(inout) :: self
    class(grid_t), intent(in), target :: grid

    if(.not. associated(self%grid)) self%grid => grid
  end subroutine set_grid_pointer

  subroutine set_cell_average_pointers(self, rho, p, lbounds)
    !< Make the cell average quantities point to the real values (taken from the fluid class).
    !< These values are used later on to reconstruct at a given P'(x,y) point needed by the
    !< Mach cones and evolution operator E0
    class(abstract_reconstruction_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in), target :: rho
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in), target :: p
    self%rho => rho
    self%p => p
  end subroutine set_cell_average_pointers

end module mod_abstract_reconstruction
