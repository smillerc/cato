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

    ! real(rk), dimension(:, :), pointer :: rho => null() !< (i,j); pointer to primitive density data
    ! real(rk), dimension(:, :), pointer :: u => null()   !< (i,j); pointer to primitive x-velocity data
    ! real(rk), dimension(:, :), pointer :: v => null()   !< (i,j); pointer to primitive y-velocity data
    ! real(rk), dimension(:, :), pointer :: p => null()   !< (i,j); pointer to primitive pressure data
    !< Pointer to the primitive variables for each cell (rho, u, v, p)

    integer(ik), public :: order = 0  !< Reconstruction order
    character(:), allocatable, public :: name  !< Name of the reconstruction scheme

    ! real(rk), dimension(:, :, :, :), allocatable :: cell_gradient
    ! !< ((rho, u ,v, p), (d/dx, d/dy), i, j); Gradient of each cell's primitive variables

    ! real(rk), dimension(:, :, :), allocatable :: cell_average
    ! !< ((rho, u ,v, p), i, j); Cell average (based on nodes and neighbor cells)

    type(slope_limiter_t), public :: limiter  !< Slope limiter (if any)

    logical :: use_post_limiter = .false. !< Use the 'a posteriori' limiter (Kitamura et al.)

    logical :: domain_has_been_reconstructed = .false.
    logical :: cell_averages_found = .false.
  contains
    procedure, public, non_overridable :: set_slope_limiter
    procedure, public, non_overridable :: set_grid_pointer
    ! procedure, public, non_overridable :: set_primitive_vars_pointer
    ! procedure, public, non_overridable :: nullify_pointer_members
    ! procedure, public, non_overridable :: find_extrema
    ! procedure, public, non_overridable :: interpolate
    procedure(initialize), public, deferred :: initialize
    ! procedure(reconstruct_point), public, deferred :: reconstruct_point
    procedure(reconstruct), public, deferred :: reconstruct
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

    subroutine reconstruct(self, primitive_var, reconstructed_var, lbounds)
      import :: abstract_reconstruction_t
      import :: ik, rk
      class(abstract_reconstruction_t), intent(inout) :: self
      integer(ik), dimension(2), intent(in) :: lbounds
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(in) :: primitive_var !< (i,j); cell primitive variable to reconstruct
      real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(out) :: reconstructed_var
      !< ((corner1:midpoint4), i, j); reconstructed variable, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint

    end subroutine reconstruct

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

end module mod_abstract_reconstruction
