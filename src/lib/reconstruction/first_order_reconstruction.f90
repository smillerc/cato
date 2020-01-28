#ifdef __DEBUG__
#define debug_write write
#else
#define debug_write ! write
#endif

module mod_first_order_reconstruction
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_grid, only: grid_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t

  implicit none

  private
  public :: first_order_reconstruction_t

  type, extends(abstract_reconstruction_t) :: first_order_reconstruction_t
  contains
    procedure, public :: initialize => init_first_order
    procedure, public :: reconstruct_domain
    procedure, public :: reconstruct_point
    procedure, public :: copy
    final :: finalize
  end type

contains

  subroutine init_first_order(self, input, grid_target)
    !< Construct the first_order_reconstruction_t type

    class(first_order_reconstruction_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid_target

    integer(ik) :: alloc_status

    call debug_print('Initializing second_order_reconstruction_t', __FILE__, __LINE__)

    self%order = 1
    self%name = 'cell_average_reconstruction'

    self%grid => grid_target

  end subroutine init_first_order

  subroutine finalize(self)
    !< Finalize the first_order_reconstruction_t type
    type(first_order_reconstruction_t), intent(inout) :: self
    integer(ik) :: alloc_status

    call debug_print('Running first_order_reconstruction_t%finalize()', __FILE__, __LINE__)

    if(associated(self%grid)) nullify(self%grid)
    if(associated(self%primitive_vars)) nullify(self%primitive_vars)
    if(allocated(self%cell_gradient)) then
      deallocate(self%cell_gradient, stat=alloc_status)
      if(alloc_status /= 0) then
        error stop "Unable to deallocate first_order_reconstruction_t%cell_gradient"
      end if
    end if
  end subroutine finalize

  subroutine copy(out_recon, in_recon)
    class(abstract_reconstruction_t), intent(in) :: in_recon
    class(first_order_reconstruction_t), intent(inout) :: out_recon

    call debug_print('Running first_order_reconstruction_t%copy()', __FILE__, __LINE__)

    if(associated(out_recon%grid)) nullify(out_recon%grid)
    out_recon%grid => in_recon%grid

    if(associated(out_recon%primitive_vars)) nullify(out_recon%primitive_vars)
    out_recon%primitive_vars => in_recon%primitive_vars

    if(allocated(out_recon%name)) deallocate(out_recon%name)
    allocate(out_recon%name, source=in_recon%name)

  end subroutine

  pure function reconstruct_point(self, xy, cell_ij) result(U_bar)
    !< Reconstruct the value of the primitive variables (U) at location (x,y)
    !< withing a cell (i,j)

    class(first_order_reconstruction_t), intent(in) :: self
    real(rk), dimension(2), intent(in) :: xy !< where should U_bar be reconstructed at?
    real(rk), dimension(4) :: U_bar  !< U_bar = reconstructed [rho, u, v, p]
    integer(ik), dimension(2), intent(in) :: cell_ij !< cell (i,j) indices to reconstruct within
    integer(ik) :: i, j

    i = cell_ij(1); j = cell_ij(2)
    U_bar = self%primitive_vars(:, i, j)

  end function reconstruct_point

  pure subroutine reconstruct_domain(self, reconstructed_domain, lbounds)
    !< Reconstruct the entire domain. Rather than do it a point at a time, this reuses some
    !< of the data necessary, like the cell average and gradient
    class(first_order_reconstruction_t), intent(inout) :: self
    integer(ik), dimension(5), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):, &
                        lbounds(4):, lbounds(5):), intent(out) :: reconstructed_domain

    !< ((rho, u ,v, p), point, node/midpoint, i, j);
    !< The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints

    integer(ik) :: i, j  !< cell i,j index
    integer(ik) :: n  !< node index -> i.e. is it a corner (1), or midpoint (2)
    integer(ik) :: p  !< point index -> which point in the grid element (corner 1-4, or midpoint 1-4)
    integer(ik) :: phi, nhi, ilo, ihi, jlo, jhi

    ! Bounds do not include ghost cells. Ghost cells get their
    ! reconstructed values and gradients from the boundary conditions
    phi = ubound(reconstructed_domain, dim=2)
    nhi = ubound(reconstructed_domain, dim=3)
    ilo = lbound(reconstructed_domain, dim=4) + 1
    ihi = ubound(reconstructed_domain, dim=4) - 1
    jlo = lbound(reconstructed_domain, dim=5) + 1
    jhi = ubound(reconstructed_domain, dim=5) - 1

    do concurrent(j=jlo:jhi)
      do concurrent(i=ilo:ihi)
        do concurrent(n=1:nhi)  ! First do corners, then to midpoints
          do concurrent(p=1:phi)  ! Loop through each point (N1-N4, and M1-M4)
            reconstructed_domain(:, p, n, i, j) = self%primitive_vars(:, i, j)

            if(reconstructed_domain(1, p, n, i, j) < 0) then
              error stop "Density <= 0 in first_order_reconstruction_t%reconstruct_domain"
            end if

            if(reconstructed_domain(4, p, n, i, j) < 0) then
              error stop "Pressure <= 0 in first_order_reconstruction_t%reconstruct_domain"
            end if

          end do
        end do
      end do
    end do

  end subroutine reconstruct_domain

end module mod_first_order_reconstruction
