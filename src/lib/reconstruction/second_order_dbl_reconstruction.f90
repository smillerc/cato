
module mod_second_order_dbl_reconstruction
  !< Summary: Reconstruct using the discontinuous bilinear recovery (DBL)
  !< References:
  !< [1] Lukacova-Medvid'ová, M., Saibertová, J., and Warnecke, G, "Finite volume evolution galerkin methods for nonlinear hyperbolic systems",
  !<     https://doi.org/10.1006/jcph.2002.7207

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_grid, only: grid_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use mod_floating_point_utils, only: equal

  implicit none

  private
  public :: second_order_dbl_reconstruction_t

  type, extends(abstract_reconstruction_t) :: second_order_dbl_reconstruction_t
    !< Implementation of a 2nd order piecewise-linear reconstruction operator. The
    !< "sgg" just means that it uses the standard Green-Gauss gradient form.
  contains
    procedure, public :: initialize
    procedure, public :: reconstruct_domain
    procedure, public :: reconstruct_point
    procedure, private, nopass :: limit
    procedure, public :: copy
    final :: finalize
  end type

contains

  subroutine initialize(self, input, grid_target)
    !< Construct the second_order_dbl_reconstruction_t type

    class(second_order_dbl_reconstruction_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid_target

    integer(ik) :: alloc_status

    call debug_print('Initializing second_order_dbl_reconstruction_t', __FILE__, __LINE__)

    self%order = 2
    self%name = 'piecewise_linear_reconstruction'

    self%grid => grid_target

    call self%set_slope_limiter(name=input%slope_limiter)

    associate(imin=>grid_target%ilo_bc_cell, imax=>grid_target%ihi_bc_cell, &
              jmin=>grid_target%jlo_bc_cell, jmax=>grid_target%jhi_bc_cell)

      allocate(self%cell_gradient(4, 2, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) then
        error stop "Unable to allocate second_order_dbl_reconstruction_t%cell_gradient"
      end if
      self%cell_gradient = 0.0_rk
    end associate

  end subroutine initialize

  subroutine finalize(self)
    !< Finalize the second_order_dbl_reconstruction_t type
    type(second_order_dbl_reconstruction_t), intent(inout) :: self
    integer(ik) :: alloc_status

    call debug_print('Running second_order_dbl_reconstruction_t%finalize()', __FILE__, __LINE__)

    if(associated(self%grid)) nullify(self%grid)
    if(associated(self%primitive_vars)) nullify(self%primitive_vars)
    if(allocated(self%cell_gradient)) then
      deallocate(self%cell_gradient, stat=alloc_status)
      if(alloc_status /= 0) then
        error stop "Unable to deallocate second_order_dbl_reconstruction_t%cell_gradient"
      end if
    end if
  end subroutine finalize

  subroutine copy(out_recon, in_recon)
    class(abstract_reconstruction_t), intent(in) :: in_recon
    class(second_order_dbl_reconstruction_t), intent(inout) :: out_recon

    call debug_print('Running second_order_dbl_reconstruction_t%copy()', __FILE__, __LINE__)

    if(associated(out_recon%grid)) nullify(out_recon%grid)
    out_recon%grid => in_recon%grid

    if(associated(out_recon%primitive_vars)) nullify(out_recon%primitive_vars)
    out_recon%primitive_vars => in_recon%primitive_vars

    if(allocated(out_recon%name)) deallocate(out_recon%name)
    allocate(out_recon%name, source=in_recon%name)

    if(allocated(out_recon%cell_gradient)) deallocate(out_recon%cell_gradient)
    allocate(out_recon%cell_gradient, source=in_recon%cell_gradient)

    out_recon%limiter = in_recon%limiter
    out_recon%domain_has_been_reconstructed = .false.
  end subroutine

  function reconstruct_point(self, xy, cell_ij) result(V_bar)
    !< Reconstruct the value of the primitive variables (U) at location (x,y)
    !< withing a cell (i,j)

    class(second_order_dbl_reconstruction_t), intent(in) :: self
    real(rk), dimension(2), intent(in) :: xy !< where should V_bar be reconstructed at?
    real(rk), dimension(4) :: V_bar  !< V_bar = reconstructed [rho, u, v, p]
    integer(ik), dimension(2), intent(in) :: cell_ij !< cell (i,j) indices to reconstruct within
    real(rk), dimension(2) :: centroid_xy !< (x,y) location of the cell centroid
    integer(ik) :: i, j

    real(rk), dimension(4) :: U_max
    real(rk), dimension(4) :: U_min
    real(rk), dimension(4) :: U_tilde
    real(rk), dimension(4) :: f_x
    real(rk), dimension(4) :: f_y
    real(rk), dimension(4) :: f_xy
    real(rk) :: h !< cell size parameter
    real(rk) :: k !< cell size parameter
    real(rk), dimension(4) :: phi_lim !< slope limiter for (rho, u, v, p)

    i = cell_ij(1); j = cell_ij(2)
    centroid_xy = self%grid%get_cell_centroid_xy(i, j)

    if(.not. self%domain_has_been_reconstructed) then
      error stop "Error in second_order_dbl_reconstruction_t%reconstruct_point(), "// &
        "domain_has_been_reconstructed is false, but should be true"
    end if

    if(i > 0 .and. j > 0 .and. i < ubound(self%primitive_vars, dim=2) .and. j < ubound(self%primitive_vars, dim=3)) then
      centroid_xy = self%grid%get_cell_centroid_xy(i=i, j=j)

      ! Note: this only works for uniform grids!!!
      associate(U=>self%primitive_vars, grid=>self%grid, v=>self%grid%cell_volume)
        h = (grid%min_dx + grid%max_dx) / 2.0_rk
        k = (grid%min_dy + grid%max_dy) / 2.0_rk
        ! Cell centered version
        f_x = (U(:, i + 1, j) - U(:, i - 1, j)) / (2.0_rk * h)
        f_y = (U(:, i, j + 1) - U(:, i, j - 1)) / (2.0_rk * k)
        f_xy = (U(:, i + 1, j + 1) - U(:, i + 1, j - 1) - U(:, i - 1, j + 1) + U(:, i - 1, j - 1)) / (4.0_rk * h * k)
      end associate

      ! call self%find_extrema(i, j, U_max, U_min)
      associate(x=>xy(1), y=>xy(2), &
                U=>self%primitive_vars, x_ij=>centroid_xy(1), y_ij=>centroid_xy(2))

        U_tilde = (x - x_ij) * f_x + (y - y_ij) * f_y + (x - x_ij) * (y - y_ij) * f_xy
        ! phi_lim = self%limit(u_ij=U(:, i, j), u_tilde=U(:, i, j) + U_tilde, u_max=U_max(:), u_min=U_min(:))
        V_bar = U(:, i, j) + U_tilde

      end associate
    else
      V_bar = self%primitive_vars(:, i, j)
    end if
  end function reconstruct_point

  subroutine reconstruct_domain(self, reconstructed_domain, lbounds)
    !< Reconstruct each corner/midpoint. This converts the cell centered conserved
    !< quantities [rho, rho u, rho v, e] to reconstructed primitive variables [rho, u, v, p]
    !< based on the chosen reconstruction order, e.g. using a piecewise-linear function based on the
    !< selected cell and it's neighbors. Rather than do it a point at a time, this reuses some
    !< of the data necessary, like the cell average and gradient
    class(second_order_dbl_reconstruction_t), intent(inout) :: self
    integer(ik), dimension(5), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):, &
                        lbounds(4):, lbounds(5):), intent(out) :: reconstructed_domain

    !< ((rho, u ,v, p), point, node/midpoint, i, j);
    !< The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints

    real(rk), dimension(2) :: centroid_xy !< (x,y) location of the cell centroid
    integer(ik) :: i, j, n, p  !< cell i,j index
    integer(ik) :: ihi, ilo, jhi, jlo
    integer(ik), parameter :: c = 1  !< corner index
    integer(ik), parameter :: m = 2  !< midpoint index
    real(rk), dimension(4, 4, 2) :: U_max
    real(rk), dimension(4, 4, 2) :: U_min
    real(rk), dimension(4) :: U_tilde
    real(rk), dimension(4) :: f_x
    real(rk), dimension(4) :: f_y
    real(rk), dimension(4) :: f_xy
    real(rk), dimension(4) :: top_right, top_left, bottom_right, bottom_left

    real(rk), dimension(2, 4, 2) :: node_pos
    real(rk) :: h !< cell size parameter
    real(rk) :: k !< cell size parameter
    real(rk), dimension(4) :: phi_lim !< slope limiter for (rho, u, v, p)

    call debug_print('Running second_order_dbl_reconstruction_t%reconstruct_domain()', __FILE__, __LINE__)
    ! Bounds do not include ghost cells. Ghost cells get their
    ! reconstructed values and gradients from the boundary conditions
    ilo = lbound(reconstructed_domain, dim=4) + 1
    ihi = ubound(reconstructed_domain, dim=4) - 1
    jlo = lbound(reconstructed_domain, dim=5) + 1
    jhi = ubound(reconstructed_domain, dim=5) - 1

    !  Reconstruction points for each cell (corners and mid-points)
    !          E3  (edge 3)
    !     C4---M3---C3
    !     |         |
    ! E4  M4   x    M2  E2
    !     |         |
    !     C1---M1---C2
    !          E1
    do j = jlo, jhi
      do i = ilo, ihi
        centroid_xy = self%grid%get_cell_centroid_xy(i=i, j=j)

        ! Note: this only works for uniform grids!!!
        associate(U=>self%primitive_vars, grid=>self%grid, v=>self%grid%cell_volume)
          h = (grid%min_dx + grid%max_dx) / 2.0_rk
          k = (grid%min_dy + grid%max_dy) / 2.0_rk
          ! Cell centered version
          f_x = (U(:, i + 1, j) - U(:, i - 1, j)) / (2.0_rk * h)
          f_y = (U(:, i, j + 1) - U(:, i, j - 1)) / (2.0_rk * k)
          f_xy = (U(:, i + 1, j + 1) - U(:, i + 1, j - 1) - U(:, i - 1, j + 1) + U(:, i - 1, j - 1)) / (4.0_rk * h * k)

        end associate

        node_pos = self%grid%cell_node_xy(:, :, :, i, j) !< (x,y), (p1,p2,p3,p4), (node/midpoint),

        call self%find_extrema(i, j, U_max, U_min)

        do n = 1, 2 ! corner / midpoint
          do p = 1, 4 ! point

            associate(x=>node_pos(1, p, n), y=>node_pos(2, p, n), &
                      U=>self%primitive_vars, x_ij=>centroid_xy(1), y_ij=>centroid_xy(2))

              U_tilde = (x - x_ij) * f_x + (y - y_ij) * f_y + (x - x_ij) * (y - y_ij) * f_xy
              phi_lim = self%limit(u_ij=U(:, i, j), u_tilde=U(:, i, j) + U_tilde, u_max=U_max(:, p, n), u_min=U_min(:, p, n))
              reconstructed_domain(:, p, n, i, j) = U(:, i, j) + phi_lim * U_tilde
            end associate
          end do
        end do
      end do
    end do

    self%domain_has_been_reconstructed = .true.
  end subroutine reconstruct_domain

  real(rk) elemental function limit(U_ij, U_tilde, U_max, U_min) result(phi_lim)
    !< Apply a limiter. Based on Eq 3.9 in Ref [1]
    real(rk), intent(in) :: U_ij
    real(rk), intent(in) :: U_tilde
    real(rk), intent(in) :: U_max
    real(rk), intent(in) :: U_min

    if(equal(U_tilde, U_ij)) then
      phi_lim = 1.0_rk
    else if(U_tilde > U_ij) then
      phi_lim = min(1.0_rk,(U_max - U_ij) / (U_tilde - U_ij))
    else if(U_tilde < U_ij) then
      phi_lim = min(1.0_rk,(U_min - U_ij) / (U_tilde - U_ij))
    end if

  end function limit

end module mod_second_order_dbl_reconstruction
