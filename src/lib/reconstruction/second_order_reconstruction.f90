module mod_second_order_reconstruction
  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_grid, only: grid_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t

  implicit none

  private
  public :: second_order_reconstruction_t

  type, extends(abstract_reconstruction_t) :: second_order_reconstruction_t
  contains
    procedure, public :: initialize => init_second_order
    procedure, public :: reconstruct_domain
    procedure, public :: reconstruct_point
    procedure, private :: estimate_gradients
    procedure, private :: estimate_single_gradient
    final :: finalize
  end type

contains

  subroutine init_second_order(self, input, grid)
    !< Construct the second_order_reconstruction_t type
    class(second_order_reconstruction_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid

    self%order = 2
    self%name = 'piecewise_linear_reconstruction'
    call self%set_slope_limiter(name=input%slope_limiter)
    self%grid => grid
  end subroutine init_second_order

  subroutine finalize(self)
    type(second_order_reconstruction_t), intent(inout) :: self
    print *, 'Finalizing second_order_reconstruction_t'
    nullify(self%grid)
  end subroutine finalize

  pure function reconstruct_point(self, conserved_vars, xy, cell_ij) result(U_bar)
    !< Reconstruct the value of the conserved variables (U) at location (x,y) based on the
    class(second_order_reconstruction_t), intent(in) :: self
    real(rk), dimension(:, 0:, 0:), intent(in) :: conserved_vars
    real(rk), dimension(2), intent(in) :: xy !< where should U_bar be reconstructed at?
    real(rk), dimension(4) :: U_bar  !< U_bar = reconstructed [rho, u, v, p]
    integer(ik), dimension(2), intent(in) :: cell_ij !< cell (i,j) indices to reconstruct within

    real(rk), dimension(4) :: cell_ave !< cell average of (rho, u, v, p)
    real(rk), dimension(2) :: centroid_xy !< (x,y) location of the cell centroid
    real(rk), dimension(2, 4) :: grad_U !< cell gradient
    integer(ik) :: i, j
    i = cell_ij(1); j = cell_ij(2)

    cell_ave = sum(conserved_vars(:, i - 1:i + 1, j - 1:j + 1)) / 5.0_rk
    grad_U = self%estimate_gradients(conserved_vars, i, j)
    centroid_xy = self%grid%get_cell_centroid_xy(i, j)

    associate(dV_dx=>grad_U(1, :), dV_dy=>grad_U(2, :), &
              x=>xy(1), y=>xy(2), &
              x_ij=>centroid_xy(1), y_ij=>centroid_xy(2))

      U_bar = cell_ave + dV_dx * (x - x_ij) + dV_dy * (y - y_ij)
    end associate

  end function reconstruct_point

  pure subroutine reconstruct_domain(self, conserved_vars, reconstructed_domain)
    !< Reconstruct the entire domain. Rather than do it a point at a time, this reuses some
    !< of the data necessary, like the cell average and gradient
    class(second_order_reconstruction_t), intent(in) :: self
    real(rk), dimension(:, 0:, 0:), intent(in) :: conserved_vars
    real(rk), dimension(:, :, :, 0:, 0:), intent(out) :: reconstructed_domain

    !< ((rho, u ,v, p), point, node/midpoint, i, j);
    !< The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints

    real(rk), dimension(4) :: cell_ave !< cell average of (rho, u, v, p)
    real(rk), dimension(2) :: centroid_xy !< (x,y) location of the cell centroid
    real(rk), dimension(2, 4) :: grad_U !< cell gradient
    integer(ik) :: i, j  !< cell i,j index
    integer(ik) :: n  !< node index -> i.e. is it a corner (1), or midpoint (2)
    integer(ik) :: p  !< point index -> which point in the grid element (corner 1-4, or midpoint 1-4)

    do concurrent(j=self%grid%jlo_cell:self%grid%jhi_cell)
      do concurrent(i=self%grid%ilo_cell:self%grid%ihi_cell)

        ! The cell average is reused for each cell
        cell_ave = 0.25_rk*[conserved_vars(:, i - 1, j) + &
                            conserved_vars(:, i + 1, j) + &
                            conserved_vars(:, i, j + 1) + &
                            conserved_vars(:, i, j - 1)]

        grad_U = self%estimate_gradients(conserved_vars, i, j)
        centroid_xy = self%grid%get_cell_centroid_xy(i=i, j=j)

        ! First do corners, then to midpoints
        do concurrent(n=1:ubound(reconstructed_domain, dim=3))

          ! Loop through each point (N1-N4, and M1-M4)
          do concurrent(p=1:ubound(reconstructed_domain, dim=2))

            associate(U_bar=>reconstructed_domain, &
                      x=>self%grid%cell_node_xy(1, p, n, i, j), &
                      y=>self%grid%cell_node_xy(2, p, n, i, j), &
                      dV_dx=>grad_U(1, :), dV_dy=>grad_U(2, :), &
                      x_ij=>centroid_xy(1), y_ij=>centroid_xy(2))

              ! reconstructed_state(rho:p, point, node/midpoint, i, j)
              U_bar(:, p, n, i, j) = cell_ave + dV_dx * (x - x_ij) + dV_dy * (y - y_ij)

            end associate

          end do
        end do
      end do
    end do

  end subroutine reconstruct_domain

  pure function estimate_gradients(self, conserved_vars, i, j) result(gradients)
    class(second_order_reconstruction_t), intent(in) :: self
    real(rk), dimension(2, 4) :: gradients !< ([x,y], [rho,u,v,p])
    real(rk), dimension(:, 0:, 0:), intent(in) :: conserved_vars
    integer(ik), intent(in) :: i, j

    ! density
    gradients(:, 1) = self%estimate_single_gradient(conserved_vars, i, j, var_idx=1)

    ! x velocity
    gradients(:, 2) = self%estimate_single_gradient(conserved_vars, i, j, var_idx=2)

    ! y velocity
    gradients(:, 3) = self%estimate_single_gradient(conserved_vars, i, j, var_idx=3)

    ! pressure
    gradients(:, 4) = self%estimate_single_gradient(conserved_vars, i, j, var_idx=4)

  end function estimate_gradients

  pure function estimate_single_gradient(self, conserved_vars, i, j, var_idx) result(grad_v)
    !< Find the gradient of a variable (v) within a cell at indices (i,j) based on the neighbor information.
    !< See Eq. 9 in https://doi.org/10.1016/j.jcp.2006.03.018. The slope limiter is set via the constructor
    !< of this derived type.

    class(second_order_reconstruction_t), intent(in) :: self
    real(rk), dimension(:, 0:, 0:), intent(in) :: conserved_vars
    integer(ik), intent(in) :: var_idx !< index of the variable to estimate the gradient
    integer(ik), intent(in) :: i, j !< cell index
    real(rk), dimension(2) :: grad_v !< (dV/dx, dV/dy) gradient of the variable

    associate(L=>self%limiter, &
              U=>conserved_vars, v=>var_idx, &
              volume=>self%grid%get_cell_volumes(i, j), &
              n1=>self%grid%cell_edge_norm_vectors(:, 1, i, j), &
              n2=>self%grid%cell_edge_norm_vectors(:, 2, i, j), &
              n3=>self%grid%cell_edge_norm_vectors(:, 3, i, j), &
              n4=>self%grid%cell_edge_norm_vectors(:, 4, i, j), &
              delta_l1=>self%grid%cell_edge_lengths(1, i, j), &
              delta_l2=>self%grid%cell_edge_lengths(2, i, j), &
              delta_l3=>self%grid%cell_edge_lengths(3, i, j), &
              delta_l4=>self%grid%cell_edge_lengths(4, i, j))

      grad_v = (1._rk / (2.0_rk * volume)) * &
               (L%limit(U(v, i + 1, j) - U(v, i, j), U(v, i, j) - U(v, i - 1, j)) * (n2 * delta_l2 - n4 * delta_l4) + &
                L%limit(U(v, i, j + 1) - U(v, i, j), U(v, i, j) - U(v, i, j - 1)) * (n3 * delta_l3 - n1 * delta_l1))
    end associate

  end function estimate_single_gradient

end module mod_second_order_reconstruction
