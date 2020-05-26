module mod_gradients
  !> Summary: Provide procedures to calculate the gradient of a particular variable
  !> Date: 05/08/2020
  !> Author: Sam Miller
  !> Notes:
  !> References:
  !      [1]

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_floating_point_utils, only: neumaier_sum
  use mod_grid, only: grid_t
  use mod_globals, only: n_ghost_layers

  implicit none
  logical, parameter :: filter_small_gradients = .true.
  real(rk), parameter :: SMALL = 1e-12_rk
  real(rk), parameter :: SMALL_GRAD = 1e-4_rk

  private
  public :: green_gauss_gradient

contains

  subroutine green_gauss_gradient(edge_vars, lbounds, grid, grad_x, grad_y)
    !< Estimate the slope-limited gradient of the primitive variables in the cell (i,j). This assumes
    !< a quadrilateral structured grid
    class(grid_t), intent(in) :: grid
    integer(ik), dimension(3), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: edge_vars
    !< ((edge 1:n), i, j); primitive variable interpolated on the cell interface/edge

    real(rk), dimension(:, :), allocatable, intent(out) :: grad_x !< (i, j); gradient in the x-direction
    real(rk), dimension(:, :), allocatable, intent(out) :: grad_y !< (i, j); gradient in the y-direction

    integer(ik) :: i, j, k
    integer(ik) :: ilo, ihi, jlo, jhi

    real(rk), dimension(4) :: edge_lengths  !< length of each face
    real(rk), dimension(4) :: n_x  !< normal vectors of each face
    real(rk), dimension(4) :: n_y  !< normal vectors of each face
    real(rk) :: d_dx, d_dy
    real(rk), dimension(4) :: edge_val_x, edge_val_y
    real(rk) :: esum_naive_x
    real(rk) :: esum_naive_y
    real(rk) :: esum_neumaier_x
    real(rk) :: esum_neumaier_y
    real(rk) :: max_edge_len, min_edge_len, diff

    ilo = lbound(edge_vars, dim=2)
    ihi = ubound(edge_vars, dim=2)
    jlo = lbound(edge_vars, dim=3)
    jhi = ubound(edge_vars, dim=3)

    allocate(grad_x(ilo:ihi, jlo:jhi))
    allocate(grad_y(ilo:ihi, jlo:jhi))

    ! print*, 'finding gradients'
    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp private(n_x, n_y, edge_lengths, diff) &
    !$omp shared(grad_x, grad_y, edge_vars, grid) &
    !$omp reduction(+:d_dx), reduction(+:d_dy) &
    !$omp reduction(max:max_edge_len), reduction(min:min_edge_len)
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi

        ! Edge (face) interface data
        edge_lengths(1) = grid%cell_edge_lengths(1, i, j - 1)  ! bottom
        edge_lengths(2) = grid%cell_edge_lengths(2, i + 1, j)  ! right
        edge_lengths(3) = grid%cell_edge_lengths(3, i, j + 1)  ! top
        edge_lengths(4) = grid%cell_edge_lengths(4, i - 1, j)  ! left

        max_edge_len = maxval(edge_lengths)
        min_edge_len = minval(edge_lengths)
        diff = max_edge_len - min_edge_len
        if(diff < 2.0_rk * epsilon(1.0_rk)) then
          edge_lengths = max_edge_len
        end if

        n_x(1) = grid%cell_edge_norm_vectors(1, 1, i, j - 1)  ! bottom
        n_y(1) = grid%cell_edge_norm_vectors(2, 1, i, j - 1)  ! bottom
        n_x(2) = grid%cell_edge_norm_vectors(1, 2, i + 1, j)  ! right
        n_y(2) = grid%cell_edge_norm_vectors(2, 2, i + 1, j)  ! right
        n_x(3) = grid%cell_edge_norm_vectors(1, 3, i, j + 1)  ! top
        n_y(3) = grid%cell_edge_norm_vectors(2, 3, i, j + 1)  ! top
        n_x(4) = grid%cell_edge_norm_vectors(1, 4, i - 1, j)  ! left
        n_y(4) = grid%cell_edge_norm_vectors(2, 4, i - 1, j)  ! left

        d_dx = sum(edge_vars(:, i, j) * n_x * edge_lengths)
        d_dy = sum(edge_vars(:, i, j) * n_y * edge_lengths)

        if(filter_small_gradients) then
          if(abs(d_dx) < epsilon(1.0_rk)) then
            grad_x(i, j) = 0.0_rk
          else
            grad_x(i, j) = d_dx / grid%cell_volume(i, j)
          end if

          if(abs(d_dy) < epsilon(1.0_rk)) then
            grad_y(i, j) = 0.0_rk
          else
            grad_y(i, j) = d_dy / grid%cell_volume(i, j)
          end if
        else
          grad_x(i, j) = d_dx / grid%cell_volume(i, j)
          grad_y(i, j) = d_dy / grid%cell_volume(i, j)
        end if

      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine green_gauss_gradient

end module mod_gradients
