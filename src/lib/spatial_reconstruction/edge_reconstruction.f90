module mod_edge_reconstruction
  !> Summary: Provide procedures to reconstruct the cell interface values
  !> Date: 05/08/2020
  !> Author: Sam Miller
  !> Notes:
  !> References:
  !>   [1] M. Berger, M. Aftosmis, S. Muman, "Analysis of Slope Limiters on Irregular Grids",
  !>       43rd AIAA Aerospace Sciences Meeting and Exhibit (2005), https://doi.org/10.2514/6.2005-490
  !>
  !>   [2] K.H. Kim, C. Kim, "Accurate, efficient and monotonic numerical methods for multi-dimensional compressible flows Part II: Multi-dimensional limiting process",
  !>       Journal of Computational Physics 208 (2005) 570–615, https://doi.org/10.1016/j.jcp.2005.02.022

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_flux_limiter, only: flux_limiter_t
  use mod_gradients

  implicit none
  private
  public :: reconstruct_edge_values

contains

  subroutine reconstruct_edge_values(q, lbounds, limiter, edge_values)
    !< Reconstruct the cell interface values, e.g. q_i-1/2, q_i+1/2. This assumes a cartesian
    !< structured square grid

    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
    !< (i,j); primitive variable to reconstruct at the edge

    real(rk), dimension(:, :, :), allocatable, intent(out) :: edge_values
    !<((bottom, right, top, left), i, j); reconstructed edge values

    type(flux_limiter_t), intent(in) :: limiter !< flux limiter used to reconstruct the edge interface

    integer(ik) :: i, j, k
    integer(ik) :: ilo, ihi, jlo, jhi
    real(rk) :: R_i !< smoothness indicator
    real(rk) :: R_j !< smoothness indicator
    real(rk) :: delta_i_minus !< q_i - q_{i-1}
    real(rk) :: delta_j_minus !< q_j - q_{j-1}
    real(rk) :: delta_i_plus  !< q_{i+1} - q_i
    real(rk) :: delta_j_plus  !< q_{j+1} - q_j

    ilo = lbound(q, dim=1)
    ihi = ubound(q, dim=1)
    jlo = lbound(q, dim=2)
    jhi = ubound(q, dim=2)

    allocate(edge_values(4, ilo:ihi, jlo:jhi))

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j, R_i, R_j) &
    !$omp private(delta_i_minus, delta_j_minus, delta_i_plus, delta_j_plus) &
    !$omp shared(q, limiter, edge_values)
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi

        delta_i_minus = q(i, i) - q(i - 1, j)
        delta_j_minus = q(i, j) - q(i, j - 1)
        delta_i_plus = q(i + 1, j) - q(i, j)
        delta_j_plus = q(i, j + 1) - q(i, j)

        ! (i-1/2, j), left
        if(abs(delta_i_plus) > 0.0_rk) then
          R_i = delta_i_minus / delta_i_plus ! this is 1/R of the normal definition
          edge_values(4, i, j) = q(i, j) - limiter%limit(R_i) * delta_i_plus
        else
          edge_values(4, i, j) = q(i, j)
        end if

        ! (i, j-1/2), bottom
        if(abs(delta_j_plus) > 0.0_rk) then
          R_j = delta_j_minus / delta_j_plus ! this is 1/R of the normal definition
          edge_values(1, i, j) = q(i, j) - limiter%limit(R_j) * delta_j_plus
        else
          edge_values(1, i, j) = q(i, j)
        end if

        ! (i+1/2, j), right
        if(abs(delta_i_minus) > 0.0_rk) then
          R_i = delta_i_plus / delta_i_minus
          edge_values(2, i, j) = q(i, j) + limiter%limit(R_i) * delta_i_minus
        else
          edge_values(2, i, j) = q(i, j)
        end if

        ! (i, j+1/2), top
        if(abs(delta_j_minus) > 0.0_rk) then
          R_j = delta_j_plus / delta_j_minus
          edge_values(3, i, j) = q(i, j) + limiter%limit(R_j) * delta_j_minus
        else
          edge_values(3, i, j) = q(i, j)
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine reconstruct_edge_values

end module mod_edge_reconstruction
