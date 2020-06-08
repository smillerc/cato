module mod_tvd_2nd_order
  !> Summary: Provide class for 2nd order TVD edge interpolation
  !> Date: 06/08/2020
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
  use mod_slope_limiter, only: slope_limiter_t
  use mod_edge_interp, only: edge_iterpolator_t

  implicit none

  type, extends(edge_iterpolator_t) :: tvd_2nd_order_t
    type(slope_limiter_t) :: limiter
  contains
    procedure, public :: initialize
    procedure, public :: reconstruct_edge_values
  end type tvd_2nd_order_t

contains
  subroutine initialize(self, limiter)
    class(tvd_2nd_order_t), intent(inout) :: self
    character(len=*), intent(in) :: limiter

    self%limiter_name = trim(limiter)
    self%limiter = slope_limiter_t(trim(limiter))

  end subroutine initialize

  subroutine reconstruct_edge_values(self, q, lbounds, edge_values)
    !< Reconstruct the cell interface values, e.g. q_i-1/2, q_i+1/2. This assumes a cartesian
    !< structured square grid

    class(tvd_2nd_order_t), intent(in) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
    !< (i,j); primitive variable to reconstruct at the edge

    real(rk), dimension(:, :, :), allocatable, intent(out) :: edge_values
    !<((bottom, right, top, left), i, j); reconstructed edge values

    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: ilo_bc, ihi_bc, jlo_bc, jhi_bc
    real(rk) :: R_i     !< R   smoothness indicator for i
    real(rk) :: R_i_inv !< 1/R smoothness indicator for i
    real(rk) :: R_j     !< R   smoothness indicator for j
    real(rk) :: R_j_inv !< 1/R smoothness indicator for j

    real(rk) :: phi_top    !< limiter for the top edge
    real(rk) :: phi_bottom !< limiter for the bottom edge
    real(rk) :: phi_left   !< limiter for the left edge
    real(rk) :: phi_right  !< limiter for the right edge

    ! call debug_print('Running reconstruct_edge_values()', __FILE__, __LINE__)

    ilo_bc = lbound(q, dim=1)
    ihi_bc = ubound(q, dim=1)
    jlo_bc = lbound(q, dim=2)
    jhi_bc = ubound(q, dim=2)

    ilo = ilo_bc ! + n_ghost_layers
    ihi = ihi_bc ! - n_ghost_layers
    jlo = jlo_bc ! + n_ghost_layers
    jhi = jhi_bc ! - n_ghost_layers

    allocate(edge_values(4, ilo:ihi, jlo:jhi))

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp private(R_i, R_i_inv, R_j, R_j_inv) &
    !$omp private(phi_bottom, phi_top, phi_left, phi_right) &
    !$omp shared(q, self, edge_values)
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi

        call self%get_smoothness(q(i - 1, j), q(i, j), q(i + 1, j), R_i, R_i_inv)
        call self%get_smoothness(q(i, j - 1), q(i, j), q(i, j + 1), R_j, R_j_inv)

        ! slope limiters
        phi_top = 0.0_rk
        phi_bottom = 0.0_rk
        phi_left = 0.0_rk
        phi_right = 0.0_rk

        ! (i, j-1/2), bottom edge
        phi_bottom = self%limiter%limit(R_j_inv)
        edge_values(1, i, j) = q(i, j) - 0.5_rk * phi_bottom * ((q(i, j + 1) - q(i, j - 1)) / 2.0_rk)

        ! (i+1/2, j), right edge
        phi_right = self%limiter%limit(R_i)
        edge_values(2, i, j) = q(i, j) + 0.5_rk * phi_right * ((q(i + 1, j) - q(i - 1, j)) / 2.0_rk)

        ! (i, j+1/2), top edge
        phi_top = self%limiter%limit(R_j)
        edge_values(3, i, j) = q(i, j) + 0.5_rk * phi_top * ((q(i, j + 1) - q(i, j - 1)) / 2.0_rk)

        ! (i-1/2, j), left edge
        phi_left = self%limiter%limit(R_i_inv)
        edge_values(4, i, j) = q(i, j) - 0.5_rk * phi_left * ((q(i + 1, j) - q(i - 1, j)) / 2.0_rk)

      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine reconstruct_edge_values
end module mod_tvd_2nd_order
