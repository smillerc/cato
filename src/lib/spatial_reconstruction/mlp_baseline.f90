module mod_mlp_baseline
  !< Summary: Provide base class MLP edge reconstruction
  !< Date: 06/08/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !<   [1] M. Berger, M. Aftosmis, S. Muman, "Analysis of Slope Limiters on Irregular Grids",
  !<       43rd AIAA Aerospace Sciences Meeting and Exhibit (2005), https://doi.org/10.2514/6.2005-490
  !<
  !<   [2] K.H. Kim, C. Kim, "Accurate, efficient and monotonic numerical methods for multi-dimensional compressible flows Part II: Multi-dimensional limiting process",
  !<       Journal of Computational Physics 208 (2005) 570–615, https://doi.org/10.1016/j.jcp.2005.02.022

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_flux_limiter, only: flux_limiter_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_edge_interp, only: edge_iterpolator_t

  implicit none

  private
  public :: mlp_baseline_t

  type, extends(edge_iterpolator_t), abstract :: mlp_baseline_t
  contains
    procedure, public :: get_alphas
    procedure, public :: get_rep_angles
  end type mlp_baseline_t

contains
  subroutine get_rep_angles(self, q, lbounds, theta_i, theta_j)

    !< Calculate the representative theta angles.
    class(mlp_baseline_t), intent(in) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
    !< (i,j); primitive variable to reconstruct at the edge

    real(rk), dimension(:, :), allocatable, intent(out) :: theta_i
    !< (i,j); theta_i = |(q(i,j+1) - q(i,j-1)) / (q(i+1,j) - q(i-1,j))|, i.e. 1/theta_j

    real(rk), dimension(:, :), allocatable, intent(out) :: theta_j
    !< (i,j); theta_j = |(q(i+1,j) - q(i-1,j)) / (q(i,j+1) - q(i,j-1))|, i.e. 1/theta_i

    integer(ik) :: i, j
    integer(ik) :: ilo, ihi
    integer(ik) :: jlo, jhi
    real(rk) :: numerator, denominator
    real(rk), parameter :: infinity = 1e20_rk

    allocate(theta_i(ilo:ihi, jlo:jhi))
    allocate(theta_j(ilo:ihi, jlo:jhi))

    !! $omp parallel default(none), &
    !! $omp firstprivate(ilo, ihi, jlo, jhi) &
    !! $omp private(i, j) &
    !! $omp do
    do j = jlo, jhi
      do i = ilo, ihi

        numerator = self%get_delta(q(i, j + 1), q(i, j - 1))
        denominator = self%get_delta(q(i + 1, j), q(i - 1, j))

        if(abs(numerator - denominator) < epsilon(1.0_rk) .or. & ! deltas are the same
           (abs(numerator) < tiny(1.0_rk) .and. abs(denominator) < tiny(1.0_rk))) then ! both are 0
          theta_i(i, j) = 1.0_rk
          theta_j(i, j) = 1.0_rk
        else if(abs(denominator) < tiny(1.0_rk)) then ! delta- is 0
          theta_i(i, j) = infinity
          theta_j(i, j) = 0.0_rk
        else if(abs(numerator) < tiny(1.0_rk)) then ! delta+ is 0
          theta_i(i, j) = 0.0_rk
          theta_j(i, j) = infinity
        else
          theta_i(i, j) = abs(numerator / denominator)
          theta_j(i, j) = abs(denominator / numerator)
        end if

      end do
    end do
    !! $omp end do
    !! $omp end parallel

  end subroutine get_rep_angles

  subroutine get_alphas(self)
    class(mlp_baseline_t), intent(in) :: self
  end subroutine get_alphas

end module mod_mlp_baseline
