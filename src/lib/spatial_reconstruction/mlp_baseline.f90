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
  use mod_edge_interpolator, only: edge_iterpolator_t
  use mod_globals, only: n_ghost_layers

  implicit none

  private
  public :: mlp_baseline_t

  type, extends(edge_iterpolator_t), abstract :: mlp_baseline_t
  contains
    procedure, nopass, non_overridable, public :: g
    procedure, public :: get_alphas
    procedure, public :: get_rep_angles
  end type mlp_baseline_t

contains
  elemental real(rk) function g(x)
    !< Simple function, see Ref[1], Eq. 64. This is mainly just to make the code look more
    !< like the math
    real(rk), intent(in) :: x
    g = max(1.0_rk, min(2.0_rk, x))
  end function

  subroutine get_rep_angles(self, q, lbounds, tan_theta_i, tan_theta_j)
    !< Calculate the representative theta angles.
    class(mlp_baseline_t), intent(in) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
    !< (i,j); primitive variable to reconstruct at the edge

    real(rk), dimension(lbounds(1) + 1:, lbounds(2) + 1:), intent(inout) :: tan_theta_i
    !< (i,j); tan_theta_i = |(q(i,j+1) - q(i,j-1)) / (q(i+1,j) - q(i-1,j))|, i.e. 1/tan_theta_j

    real(rk), dimension(lbounds(1) + 1:, lbounds(2) + 1:), intent(inout) :: tan_theta_j
    !< (i,j); tan_theta_j = |(q(i+1,j) - q(i-1,j)) / (q(i,j+1) - q(i,j-1))|, i.e. 1/tan_theta_i

    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: ilo_bc, ihi_bc, jlo_bc, jhi_bc
    real(rk) :: numerator, denominator
    real(rk), parameter :: infinity = 1e20_rk

    ilo_bc = lbound(q, dim=1)
    ihi_bc = ubound(q, dim=1)
    jlo_bc = lbound(q, dim=2)
    jhi_bc = ubound(q, dim=2)

    ilo = ilo_bc + n_ghost_layers
    ihi = ihi_bc - n_ghost_layers
    jlo = jlo_bc + n_ghost_layers
    jhi = jhi_bc - n_ghost_layers

    !! $omp parallel default(none), &
    !! $omp firstprivate(ilo, ihi, jlo, jhi) &
    !! $omp private(i, j) &
    !! $omp do
    do j = jlo + 1, jhi - 1
      do i = ilo + 1, ihi - 1

        numerator = self%get_delta(q(i, j + 1), q(i, j - 1))
        denominator = self%get_delta(q(i + 1, j), q(i - 1, j))

        if(abs(numerator - denominator) < epsilon(1.0_rk) .or. & ! deltas are the same
           (abs(numerator) < tiny(1.0_rk) .and. abs(denominator) < tiny(1.0_rk))) then ! both are 0
          tan_theta_i(i, j) = 1.0_rk
          tan_theta_j(i, j) = 1.0_rk
        else if(abs(denominator) < tiny(1.0_rk)) then ! delta- is 0
          tan_theta_i(i, j) = infinity
          tan_theta_j(i, j) = 0.0_rk
        else if(abs(numerator) < tiny(1.0_rk)) then ! delta+ is 0
          tan_theta_i(i, j) = 0.0_rk
          tan_theta_j(i, j) = infinity
        else
          tan_theta_i(i, j) = abs(numerator / denominator)
          tan_theta_j(i, j) = abs(denominator / numerator)
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
