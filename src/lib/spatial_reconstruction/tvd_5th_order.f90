module mod_tvd_5th_order
  !> Summary: Provide class for 5th order TVD edge interpolation
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

  type, extends(edge_iterpolator_t) :: tvd_5th_order_t
  contains
    procedure, public :: reconstruct_edge_values
  end type tvd_5th_order_t

contains
  subroutine reconstruct_edge_values(self, q, lbounds, limiter, edge_values)
    class(tvd_5th_order_t), intent(in) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
    !< (i,j); primitive variable to reconstruct at the edge
    type(slope_limiter_t), intent(in) :: limiter !< slope limiter used to reconstruct the edge interface
    real(rk), dimension(:, :, :), allocatable, intent(out) :: edge_values
    !<((bottom, right, top, left), i, j); reconstructed edge values
  end subroutine

end module mod_tvd_5th_order
