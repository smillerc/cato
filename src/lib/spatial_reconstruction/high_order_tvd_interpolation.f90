module mod_high_order_interpolator
  !> Summary: Provide base class for high order TVD edge interpolation
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

  private
  public :: high_order_interpolator_t

  type, extends(edge_iterpolator_t), abstract :: high_order_interpolator_t
  contains
    procedure(get_beta), deferred, nopass :: get_beta
  end type high_order_interpolator_t

  abstract interface
    pure subroutine get_beta(r, beta)
      import :: rk
      real(rk), dimension(:), intent(in) :: r
      real(rk), dimension(:), intent(out) :: beta
    end subroutine
  end interface

contains

end module mod_high_order_interpolator
