module mod_mlp_3rd_order
  !> Summary: Provide class for 3rd order MLP edge interpolation
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
  use mod_mlp_baseline, only: mlp_baseline_t

  implicit none

  type, extends(mlp_baseline_t) :: mlp_3rd_order_t
  contains
    procedure, public :: reconstruct_edge_values
    procedure, public, nopass :: get_beta
  end type mlp_3rd_order_t

contains

  function mlp_3rd_order_constructor() result(mlp_3rd_order)
    type(mlp_3rd_order_t) :: mlp_3rd_order

    ! allocate(mlp_3rd_order)
    mlp_3rd_order%order = 3
    mlp_3rd_order%name = 'MLP3'

  end function mlp_3rd_order_constructor

  pure subroutine get_beta(r, beta)
    real(rk), dimension(:), intent(in) :: r
    real(rk), dimension(:), intent(out) :: beta
    integer(ik) :: i, ilo, ihi

    ilo = lbound(r, dim=1)
    ihi = ubound(r, dim=1)
    do i = ilo, ihi
      beta(i) = (1.0_rk + 2.0_rk * r(i)) / 3.0_rk
    end do
  end subroutine

  subroutine reconstruct_edge_values(self, q, lbounds, edge_values)
    class(mlp_3rd_order_t), intent(in) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
    !< (i,j); primitive variable to reconstruct at the edge
    real(rk), dimension(:, :, :), allocatable, intent(out) :: edge_values
    !<((bottom, right, top, left), i, j); reconstructed edge values
  end subroutine

end module mod_mlp_3rd_order
