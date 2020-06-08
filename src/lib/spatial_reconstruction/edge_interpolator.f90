module mod_edge_interp
  !> Summary: Provide baseline class for edge interpolation schemes
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

  implicit none

  type, abstract :: edge_iterpolator_t
    private
    integer(ik), public :: order = 2 !< interpolation order
    character(len=32), public :: name = ''
  contains
    procedure, non_overridable, nopass :: get_delta
    procedure, non_overridable, nopass :: get_smoothness
    procedure(basic_interface), deferred, public :: reconstruct_edge_values
  end type edge_iterpolator_t

  abstract interface
    subroutine basic_interface(self, q, lbounds, limiter, edge_values)
      import :: edge_iterpolator_t, ik, rk, slope_limiter_t
      class(edge_iterpolator_t), intent(in) :: self
      integer(ik), dimension(2), intent(in) :: lbounds
      real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
      !< (i,j); primitive variable to reconstruct at the edge
      type(slope_limiter_t), intent(in) :: limiter !< slope limiter used to reconstruct the edge interface
      real(rk), dimension(:, :, :), allocatable, intent(out) :: edge_values
      !<((bottom, right, top, left), i, j); reconstructed edge values
    end subroutine
  end interface

contains
  real(rk) function get_delta(a, b) result(delta)
    !< Find the delta in the solution, e.g. delta = a - b. This checks for numbers
    !< near 0 and when a and b are very similar in magnitude. The aim is to avoid
    !< catastrophic cancellation and very small numbers that are essentially 0 for this scenario

    real(rk), intent(in) :: a, b
    real(rk), parameter :: rel_tol = 1e-12_rk     !< relative error tolerance
    real(rk), parameter :: abs_tol = tiny(1.0_rk) !< absolute error tolerance
    real(rk) :: abs_err !< absolute error

    delta = a - b
    abs_err = abs_tol + rel_tol * max(abs(a), abs(b))

    if(abs(delta) < epsilon(1.0_rk)) then
      delta = 0.0_rk
      return
    end if

    if(abs(a) < tiny(1.0_rk) .and. abs(b) < tiny(1.0_rk)) then
      delta = 0.0_rk
      return
    end if

    if(abs(delta) < abs_err) then
      delta = 0.0_rk
      return
    end if
  end function get_delta

  subroutine get_smoothness(minus, current, plus, R, R_inv)
    !< Calculate R = q(i+1) - q(i) / q(i) - q(i-1)
    real(rk), intent(in) :: minus   !< q(i-1)
    real(rk), intent(in) :: current !< q(i)
    real(rk), intent(in) :: plus    !< q(i+1)
    real(rk), intent(out) :: R      !< R
    real(rk), intent(out) :: R_inv  !< 1/R

    real(rk) :: delta_plus, delta_minus
    real(rk), parameter :: infinity = 1e20_rk

    delta_minus = get_delta(current, minus) ! q(i) - q(i-1)
    delta_plus = get_delta(plus, current)   ! q(i+1) - q(i)

    if(abs(delta_plus - delta_minus) < epsilon(1.0_rk) .or. & ! deltas are the same
       (abs(delta_plus) < tiny(1.0_rk) .and. abs(delta_minus) < tiny(1.0_rk)) & ! both are 0
       ) then
      R = 1.0_rk
      R_inv = 1.0_rk
    else if(abs(delta_minus) < tiny(1.0_rk)) then ! delta- is 0
      R = infinity
      R_inv = 0.0_rk
    else if(abs(delta_plus) < tiny(1.0_rk)) then ! delta+ is 0
      R = 0.0_rk
      R_inv = infinity
    else
      R = delta_plus / delta_minus
      R_inv = 1.0_rk / R
    end if
  end subroutine get_smoothness
end module mod_edge_interp
