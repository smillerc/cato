! MIT License
! Copyright (c) 2019 Sam Miller
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module mod_flux_limiter
  !< Summary: Provide various flux limiters
  !< Date: 05/07/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !<   [1] M. Berger, M. Aftosmis, S. Muman, "Analysis of Slope Limiters on Irregular Grids",
  !<       43rd AIAA Aerospace Sciences Meeting and Exhibit (2005), https://doi.org/10.2514/6.2005-490
  !<
  !<   [2] K.H. Kim, C. Kim, "Accurate, efficient and monotonic numerical methods for multi-dimensional compressible flows Part II: Multi-dimensional limiting process",
  !<       Journal of Computational Physics 208 (2005) 570–615, https://doi.org/10.1016/j.jcp.2005.02.022

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_floating_point_utils, only: equal
  use math_constants, only: pi

  implicit none

  private
  public :: flux_limiter_t, smoothness, delta

  real(rk), parameter :: EPS = 1e-30_rk !< small number epsilon, mainly used in the smoothness function
  real(rk), parameter :: PHI_EPS = 1e-5_rk !< small number epsilon, mainly used in the smoothness function

  type :: flux_limiter_t
    character(:), allocatable :: name
    procedure(limit), pointer, nopass :: limit => null()
  contains
    final :: finalize
  end type

  ! These flux limiters are to be used in the form of the following:
  ! See [1], Eq 4.
  ! To compute the edge quantities:
  !
  !     q_{i+1/2} = q_i + 0.5 * psi(R)   * (q_i - q_{i-1})
  !     q_{i-1/2} = q_i - 0.5 * psi(1/R) * (q_{i+1} - q_i)
  !
  !     The relationship between flux (psi(R)) and slope (phi(R)) limiters is the following:
  !     psi(R) = phi(R)((R+1)/2)
  !
  !     psi(R): the flux limiters given in this module
  !     q_i: primitive cell-averaged quantity
  !     q_{i+1/2}: quantity at the cell interface
  !     R = (q_{i+1} - q_i) / (q_i - q_{i-1}); difference ratio, e.g. smoothness indicator

  interface flux_limiter_t
    module procedure constructor
  end interface

  interface
    pure real(rk) function limit(R) result(psi_lim)
      import :: rk
      real(rk), intent(in) :: R
      !< u_{i+1} - u_i / u_i - u_{i-1}; ratio of forward to backward differences in the solution
    end function
  end interface

contains
  type(flux_limiter_t) function constructor(name) result(limiter)
    !< Flux limiter factory/constructor
    character(len=*), intent(in) :: name

    limiter%name = trim(name)

    select case(trim(name))
    case('none')
      limiter%limit => none
    case('minmod')
      limiter%limit => minmod
    case('van_leer')
      limiter%limit => van_leer
    case('superbee')
      limiter%limit => superbee
    case default
      error stop "Error in flux_limiter_t%constructor(): Unknown flux limiter name"
    end select
  end function

  subroutine finalize(self)
    type(flux_limiter_t), intent(inout) :: self
    if(associated(self%limit)) nullify(self%limit)
    if(allocated(self%name)) deallocate(self%name)
  end subroutine finalize

  pure real(rk) function smoothness(plus, current, minus) result(r)
    real(rk), intent(in) :: plus, current, minus
    real(rk) :: delta_plus, delta_minus
    delta_plus = plus - current
    if(abs(delta_plus) < EPS) delta_plus = EPS

    delta_minus = current - minus
    if(abs(delta_minus) < EPS) delta_minus = EPS

    r = (delta_minus + EPS) / (delta_plus + EPS)
  end function smoothness

  elemental real(rk) function delta(a, b)
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
  end function delta

  pure real(rk) function none(R) result(psi_lim)
    !< Unlimited slope
    real(rk), intent(in) :: R !< smoothness indicator

    psi_lim = 1.0_rk
  end function none

  pure real(rk) function minmod(R) result(psi_lim)
    !< Min-mod flux limiter. See Eq. 8 in [2]
    real(rk), intent(in) :: R !< smoothness indicator

    if(R < 0.0_rk .or. abs(R) < PHI_EPS) then
      psi_lim = 0.0_rk
    else if(R > 1.0_rk) then
      psi_lim = 1.0_rk
    else
      psi_lim = max(0.0_rk, min(R, 1.0_rk))
    end if
  end function minmod

  pure real(rk) function van_leer(R) result(psi_lim)
    !< van Leer flux limiter. See Eq. 9 in [2]
    real(rk), intent(in) :: R !< smoothness indicator

    if(R < 0.0_rk .or. abs(R) < PHI_EPS) then
      psi_lim = 0.0_rk
    else
      psi_lim = (R + abs(R)) / (1.0_rk + abs(R))
    end if
  end function van_leer

  pure real(rk) function superbee(R) result(psi_lim)
    !< Superbee flux limiter. See Eq. 3 in [1]
    real(rk), intent(in) :: R !< smoothness indicator

    if(R < 0.0_rk .or. abs(R) < PHI_EPS) then
      psi_lim = 0.0_rk
    else if(R > 2.0_rk) then
      psi_lim = 2.0_rk
    else
      psi_lim = max(0.0_rk, min(2.0_rk * R, 1.0_rk), min(R, 2.0_rk))
    end if
  end function superbee

end module mod_flux_limiter
