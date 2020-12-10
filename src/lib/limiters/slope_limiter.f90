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

module mod_slope_limiter
  !< Summary: Define the various slope limiters
  !< Date: 05/07/2020
  !< Author: Sam Miller
  !< Notes: The slope limiters in the form listed below is based on [1] the slope (not flux)
  !<       limiters. Slope and flux limiters are very similar, but have some subtle differences.
  !<       See the reference for more info. The formulas are based on Eq 10.
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
  public :: slope_limiter_t

  type :: slope_limiter_t
    character(:), allocatable :: name
    procedure(limit), pointer, nopass :: limit => null()
  contains
    final :: finalize
  endtype

  ! These slope limiters are to be used in the form of the following:
  ! See [1], Eq 5.
  ! To compute the edge quantities:
  !
  !     q_{i+1/2} = q_i + 0.5 * phi(R)   * (q_{i+1} - q_{i-1}) / 2
  !     q_{i-1/2} = q_i - 0.5 * phi(1/R) * (q_{i+1} - q_{i-1}) / 2
  !
  !     The relationship between flux (psi(R)) and slope (phi(R)) limiters is the following:
  !     psi(R) = phi(R)((R+1)/2)
  !
  !     phi(R): the slope limiters given in this module
  !     q_i: primitive cell-averaged quantity
  !     q_{i+1/2}: quantity at the cell interface
  !     R = (q_{i+1} - q_i) / (q_i - q_{i-1}); difference ratio, e.g. smoothness indicator
  !
  ! A symmetric form, using "f" instead of R can also be used, where
  ! f is defined as f = (q_{i} - q_{i-1}) / (q_{i+1} - q_{i})

  interface slope_limiter_t
    module procedure constructor
  endinterface

  interface
    real(rk) function limit(R) result(phi_lim)
      import :: rk
      real(rk), intent(in) :: R
    endfunction
  endinterface

contains

  type(slope_limiter_t) function constructor(name, use_symmetric_form) result(limiter)
    !< Slope limiter factory/constructor
    character(len=*), intent(in) :: name
    logical, optional, intent(in) :: use_symmetric_form !< use the f form?
    logical :: use_sym_form

    use_sym_form = .false.

    if(present(use_symmetric_form)) then
      use_sym_form = use_symmetric_form
    endif

    limiter%name = trim(name)

    if(use_sym_form) then
      select case(trim(name))
      case('none')
        write(*, '(a)') "Selecting 'none' as the slope limiter"
        limiter%limit => none
      case('barth_jesperson')
        write(*, '(a)') "Selecting 'barth_jesperson_sym_form' as the slope limiter"
        limiter%limit => barth_jespersen_sym_form
      case('minmod')
        write(*, '(a)') "Selecting 'minmod_sym_form' as the slope limiter"
        limiter%limit => minmod_sym_form
      case('van_leer')
        write(*, '(a)') "Selecting 'van_leer_sym_form' as the slope limiter"
        limiter%limit => van_leer_sym_form
      case default
        error stop "Error in slope_limiter_t%constructor(): Unknown slope limiter name"
      endselect
    else
      select case(trim(name))
      case('none')
        write(*, '(a)') "Selecting 'none' as the slope limiter"
        limiter%limit => none
      case('barth_jesperson')
        write(*, '(a)') "Selecting 'barth_jespersen' as the slope limiter"
        limiter%limit => barth_jespersen
      case('minmod')
        write(*, '(a)') "Selecting 'minmod' as the slope limiter"
        limiter%limit => minmod
      case('van_leer')
        write(*, '(a)') "Selecting 'van_leer' as the slope limiter"
        limiter%limit => van_leer
      case default
        error stop "Error in slope_limiter_t%constructor(): Unknown symmetric slope limiter name"
      endselect
    endif

  endfunction

  subroutine finalize(self)
    type(slope_limiter_t), intent(inout) :: self
    if(associated(self%limit)) nullify(self%limit)
    if(allocated(self%name)) deallocate(self%name)
  endsubroutine finalize

  pure real(rk) function none(R) result(phi_lim)
    !< Unlimited slope
    real(rk), intent(in) :: R !< R  = (u_{i+1} - u_i) / (u_i - u_{i-1})
    phi_lim = 1.0_rk
  endfunction none

  pure real(rk) function barth_jespersen(R) result(phi_lim)
    !< Barth Jesperson slope limiter. See Eq. 8 in [1]
    real(rk), intent(in) :: R !< smoothness indicator

    if(R < 0.0_rk) then
      phi_lim = 0.0_rk
    else
      phi_lim = min(1.0_rk, 4.0_rk / (R + 1.0_rk),(4.0_rk * R) / (R + 1.0_rk))
    endif
  endfunction barth_jespersen

  pure real(rk) function van_leer(R) result(phi_lim)
    !< van Leer slope limiter. See Eq. 8 in [1]
    real(rk), intent(in) :: R !< smoothness indicator

    if(R < 0.0_rk) then
      phi_lim = 0.0_rk
    else
      phi_lim = 4.0_rk * R / (R + 1.0_rk)**2
    endif
  endfunction van_leer

  pure real(rk) function minmod(R) result(phi_lim)
    !< Min-mod slope limiter. See Eq. 8 in [1]
    real(rk), intent(in) :: R !< smoothness indicator

    if(R < 0.0_rk) then
      phi_lim = 0.0_rk
    else
      phi_lim = min(2.0_rk / (1.0_rk + R),(2.0_rk * R) / (1.0_rk + R))
    endif
  endfunction minmod

  pure real(rk) function barth_jespersen_sym_form(f) result(phi_lim)
    !< Barth Jesperson slope limiter. See Eq. 8 in [1]
    real(rk), intent(in) :: f !< symmetric smoothness indicator

    phi_lim = min(1.0_rk, 4.0_rk * f, 4.0_rk * (1.0_rk - f))
  endfunction barth_jespersen_sym_form

  pure real(rk) function van_leer_sym_form(f) result(phi_lim)
    !< van Leer slope limiter. See Eq. 8 in [1]
    real(rk), intent(in) :: f !< symmetric smoothness indicator

    phi_lim = 4.0_rk * f * (1.0_rk - f)
  endfunction van_leer_sym_form

  pure real(rk) function minmod_sym_form(f) result(phi_lim)
    !< Min-mod slope limiter. See Eq. 8 in [1]
    real(rk), intent(in) :: f !< symmetric smoothness indicator

    phi_lim = min(2.0_rk * f, 2.0_rk * (1.0_rk - f))
  endfunction minmod_sym_form
endmodule mod_slope_limiter
