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

module mod_gnoffo_aux_limiter

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use math_constants, only: pi
  implicit none

  private
  public :: gnoffo_aux_limiter
contains

  real(rk) elemental function gnoffo_aux_limiter(p_unlim, p_i, p_j) result(phi_g)
    !< The Gnoffo Auxiliary limiter as described in [1]
    !< References
    !< [1] Keiichi Kitamura, Atsushi Hashimoto, Simple a posteriori slope limiter (Post Limiter) for high resolution and efficient flow computations
    !<     https://doi.org/10.1016/j.jcp.2017.04.002
    !< [2] P.A. Gnoffo, Updates to Multi-Dimensional Flux Reconstruction for Hypersonic Simulations on Tetrahedral Grids
    !<     https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20100003598.pdf

    real(rk), intent(in) :: p_unlim
    real(rk), intent(in) :: p_i
    real(rk), intent(in) :: p_j
    real(rk) :: phi
    real(rk) :: psi
    real(rk), parameter :: psi_max = 3.0_rk
    real(rk), parameter :: psi_min = 2.0_rk

    real(rk) :: p_max, p_min

    ! Eq. 5e
    p_max = max(p_unlim, p_i, p_j)

    ! Eq. 5d
    p_min = min(p_unlim, p_i, p_j)

    ! Eq. 5c
    psi = p_max / p_min

    ! Eq. 5b
    phi = min(1.0_rk, max(0.0_rk,(psi_max - psi) / (psi_max - psi_min)))

    ! Eq. 5a
    phi_g = (1.0_rk - cos(phi * pi)) / 2.0_rk

  endfunction

endmodule mod_gnoffo_aux_limiter
