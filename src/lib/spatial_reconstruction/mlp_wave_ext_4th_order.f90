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

module mod_mlp_wave_ext_4th_order

  !< Summary: Provide class for 4th order wave-extended MLP interpolation
  !< Date: 07/09/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !<   [1]

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use math_constants, only: pi
  use mod_flux_limiter, only: flux_limiter_t, smoothness, delta
  use mod_mlp_baseline, only: mlp_baseline_t
  use mod_globals, only: n_ghost_layers, debug_print

  implicit none
  private

  real(rk), parameter :: c20 = 0.0_rk !< Coefficient C_20 in Table 1 in Ref [1]
  real(rk), parameter :: c30 = 0.0_rk !< Coefficient C_30 in Table 1 in Ref [1]
  real(rk), parameter :: c40 = 0.0_rk !< Coefficient C_40 in Table 1 in Ref [1]

  real(rk), parameter :: c21 = 0.0_rk !< Coefficient C_21 in Table 1 in Ref [1]
  real(rk), parameter :: c31 = 0.0_rk !< Coefficient C_31 in Table 1 in Ref [1]
  real(rk), parameter :: c41 = 0.0_rk !< Coefficient C_41 in Table 1 in Ref [1]

  real(rk), parameter :: c22 = 0.0_rk !< Coefficient C_22 in Table 1 in Ref [1]
  real(rk), parameter :: c32 = 0.0_rk !< Coefficient C_32 in Table 1 in Ref [1]
  real(rk), parameter :: c42 = 0.0_rk !< Coefficient C_42 in Table 1 in Ref [1]

  real(rk), parameter :: c23 = 0.0_rk !< Coefficient C_23 in Table 1 in Ref [1]
  real(rk), parameter :: c33 = 0.0_rk !< Coefficient C_33 in Table 1 in Ref [1]
  real(rk), parameter :: c43 = 0.0_rk !< Coefficient C_43 in Table 1 in Ref [1]

  real(rk), parameter :: c24 = 0.0_rk !< Coefficient C_24 in Table 1 in Ref [1]
  real(rk), parameter :: c34 = 0.0_rk !< Coefficient C_34 in Table 1 in Ref [1]
  real(rk), parameter :: c44 = 0.0_rk !< Coefficient C_44 in Table 1 in Ref [1]

  real(rk), parameter :: eta2 = 0.355_rk * pi !< Coefficient eta 2 in Table 1 in Ref [1]
  real(rk), parameter :: eta3 = 0.361_rk * pi !< Coefficient eta 2 in Table 1 in Ref [1]
  real(rk), parameter :: eta4 = 0.367_rk * pi !< Coefficient eta 2 in Table 1 in Ref [1]

end module mod_mlp_wave_ext_4th_order
