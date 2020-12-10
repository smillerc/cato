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

module math_constants
  !< Summary: Provide global math constants
  !< Author: Sam Miller

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  real(rk), parameter :: pi = 4.0_rk * atan(1.0_rk)
  real(rk), parameter :: universal_gas_const = 8.31446261815324_rk

contains
  real(rk) elemental function rad2deg(rad) result(deg)
    real(rk), intent(in) :: rad
    deg = rad * 180.0_rk / pi
  endfunction

  real(rk) elemental function deg2rad(deg) result(rad)
    real(rk), intent(in) :: deg
    rad = deg * pi / 180.0_rk
  endfunction
endmodule math_constants
