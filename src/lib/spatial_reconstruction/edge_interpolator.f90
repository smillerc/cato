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

module mod_edge_interpolator
  !< Summary: Provide baseline class for edge interpolation schemes
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

  implicit none

  type, abstract :: edge_iterpolator_t
    private
    integer(ik), public :: order = 2 !< interpolation order
    character(len=32), public :: name = ''
    character(len=32), public :: limiter_name = ''
  contains
    procedure(init), deferred, public ::  initialize
    procedure(basic_interface), deferred, public :: interpolate_edge_values
  end type edge_iterpolator_t

  abstract interface
    subroutine init(self, limiter)
      import :: edge_iterpolator_t
      class(edge_iterpolator_t), intent(inout) :: self
      character(len=*), intent(in) :: limiter
    end subroutine init

    subroutine basic_interface(self, q, lbounds, edge_values)
      import :: edge_iterpolator_t, ik, rk
      class(edge_iterpolator_t), intent(in) :: self
      integer(ik), dimension(2), intent(in) :: lbounds
      real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
      !< (i,j); primitive variable to reconstruct at the edge
      real(rk), dimension(:, :, :), allocatable, intent(out) :: edge_values
      !<((bottom, right, top, left), i, j); reconstructed edge values
    end subroutine
  end interface

contains

end module mod_edge_interpolator
