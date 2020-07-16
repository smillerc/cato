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

module mod_edge_interpolator_factory
  !< Summary: Provide a factory method to create an edge interpolation scheme operator
  !< Date: 06/24/2020
  !< Author: Sam Miller
  !< Notes: The factory can be used with the following syntax:
  !<
  !<  use mod_edge_interpolator_factory, only: edge_interpolator_factory
  !<  class(edge_iterpolator_t), pointer :: edge_interpolator => null()
  !<  edge_interpolator => edge_interpolator_factory(input)
  !<  ... do stuff with the interpolator
  !<  deallocate(edge_interpolator)
  !<

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use mod_globals, only: debug_print
  use mod_input, only: input_t
  use mod_edge_interpolator, only: edge_iterpolator_t
  use mod_tvd_2nd_order, only: tvd_2nd_order_t, new_tvd_2nd_order_t
  use mod_tvd_3rd_order, only: tvd_3rd_order_t, new_tvd_3rd_order_t
  use mod_tvd_5th_order, only: tvd_5th_order_t, new_tvd_5th_order_t
  use mod_mlp_3rd_order, only: mlp_3rd_order_t
  use mod_mlp_5th_order, only: mlp_5th_order_t

  use mod_grid, only: grid_t

  implicit none

  private
  public :: edge_interpolator_factory

contains

  function edge_interpolator_factory(input, limiter_name) result(interpolator)
    class(input_t), intent(in) :: input
    class(edge_iterpolator_t), pointer :: interpolator
    character(len=*), optional :: limiter_name

    select case(trim(input%edge_interpolation_scheme))
    case('TVD2')
      interpolator => new_tvd_2nd_order_t(limiter=input%limiter)
    case('TVD3')
      interpolator => new_tvd_3rd_order_t(limiter=input%limiter)
    case('TVD5')
      interpolator => new_tvd_5th_order_t(limiter=input%limiter)
      ! case('MLP3')
      !   interpolator => new_mlp_3rd_order_t(limiter=input%limiter)
      ! case('MLP5')
      !   interpolator => new_mlp_5th_order_t(limiter=input%limiter)
    case default
      write(std_err, '(a)') "Unknown edge interpolation scheme, must be one of the following: "// &
        "'TVD2', 'TVD3', 'TVD5', 'MLP3', or 'MLP5', the input value was '"//trim(input%edge_interpolation_scheme)//"'"

      error stop "Unknown edge interpolation scheme, must be one of the following: "// &
        "'TVD2', 'TVD3', 'TVD5', 'MLP3', or 'MLP5'"
    end select

    call interpolator%initialize(limiter=input%limiter)

  end function

end module mod_edge_interpolator_factory
