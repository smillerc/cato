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
  use mod_error, only: error_msg
  use mod_globals, only: debug_print
  use mod_input, only: input_t
  use mod_edge_interpolator, only: edge_iterpolator_t
  use mod_tvd_2nd_order, only: tvd_2nd_order_t, new_tvd_2nd_order_t
  use mod_tvd_3rd_order, only: tvd_3rd_order_t, new_tvd_3rd_order_t
  use mod_tvd_5th_order, only: tvd_5th_order_t, new_tvd_5th_order_t
  use mod_grid, only: grid_t

  implicit none

  private
  public :: edge_interpolator_factory

contains

  function edge_interpolator_factory(input, limiter_name) result(interpolator)
    class(input_t), intent(in) :: input
    class(edge_iterpolator_t), pointer :: interpolator
    character(len=*), optional :: limiter_name !< option to override the input limiter type
    character(len=:), allocatable :: limiter

    ! Some of the flux solvers use more than one limiter at the same time, so
    ! this is an optional argument. M-AUSMPW+ used both superbee and the user's limiter
    ! of choice and interpolated between the two versions, like superbee + minmod.
    if(present(limiter_name)) then
      limiter = trim(limiter_name)
    else
      limiter = input%limiter
    end if

    select case(trim(input%spatial_reconstruction))
    case('TVD2')
      interpolator => new_tvd_2nd_order_t(limiter=limiter)
    case('TVD3')
      interpolator => new_tvd_3rd_order_t(limiter=limiter)
    case('TVD5')
      interpolator => new_tvd_5th_order_t(limiter=limiter)
    case default
      call error_msg(module='mod_edge_interpolator_factory', procedure='edge_interpolator_factory', &
          message="Unknown edge interpolation scheme '"//trim(input%spatial_reconstruction)//"', must be one of the following: "// &
                     "'TVD2', 'TVD3', 'TVD5', 'MLP3', or 'MLP5'", &
                     file_name=__FILE__, line_number=__LINE__)
    end select

    call interpolator%initialize(limiter=limiter)
    deallocate(limiter)

  end function edge_interpolator_factory

end module mod_edge_interpolator_factory
