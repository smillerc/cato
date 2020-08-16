! MIT License
! Copyright (c) 2020 Sam Miller
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

module mod_muscl_interpolator_factory
  !< Summary: Provide a factory method to create an edge interpolation scheme based on MUSCL
  !< Date: 08/04/2020
  !< Author: Sam Miller
  !< Notes: The factory can be used with the following syntax:
  !<
  !<  use mod_muscl_interpolator_factory, only: muscl_interpolator_factory
  !<  class(muscl_interpolation_t), pointer :: edge_interpolator => null()
  !<  edge_interpolator => muscl_interpolator_factory(input)
  !<  ... do stuff with the interpolator
  !<  deallocate(edge_interpolator)
  !<

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use mod_globals, only: debug_print
  use mod_error, only: error_msg
  use mod_input, only: input_t
  use mod_muscl_interpolation, only: muscl_interpolation_t
  use mod_muscl_tvd2, only: muscl_tvd2_t, new_muscl_tvd2
  use mod_muscl_tvd3, only: muscl_tvd3_t, new_muscl_tvd3
  use mod_muscl_tvd5, only: muscl_tvd5_t, new_muscl_tvd5
  use mod_muscl_mlp, only: muscl_mlp_t, new_muscl_mlp
  use mod_muscl_e_mlp, only: muscl_e_mlp_t, new_muscl_e_mlp

  use mod_grid, only: grid_t

  implicit none

  private
  public :: muscl_interpolator_factory

contains

  function muscl_interpolator_factory(input, limiter_name) result(interpolator)
    class(input_t), intent(in) :: input
    class(muscl_interpolation_t), pointer :: interpolator
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

    select case(trim(input%limiter))
    case('superbee', 'van_leer', 'minmod')
      interpolator => new_muscl_tvd2(limiter=limiter)
    case('TVD3')
      interpolator => new_muscl_tvd3(limiter=limiter)
    case('TVD5')
      interpolator => new_muscl_tvd5(limiter=limiter)
    case('MLP3')
      interpolator => new_muscl_mlp(limiter=limiter, order=3)
    case('MLP5')
      interpolator => new_muscl_mlp(limiter=limiter, order=5)
    case('e-MLP3')
      interpolator => new_muscl_e_mlp(limiter=limiter, order=3)
    case('e-MLP5')
      interpolator => new_muscl_e_mlp(limiter=limiter, order=5)
    case default
      call error_msg(module='mod_muscl_interpolator_factory', procedure='muscl_interpolator_factory', &
                     message="Unknown edge interpolation scheme, must be one of the following: "// &
                     "['TVD2', '(e-)MLP3', or '(e-)MLP5]', the input value was '"// &
                     trim(input%spatial_reconstruction)//"'", &
                     file_name=__FILE__, line_number=__LINE__)
    end select

    deallocate(limiter)

  end function

end module mod_muscl_interpolator_factory
