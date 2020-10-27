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

module mod_grid_factory

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print, set_domain_dimensionality
  use mod_input, only: input_t
  use mod_grid_block, only: grid_block_t
  use mod_grid_block_1d, only: new_1d_grid_block
  use mod_grid_block_2d, only: new_2d_grid_block
  use mod_grid_block_3d, only: new_3d_grid_block
  implicit none

  private
  public :: grid_factory

contains

  function grid_factory(input) result(grid)
    !! Factory function to create a grid object
    class(input_t), intent(in) :: input
    class(grid_block_t), pointer :: grid

    call debug_print('Creating a grid in grid_factory', __FILE__, __LINE__)

    ! select case(trim(input%grid_type))
    ! case('X')
    !   grid => new_1d_grid_block(input)
    ! case('XY')
    !   grid => new_2d_grid_block(input)
    !   ! call set_domain_dimensionality(dimensionality='2D_XY', &
    !   !                                grid_orthogonality=.true., num_ghost_layers=grid%n_halo_cells)
    ! ! case('RZ')
    ! case('XYZ')
    !   grid => new_3d_grid_block(input)
    ! case default
    !   error stop 'Unsupported grid type in grid_factory'
    ! end select
  end function grid_factory
end module mod_grid_factory
