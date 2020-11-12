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

module mod_bc_factory
  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t
  use mod_error, only: error_msg
  use mod_grid_block_2d, only: grid_block_2d_t
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_periodic_bc, only: periodic_bc_t, periodic_bc_constructor
  use mod_symmetry_bc, only: symmetry_bc_t, symmetry_bc_constructor
  use mod_pressure_input_bc, only: pressure_input_bc_t, pressure_input_bc_constructor
  use mod_zero_gradient_bc, only: zero_gradient_bc_t, zero_gradient_bc_constructor
  ! use mod_vacuum_bc, only: vacuum_bc_t, vacuum_bc_constructor
  implicit none

  private
  public :: bc_factory
contains

  function bc_factory(bc_type, location, input, grid, time) result(bc)
    !< Factory function to create a boundary condition object
    character(len=*), intent(in) :: bc_type
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(boundary_condition_t), pointer :: bc
    class(input_t), intent(in) :: input
    class(grid_block_2d_t), intent(in) :: grid
    real(rk), intent(in) :: time

    select case(trim(bc_type))
    case('periodic')
      bc => periodic_bc_constructor(location, input, grid)
      bc%priority = 0
    case('symmetry')
      bc => symmetry_bc_constructor(location, input, grid)
      bc%priority = 1
    case('pressure_input')
      bc => pressure_input_bc_constructor(location, input, grid, time)
      bc%priority = 2
    case('zero_gradient')
      bc => zero_gradient_bc_constructor(location, input, grid)
      bc%priority = 1
    case default
      call error_msg(module_name='mod_bc_factory', &
                     procedure_name='bc_factory', &
                     message="Unsupported boundary condition type in bc_factory: '" // trim(bc_type) // "'", &
                     file_name=__FILE__, line_number=__LINE__)
    end select

  end function bc_factory

end module mod_bc_factory
