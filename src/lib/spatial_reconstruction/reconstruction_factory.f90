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

module mod_reconstruction_factory

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_error, only: error_msg
  use mod_globals, only: debug_print
  use mod_input, only: input_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_piecewise_constant_reconstruction, only: piecewise_constant_reconstruction_t
  use mod_piecewise_linear_reconstruction, only: piecewise_linear_reconstruction_t

  use mod_grid_block, only: grid_block_t

  implicit none

  private
  public :: reconstruction_factory

contains

  function reconstruction_factory(input, grid_target) result(operator)
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid_target
    class(abstract_reconstruction_t), pointer :: operator

    character(len=:), allocatable :: recon_type

    recon_type = trim(input%spatial_reconstruction)
    call debug_print('Making a "'//recon_type//'" reconstruction operator', __FILE__, __LINE__)

    select case(recon_type)
    case('piecewise_linear')
      allocate(piecewise_linear_reconstruction_t :: operator)
      call operator%initialize(input=input, grid_target=grid_target)
    case('piecewise_constant', 'cell_average')
      allocate(piecewise_constant_reconstruction_t :: operator)
      call operator%initialize(input=input, grid_target=grid_target)
    case default
      call error_msg(module_name='mod_reconstruction_factory', procedure_name='reconstruction_factory', &
                     message="Unknown reconstruction type '"//recon_type//"'", &
                     file_name=__FILE__, line_number=__LINE__)
    endselect
  endfunction

endmodule mod_reconstruction_factory
