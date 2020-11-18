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

module mod_symmetry_bc
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: enable_debug_print, debug_print
  use mod_field, only: field_2d_t
  use mod_error, only: error_msg
  use mod_grid_block_2d, only: grid_block_2d_t
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_input, only: input_t

  implicit none

  private
  public :: symmetry_bc_t, symmetry_bc_constructor

  type, extends(boundary_condition_t) :: symmetry_bc_t
  contains
    procedure, public :: apply => apply_symmetry_primitive_var_bc
    final :: finalize
  end type

contains

  function symmetry_bc_constructor(location, input, grid) result(bc)
    type(symmetry_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input
    class(grid_block_2d_t), intent(in) :: grid

    allocate(bc)
    bc%name = 'symmetry'
    bc%location = location
    call bc%set_indices(grid)
  end function symmetry_bc_constructor

  subroutine apply_symmetry_primitive_var_bc(self, rho, u, v, p)

    class(symmetry_bc_t), intent(inout) :: self
    class(field_2d_t), intent(inout) :: rho
    class(field_2d_t), intent(inout) :: u
    class(field_2d_t), intent(inout) :: v
    class(field_2d_t), intent(inout) :: p

    integer(ik) :: g

    associate(left => self%ilo, right => self%ihi, bottom => self%jlo, top => self%jhi, &
              left_ghost => self%ilo_ghost, right_ghost => self%ihi_ghost, &
              bottom_ghost => self%jlo_ghost, top_ghost => self%jhi_ghost)

      select case(self%location)
      case('+x')
        if(rho%on_ihi_bc) then
          if(enable_debug_print) call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() +x', __FILE__, __LINE__)

          ! ghost layer indexing is always innermost to outermost, e.g. g=1 is right next to the real domain
          do g = 1, self%n_ghost_layers
            rho%data(right_ghost(g), :) = rho%data(right - (g - 1), :)
            u%data(right_ghost(g), :) = -u%data(right - (g - 1), :)
            v%data(right_ghost(g), :) = v%data(right - (g - 1), :)
            p%data(right_ghost(g), :) = p%data(right - (g - 1), :)
          end do
        endif

      case('-x')
        if(rho%on_ilo_bc) then
          if(enable_debug_print) call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() -x', __FILE__, __LINE__)

          ! ghost layer indexing is always innermost to outermost, e.g. g=1 is right next to the real domain
          do g = 1, self%n_ghost_layers
            rho%data(left_ghost(g), :) = rho%data(left + (g - 1), :)
            u%data(left_ghost(g), :) = -u%data(left + (g - 1), :)
            v%data(left_ghost(g), :) = v%data(left + (g - 1), :)
            p%data(left_ghost(g), :) = p%data(left + (g - 1), :)
          end do
        end if
      case('+y')
        if(rho%on_jhi_bc) then
          if(enable_debug_print) call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() +y', __FILE__, __LINE__)

          ! ghost layer indexing is always innermost to outermost, e.g. g=1 is right next to the real domain
          do g = 1, self%n_ghost_layers
            rho%data(:, top_ghost(g)) = rho%data(:, top - (g - 1))
            u%data(:, top_ghost(g)) = u%data(:, top - (g - 1))
            v%data(:, top_ghost(g)) = -v%data(:, top - (g - 1))
            p%data(:, top_ghost(g)) = p%data(:, top - (g - 1))
          end do
        end if
      case('-y')
        if(rho%on_jlo_bc) then
          if(enable_debug_print) call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() -y', __FILE__, __LINE__)

          ! ghost layer indexing is always innermost to outermost, e.g. g=1 is right next to the real domain
          do g = 1, self%n_ghost_layers
            rho%data(:, bottom_ghost(g)) = rho%data(:, bottom + (g - 1))
            u%data(:, bottom_ghost(g)) = u%data(:, bottom + (g - 1))
            v%data(:, bottom_ghost(g)) = -v%data(:, bottom + (g - 1))
            p%data(:, bottom_ghost(g)) = p%data(:, bottom + (g - 1))
          end do
        endif
      case default
        call error_msg(module_name='mod_symmetry_bc', &
                       class_name='symmetry_bc_t', &
                       procedure_name='apply_symmetry_primitive_var_bc', &
                       message="Unsupported BC location", &
                       file_name=__FILE__, line_number=__LINE__)
      end select
    end associate

  end subroutine apply_symmetry_primitive_var_bc

  subroutine finalize(self)
    type(symmetry_bc_t), intent(inout) :: self
    if(enable_debug_print) call debug_print('Running symmetry_bc_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%ilo_ghost)) deallocate(self%ilo_ghost)
    if(allocated(self%ihi_ghost)) deallocate(self%ihi_ghost)
    if(allocated(self%jlo_ghost)) deallocate(self%jlo_ghost)
    if(allocated(self%jhi_ghost)) deallocate(self%jhi_ghost)
  end subroutine finalize
end module mod_symmetry_bc
