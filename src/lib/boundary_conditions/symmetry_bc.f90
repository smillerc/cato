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
  endtype

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
  endfunction symmetry_bc_constructor

  subroutine apply_symmetry_primitive_var_bc(self, rho, u, v, p)

    class(symmetry_bc_t), intent(inout) :: self
    class(field_2d_t), intent(inout) :: rho
    class(field_2d_t), intent(inout) :: u
    class(field_2d_t), intent(inout) :: v
    class(field_2d_t), intent(inout) :: p

    integer(ik) :: g

    ! left => self%ilo
    ! right => self%ihi
    ! bottom => self%jlo
    ! top => self%jhi
    ! left_ghost => self%ilo_ghost
    ! right_ghost => self%ihi_ghost
    ! bottom_ghost => self%jlo_ghost
    ! top_ghost => self%jhi_ghost

    select case(self%location)
    case('+x')
      if(rho%on_ihi_bc) then
        if(enable_debug_print) call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() +x', __FILE__, __LINE__)

        ! ghost layer indexing is always innermost to outermost, e.g. g=1 is right next to the real domain
        do g = 1, self%n_ghost_layers
          rho%data(self%ihi_ghost(g), :) = rho%data(self%ihi - (g - 1), :)
          u%data(self%ihi_ghost(g), :) = -u%data(self%ihi - (g - 1), :)
          v%data(self%ihi_ghost(g), :) = v%data(self%ihi - (g - 1), :)
          p%data(self%ihi_ghost(g), :) = p%data(self%ihi - (g - 1), :)
        enddo
      endif

    case('-x')
      if(rho%on_ilo_bc) then
        if(enable_debug_print) call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() -x', __FILE__, __LINE__)

        ! ghost layer indexing is always innermost to outermost, e.g. g=1 is right next to the real domain
        do g = 1, self%n_ghost_layers
          rho%data(self%ilo_ghost(g), :) = rho%data(self%ilo + (g - 1), :)
          u%data(self%ilo_ghost(g), :) = -u%data(self%ilo + (g - 1), :)
          v%data(self%ilo_ghost(g), :) = v%data(self%ilo + (g - 1), :)
          p%data(self%ilo_ghost(g), :) = p%data(self%ilo + (g - 1), :)
        enddo
      endif
    case('+y')
      if(rho%on_jhi_bc) then
        if(enable_debug_print) call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() +y', __FILE__, __LINE__)

        ! ghost layer indexing is always innermost to outermost, e.g. g=1 is right next to the real domain
        do g = 1, self%n_ghost_layers
          rho%data(:, self%jhi_ghost(g)) = rho%data(:, self%jhi - (g - 1))
          u%data(:, self%jhi_ghost(g)) = u%data(:, self%jhi - (g - 1))
          v%data(:, self%jhi_ghost(g)) = -v%data(:, self%jhi - (g - 1))
          p%data(:, self%jhi_ghost(g)) = p%data(:, self%jhi - (g - 1))
        enddo
      endif
    case('-y')
      if(rho%on_jlo_bc) then
        if(enable_debug_print) call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() -y', __FILE__, __LINE__)

        ! ghost layer indexing is always innermost to outermost, e.g. g=1 is right next to the real domain
        do g = 1, self%n_ghost_layers
          rho%data(:, self%jlo_ghost(g)) = rho%data(:, self%jlo + (g - 1))
          u%data(:, self%jlo_ghost(g)) = u%data(:, self%jlo + (g - 1))
          v%data(:, self%jlo_ghost(g)) = -v%data(:, self%jlo + (g - 1))
          p%data(:, self%jlo_ghost(g)) = p%data(:, self%jlo + (g - 1))
        enddo
      endif
    case default
      call error_msg(module_name='mod_symmetry_bc', &
                     class_name='symmetry_bc_t', &
                     procedure_name='apply_symmetry_primitive_var_bc', &
                     message="Unsupported BC location", &
                     file_name=__FILE__, line_number=__LINE__)
    endselect

  endsubroutine apply_symmetry_primitive_var_bc

  subroutine finalize(self)
    type(symmetry_bc_t), intent(inout) :: self
    if(enable_debug_print) call debug_print('Running symmetry_bc_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%ilo_ghost)) deallocate(self%ilo_ghost)
    if(allocated(self%ihi_ghost)) deallocate(self%ihi_ghost)
    if(allocated(self%jlo_ghost)) deallocate(self%jlo_ghost)
    if(allocated(self%jhi_ghost)) deallocate(self%jhi_ghost)
  endsubroutine finalize
endmodule mod_symmetry_bc
