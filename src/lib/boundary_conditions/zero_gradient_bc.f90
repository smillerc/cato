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

module mod_zero_gradient_bc
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: enable_debug_print, debug_print
  use mod_field, only: field_2d_t
  use mod_grid_block_2d, only: grid_block_2d_t
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_input, only: input_t

  implicit none

  private
  public :: zero_gradient_bc_t, zero_gradient_bc_constructor

  type, extends(boundary_condition_t) :: zero_gradient_bc_t
  contains
    procedure, public :: apply => apply_zero_gradient_primitive_var_bc
    final :: finalize
  endtype
contains

  function zero_gradient_bc_constructor(location, input, grid) result(bc)
    type(zero_gradient_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input
    class(grid_block_2d_t), intent(in) :: grid

    allocate(bc)
    bc%name = 'zero_gradient'
    bc%location = location
    bc%priority = 1
    call bc%set_indices(grid)
  endfunction zero_gradient_bc_constructor

  subroutine apply_zero_gradient_primitive_var_bc(self, rho, u, v, p)

    class(zero_gradient_bc_t), intent(inout) :: self
    class(field_2d_t), intent(inout) :: rho
    class(field_2d_t), intent(inout) :: u
    class(field_2d_t), intent(inout) :: v
    class(field_2d_t), intent(inout) :: p

    integer(ik) :: i

    select case(self%location)
    case('+x')
      if(rho%on_ihi_bc) then
        if(enable_debug_print) then
          call debug_print('Running zero_gradient_bc_t%apply_zero_gradient_primitive_var_bc() +x', &
                           __FILE__, __LINE__)
        endif
        do i = 1, self%n_ghost_layers
          rho%data(self%ihi_ghost(i), :) = rho%data(self%ihi, :)
          u%data(self%ihi_ghost(i), :) = u%data(self%ihi, :)
          v%data(self%ihi_ghost(i), :) = v%data(self%ihi, :)
          p%data(self%ihi_ghost(i), :) = p%data(self%ihi, :)
        enddo
      endif
    case('-x')
      if(rho%on_ilo_bc) then
        if(enable_debug_print) then
          call debug_print('Running zero_gradient_bc_t%apply_zero_gradient_primitive_var_bc() -x', &
                           __FILE__, __LINE__)
        endif
        do i = 1, self%n_ghost_layers
          rho%data(self%ilo_ghost(i), :) = rho%data(self%ilo, :)
          u%data(self%ilo_ghost(i), :) = u%data(self%ilo, :)
          v%data(self%ilo_ghost(i), :) = v%data(self%ilo, :)
          p%data(self%ilo_ghost(i), :) = p%data(self%ilo, :)
        enddo
      endif
    case('+y')
      if(rho%on_jhi_bc) then
        if(enable_debug_print) then
          call debug_print('Running zero_gradient_bc_t%apply_zero_gradient_primitive_var_bc() +y', &
                           __FILE__, __LINE__)
        endif
        do i = 1, self%n_ghost_layers
          rho%data(:, self%jhi_ghost(i)) = rho%data(:, self%jhi)
          u%data(:, self%jhi_ghost(i)) = u%data(:, self%jhi)
          v%data(:, self%jhi_ghost(i)) = v%data(:, self%jhi)
          p%data(:, self%jhi_ghost(i)) = p%data(:, self%jhi)
        enddo
      endif
    case('-y')
      if(rho%on_jlo_bc) then
        if(enable_debug_print) then
          call debug_print('Running zero_gradient_bc_t%apply_zero_gradient_primitive_var_bc() -y', &
                           __FILE__, __LINE__)
        endif
        do i = 1, self%n_ghost_layers
          rho%data(:, self%jlo_ghost(i)) = rho%data(:, self%jlo)
          u%data(:, self%jlo_ghost(i)) = u%data(:, self%jlo)
          v%data(:, self%jlo_ghost(i)) = v%data(:, self%jlo)
          p%data(:, self%jlo_ghost(i)) = p%data(:, self%jlo)
        enddo
      endif
    case default
      error stop "Unsupported location to apply the bc at in zero_gradient_bc_t%apply_zero_gradient_cell_gradient_bc()"
    endselect

  endsubroutine apply_zero_gradient_primitive_var_bc

  subroutine finalize(self)
    type(zero_gradient_bc_t), intent(inout) :: self
    if(enable_debug_print) call debug_print('Running zero_gradient_bc_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%ilo_ghost)) deallocate(self%ilo_ghost)
    if(allocated(self%ihi_ghost)) deallocate(self%ihi_ghost)
    if(allocated(self%jlo_ghost)) deallocate(self%jlo_ghost)
    if(allocated(self%jhi_ghost)) deallocate(self%jhi_ghost)
  endsubroutine finalize

endmodule mod_zero_gradient_bc
