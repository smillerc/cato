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
  end type
contains

  function zero_gradient_bc_constructor(location, input, grid) result(bc)
    type(zero_gradient_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input
    class(grid_block_2d_t), intent(in) :: grid

    allocate(bc)
    bc%name = 'zero_gradient'
    bc%location = location
    call bc%set_indices(grid)
  end function zero_gradient_bc_constructor

  subroutine apply_zero_gradient_primitive_var_bc(self, rho, u, v, p)

    class(zero_gradient_bc_t), intent(inout) :: self
    class(field_2d_t), intent(inout) :: rho
    class(field_2d_t), intent(inout) :: u
    class(field_2d_t), intent(inout) :: v
    class(field_2d_t), intent(inout) :: p

    integer(ik) :: i

    ! associate(ilo => rho%lbounds(1), ihi => rho%ubounds(1), &
    !           jlo => rho%lbounds(2), jhi => rho%ubounds(2), &
    !           ilo_halo => rho%lbounds_halo(1), ihi_halo => rho%ubounds_halo(1), &
    !           jlo_halo => rho%lbounds_halo(2), jhi_halo => rho%ubounds_halo(2), &
    !           nh => rho%n_halo_cells, &
    !           neighbors => rho%neighbors)

    associate(left => self%ilo, right => self%ihi, bottom => self%jlo, top => self%jhi, &
              left_ghost => self%ilo_ghost, right_ghost => self%ihi_ghost, &
              bottom_ghost => self%jlo_ghost, top_ghost => self%jhi_ghost)

      select case(self%location)
      case('+x')
        if(rho%on_ihi_bc) then
          if(enable_debug_print) then
            call debug_print('Running zero_gradient_bc_t%apply_zero_gradient_primitive_var_bc() +x', &
                            __FILE__, __LINE__)
          end if
          do i = 1, self%n_ghost_layers
            rho%data(right_ghost(i), :) = rho%data(right, :)
            u%data(right_ghost(i), :) = u%data(right, :)
            v%data(right_ghost(i), :) = v%data(right, :)
            p%data(right_ghost(i), :) = p%data(right, :)
          end do
        end if
      case('-x')
        if(rho%on_ilo_bc) then
          if(enable_debug_print) then
            call debug_print('Running zero_gradient_bc_t%apply_zero_gradient_primitive_var_bc() -x', &
                            __FILE__, __LINE__)
          end if
          do i = 1, self%n_ghost_layers
            rho%data(left_ghost(i), :) = rho%data(left, :)
            u%data(left_ghost(i), :) = u%data(left, :)
            v%data(left_ghost(i), :) = v%data(left, :)
            p%data(left_ghost(i), :) = p%data(left, :)
          end do
        end if
      case('+y')
        if(rho%on_jhi_bc) then
          if(enable_debug_print) then
            call debug_print('Running zero_gradient_bc_t%apply_zero_gradient_primitive_var_bc() +y', &
                            __FILE__, __LINE__)
          end if
          do i = 1, self%n_ghost_layers
            rho%data(:, top_ghost(i)) = rho%data(:, top)
            u%data(:, top_ghost(i)) = u%data(:, top)
            v%data(:, top_ghost(i)) = v%data(:, top)
            p%data(:, top_ghost(i)) = p%data(:, top)
          end do
        end if
      case('-y')
        if(rho%on_jlo_bc) then
          if(enable_debug_print) then
            call debug_print('Running zero_gradient_bc_t%apply_zero_gradient_primitive_var_bc() -y', &
                            __FILE__, __LINE__)
          end if
          do i = 1, self%n_ghost_layers
            rho%data(:, bottom_ghost(i)) = rho%data(:, bottom)
            u%data(:, bottom_ghost(i)) = u%data(:, bottom)
            v%data(:, bottom_ghost(i)) = v%data(:, bottom)
            p%data(:, bottom_ghost(i)) = p%data(:, bottom)
          end do
        end if
      case default
        error stop "Unsupported location to apply the bc at in zero_gradient_bc_t%apply_zero_gradient_cell_gradient_bc()"
      end select
    end associate

  end subroutine apply_zero_gradient_primitive_var_bc

  subroutine finalize(self)
    type(zero_gradient_bc_t), intent(inout) :: self
    if(enable_debug_print) call debug_print('Running zero_gradient_bc_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%ilo_ghost)) deallocate(self%ilo_ghost)
    if(allocated(self%ihi_ghost)) deallocate(self%ihi_ghost)
    if(allocated(self%jlo_ghost)) deallocate(self%jlo_ghost)
    if(allocated(self%jhi_ghost)) deallocate(self%jhi_ghost)
  end subroutine finalize

end module mod_zero_gradient_bc
