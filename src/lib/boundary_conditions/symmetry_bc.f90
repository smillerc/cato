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
  use mod_grid, only: grid_t
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

  ! These are not type-bound, b/c they are no-pass functions, e.g. they don't use "self"
  interface operator(.iflip.)
    procedure :: flip_on_right_left_edge
  end interface

  interface operator(.jflip.)
    procedure :: flip_on_top_bottom_edge
  end interface

contains

  pure function flip_on_right_left_edge(recon_data) result(flipped)
    !< Helper function to flip the reconstructed cell data along the right (+i) or left (-i) edge
    real(rk), dimension(:), intent(in) :: recon_data !< reconstructed cell data
    real(rk), dimension(size(recon_data)) :: flipped

    select case(size(recon_data))
    case(4)

      !     recon_data              flipped
      !  |-----E3-----|    |    |-----E3-----|
      !  |            |    |    |            |
      !  E4          E2   ->    E2          E4
      !  |            |    |    |            |
      !  |-----E1-----|    |    |-----E1-----|
      !                    |
      !       flip (mirror) on this edge
      !
      !                    or
      !
      !      flipped              recon_data
      !  |-----E3-----|    |    |-----E3-----|
      !  |            |    |    |            |
      !  E2          E4   <-    E4          E2
      !  |            |    |    |            |
      !  |-----E1-----|    |    |-----E1-----|
      !                    |
      !       flip (mirror) on this edge

      flipped(1) = recon_data(1)
      flipped(2) = recon_data(4)
      flipped(3) = recon_data(3)
      flipped(4) = recon_data(2)
    case(8)

      !     recon_data              flipped
      !  P7----P6----P5    |    P5----P6----P7
      !  |            |    |    |            |
      !  P8          P4   ->    P4          P8
      !  |            |    |    |            |
      !  P1----P2----P3    |    P3----P2----P1
      !                    |
      !       flip (mirror) on this edge
      !
      !                   or
      !
      !      flipped              recon_data
      !  P5----P6----P7    |    P7----P6----P5
      !  |            |    |    |            |
      !  P4          P8   <-    P8          P4
      !  |            |    |    |            |
      !  P3----P2----P1    |    P1----P2----P3
      !                    |
      !       flip (mirror) on this edge

      flipped(1) = recon_data(3)
      flipped(2) = recon_data(2)
      flipped(3) = recon_data(1)
      flipped(4) = recon_data(8)
      flipped(5) = recon_data(7)
      flipped(6) = recon_data(6)
      flipped(7) = recon_data(5)
      flipped(8) = recon_data(4)
    case default
      error stop "Error in mod_symmetry_bc%flip_on_right_edge(); only reconstructed data that is of size 4 or 8 is handled"
    end select
  end function flip_on_right_left_edge

  pure function flip_on_top_bottom_edge(recon_data) result(flipped)
    !< Helper function to flip the reconstructed cell data along the top (+j) or bottom (-j) edge
    real(rk), dimension(:), intent(in) :: recon_data !< reconstructed cell data
    real(rk), dimension(size(recon_data)) :: flipped

    select case(size(recon_data))
    case(4)

      !
      !      flipped
      !  |-----E1-----|
      !  |            |
      !  E4          E2
      !  |            |
      !  |-----E3-----|
      !
      !  -------------- flip on this edge
      !
      !  |-----E3-----|
      !  |            |
      !  E4          E2
      !  |            |
      !  |-----E1-----|
      !    recon_data
      !

      flipped(1) = recon_data(3)
      flipped(2) = recon_data(2)
      flipped(3) = recon_data(1)
      flipped(4) = recon_data(4)
    case(8)

      !
      !     flipped
      !  P1----P2----P3
      !  |            |
      !  P8          P4
      !  |            |
      !  P7----P6----P5
      !
      !  -------------- flip on this edge
      !
      !  P7----P6----P5
      !  |            |
      !  P8          P4
      !  |            |
      !  P1----P2----P3
      !    recon_data
      !

      flipped(1) = recon_data(7)
      flipped(2) = recon_data(6)
      flipped(3) = recon_data(5)
      flipped(4) = recon_data(4)
      flipped(5) = recon_data(3)
      flipped(6) = recon_data(2)
      flipped(7) = recon_data(1)
      flipped(8) = recon_data(8)
    case default
      error stop "Error in mod_symmetry_bc%flip_on_top_bottom_edge(); only reconstructed data that is of size 4 or 8 is handled"
    end select
  end function flip_on_top_bottom_edge

  function symmetry_bc_constructor(location, input, grid) result(bc)
    type(symmetry_bc_t), pointer :: bc
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(input_t), intent(in) :: input
    class(grid_t), intent(in) :: grid

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
        if(enable_debug_print) call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() +x', __FILE__, __LINE__)

        ! ghost layer indexing is always innermost to outermost, e.g. g=1 is right next to the real domain
        do g = 1, self%n_ghost_layers
          rho%data(right_ghost(g), :) = rho%data(right - (g - 1), :)
          u%data(right_ghost(g), :) = -u%data(right - (g - 1), :)
          v%data(right_ghost(g), :) = v%data(right - (g - 1), :)
          p%data(right_ghost(g), :) = p%data(right - (g - 1), :)
        end do

      case('-x')
        if(enable_debug_print) call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() -x', __FILE__, __LINE__)

        ! ghost layer indexing is always innermost to outermost, e.g. g=1 is right next to the real domain
        do g = 1, self%n_ghost_layers
          rho%data(left_ghost(g), :) = rho%data(left + (g - 1), :)
          u%data(left_ghost(g), :) = -u%data(left + (g - 1), :)
          v%data(left_ghost(g), :) = v%data(left + (g - 1), :)
          p%data(left_ghost(g), :) = p%data(left + (g - 1), :)
        end do

      case('+y')
        if(enable_debug_print) call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() +y', __FILE__, __LINE__)

        ! ghost layer indexing is always innermost to outermost, e.g. g=1 is right next to the real domain
        do g = 1, self%n_ghost_layers
          rho%data(:, top_ghost(g)) = rho%data(:, top - (g - 1))
          u%data(:, top_ghost(g)) = u%data(:, top - (g - 1))
          v%data(:, top_ghost(g)) = -v%data(:, top - (g - 1))
          p%data(:, top_ghost(g)) = p%data(:, top - (g - 1))
        end do

      case('-y')
        if(enable_debug_print) call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() -y', __FILE__, __LINE__)

        ! ghost layer indexing is always innermost to outermost, e.g. g=1 is right next to the real domain
        do g = 1, self%n_ghost_layers
          rho%data(:, bottom_ghost(g)) = rho%data(:, bottom + (g - 1))
          u%data(:, bottom_ghost(g)) = u%data(:, bottom + (g - 1))
          v%data(:, bottom_ghost(g)) = -v%data(:, bottom + (g - 1))
          p%data(:, bottom_ghost(g)) = p%data(:, bottom + (g - 1))
        end do

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
