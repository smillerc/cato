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
  use mod_globals, only: debug_print
  use mod_grid, only: grid_t
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_input, only: input_t

  implicit none

  private
  public :: symmetry_bc_t, symmetry_bc_constructor

  type, extends(boundary_condition_t) :: symmetry_bc_t
  contains
    procedure, public :: apply_primitive_var_bc => apply_symmetry_primitive_var_bc
    procedure, public :: apply_reconstructed_state_bc => apply_symmetry_reconstructed_state_bc
    procedure, public :: apply_gradient_bc
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

  subroutine apply_symmetry_primitive_var_bc(self, rho, u, v, p, lbounds)

    class(symmetry_bc_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: rho
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: u
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: v
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: p

    integer(ik) :: g

    ! error stop "FixME symmetry_bc_t"

    associate(left => self%ilo, right => self%ihi, bottom => self%jlo, top => self%jhi, &
              left_ghost => self%ilo_ghost, right_ghost => self%ihi_ghost, &
              bottom_ghost => self%jlo_ghost, top_ghost => self%jhi_ghost)

      select case(self%location)
      case('+x')
        call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() +x', __FILE__, __LINE__)

        do g = 1, self%n_ghost_layers
          rho(right_ghost(g), :) = rho(right + (g - 1), :)
          u(right_ghost(g), :) = -u(right + (g - 1), :)
          v(right_ghost(g), :) = v(right + (g - 1), :)
          p(right_ghost(g), :) = p(right + (g - 1), :)
        end do

      case('-x')
        call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() -x', __FILE__, __LINE__)

        do g = 1, self%n_ghost_layers
          rho(left_ghost(g), :) = rho(left + (g - 1), :)
          u(left_ghost(g), :) = -u(left + (g - 1), :)
          v(left_ghost(g), :) = v(left + (g - 1), :)
          p(left_ghost(g), :) = p(left + (g - 1), :)
        end do

      case('+y')
        call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() +y', __FILE__, __LINE__)

        do g = 1, self%n_ghost_layers
          rho(:, top_ghost(g)) = rho(:, top + (g - 1))
          u(:, top_ghost(g)) = u(:, top + (g - 1))
          v(:, top_ghost(g)) = -v(:, top + (g - 1))
          p(:, top_ghost(g)) = p(:, top + (g - 1))
        end do

      case('-y')
        call debug_print('Running symmetry_bc_t%apply_symmetry_primitive_var_bc() -y', __FILE__, __LINE__)
        do g = 1, self%n_ghost_layers
          rho(:, bottom_ghost(g)) = rho(:, bottom + (g - 1))
          u(:, bottom_ghost(g)) = u(:, bottom + (g - 1))
          v(:, bottom_ghost(g)) = -v(:, bottom + (g - 1))
          p(:, bottom_ghost(g)) = p(:, bottom + (g - 1))
        end do

      case default
        error stop "Unsupported location to apply the bc at in symmetry_bc_t%apply_symmetry_primitive_var_bc()"
      end select
    end associate

  end subroutine apply_symmetry_primitive_var_bc

  subroutine apply_symmetry_reconstructed_state_bc(self, recon_rho, recon_u, recon_v, recon_p, lbounds)
    !< Apply zero gradient boundary conditions to the reconstructed state vector field

    class(symmetry_bc_t), intent(inout) :: self

    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_rho
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_u
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_v
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_p

    integer(ik) :: i, j, g

    associate(left => self%ilo, right => self%ihi, bottom => self%jlo, top => self%jhi, &
              left_ghost => self%ilo_ghost, right_ghost => self%ihi_ghost, &
              bottom_ghost => self%jlo_ghost, top_ghost => self%jhi_ghost)

      select case(self%location)
      case('+x')
        call debug_print('Running periodic_bc_t%apply_periodic_reconstructed_state_bc() +x', __FILE__, __LINE__)

        do j = self%jlo, self%jhi
          do g = 1, self%n_ghost_layers
            recon_rho(:, right_ghost(g), j) = .iflip.recon_rho(:, right + (g - 1), j)
            recon_u(:, right_ghost(g), j) = -.iflip.recon_u(:, right + (g - 1), j)
            recon_v(:, right_ghost(g), j) = .iflip.recon_v(:, right + (g - 1), j)
            recon_p(:, right_ghost(g), j) = .iflip.recon_p(:, right + (g - 1), j)
          end do
        end do

      case('-x')
        call debug_print('Running periodic_bc_t%apply_periodic_reconstructed_state_bc() -x', __FILE__, __LINE__)

        do j = self%jlo, self%jhi
          do g = 1, self%n_ghost_layers
            recon_rho(:, left_ghost(g), j) = .iflip.recon_rho(:, left + (g - 1), j)
            recon_u(:, left_ghost(g), j) = -.iflip.recon_u(:, left + (g - 1), j)
            recon_v(:, left_ghost(g), j) = .iflip.recon_v(:, left + (g - 1), j)
            recon_p(:, left_ghost(g), j) = .iflip.recon_p(:, left + (g - 1), j)
          end do
        end do
      case('+y')
        call debug_print('Running periodic_bc_t%apply_periodic_reconstructed_state_bc() +y', __FILE__, __LINE__)

        do i = self%ilo, self%ihi
          do g = 1, self%n_ghost_layers
            recon_rho(:, i, top_ghost(g)) = .jflip.recon_rho(:, i, top + (g - 1))
            recon_u(:, i, top_ghost(g)) = .jflip.recon_u(:, i, top + (g - 1))
            recon_v(:, i, top_ghost(g)) = -.jflip.recon_v(:, i, top + (g - 1))
            recon_p(:, i, top_ghost(g)) = .jflip.recon_p(:, i, top + (g - 1))
          end do
        end do
      case('-y')
        call debug_print('Running periodic_bc_t%apply_periodic_reconstructed_state_bc() -y', __FILE__, __LINE__)

        do i = self%ilo, self%ihi
          do g = 1, self%n_ghost_layers
            recon_rho(:, i, bottom_ghost(g)) = .jflip.recon_rho(:, i, bottom + (g - 1))
            recon_u(:, i, bottom_ghost(g)) = .jflip.recon_u(:, i, bottom + (g - 1))
            recon_v(:, i, bottom_ghost(g)) = -.jflip.recon_v(:, i, bottom + (g - 1))
            recon_p(:, i, bottom_ghost(g)) = .jflip.recon_p(:, i, bottom + (g - 1))
          end do
        end do
      case default
        error stop "Unsupported location to apply the bc at in periodic_bc_t%apply_periodic_reconstructed_state_bc()"
      end select
    end associate

  end subroutine apply_symmetry_reconstructed_state_bc

  subroutine apply_gradient_bc(self, grad_x, grad_y, lbounds)
    class(symmetry_bc_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: grad_x
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: grad_y
    integer(ik) :: g

    associate(left => self%ilo, right => self%ihi, bottom => self%jlo, top => self%jhi, &
              left_ghost => self%ilo_ghost, right_ghost => self%ihi_ghost, &
              bottom_ghost => self%jlo_ghost, top_ghost => self%jhi_ghost)

      select case(self%location)
      case('+x')
        call debug_print('Running periodic_bc_t%apply_gradient_bc() +x', __FILE__, __LINE__)
        do g = 1, self%n_ghost_layers
          grad_x(right_ghost(g), :) = -grad_x(left + (g - 1), :)
          grad_y(right_ghost(g), :) = grad_y(left + (g - 1), :)
        end do
      case('-x')
        call debug_print('Running periodic_bc_t%apply_gradient_bc() -x', __FILE__, __LINE__)
        do g = 1, self%n_ghost_layers
          grad_x(left_ghost(g), :) = -grad_x(right + (g - 1), :)
          grad_y(left_ghost(g), :) = grad_y(right + (g - 1), :)
        end do
      case('+y')
        call debug_print('Running periodic_bc_t%apply_gradient_bc() +y', __FILE__, __LINE__)
        do g = 1, self%n_ghost_layers
          grad_x(:, top_ghost(g)) = grad_x(:, bottom + (g - 1))
          grad_y(:, top_ghost(g)) = -grad_y(:, bottom + (g - 1))
        end do
      case('-y')
        call debug_print('Running periodic_bc_t%apply_gradient_bc() -y', __FILE__, __LINE__)
        do g = 1, self%n_ghost_layers
          grad_x(:, bottom_ghost(g)) = grad_x(:, top + (g - 1))
          grad_y(:, bottom_ghost(g)) = -grad_y(:, top + (g - 1))
        end do
      case default
        error stop "Unsupported location to apply the bc at in periodic_bc_t%apply_gradient_bc()"
      end select
    end associate
  end subroutine apply_gradient_bc

  subroutine finalize(self)
    type(symmetry_bc_t), intent(inout) :: self
    call debug_print('Running symmetry_bc_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%ilo_ghost)) deallocate(self%ilo_ghost)
    if(allocated(self%ihi_ghost)) deallocate(self%ihi_ghost)
    if(allocated(self%jlo_ghost)) deallocate(self%jlo_ghost)
    if(allocated(self%jhi_ghost)) deallocate(self%jhi_ghost)
  end subroutine finalize
end module mod_symmetry_bc
