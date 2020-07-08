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

module mod_boundary_conditions
  !< Define the based boundary condition class

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_functional, only: operator(.reverse.)
  use mod_input, only: input_t
  use mod_grid, only: grid_t

  implicit none

  private
  public :: boundary_condition_t

  type, abstract :: boundary_condition_t
    character(len=32) :: name = ''
    character(len=2) :: location = '' !< Location (+x, -x, +y, or -y)
    real(rk), private :: time = 0.0_rk ! Solution time (for time dependent bc's)
    real(rk) :: max_time = 0.0_rk !< Max time in source (e.g. stop after this)
    integer(ik) :: priority = 0 !< Certain b.c.'s take priority over others in terms of when they are applied (periodic is last)
    integer(ik) :: io_unit = 0

    integer(ik) :: n_ghost_layers = 0
    integer(ik), dimension(:), allocatable :: ilo_ghost !< ilo ghost boundary cell indices
    integer(ik), dimension(:), allocatable :: ihi_ghost !< ihi ghost boundary cell indices
    integer(ik), dimension(:), allocatable :: jlo_ghost !< jlo ghost boundary cell indices
    integer(ik), dimension(:), allocatable :: jhi_ghost !< jhi ghost boundary cell indices
    integer(ik) :: ilo = 0 !< ilo real domain cell index
    integer(ik) :: ihi = 0 !< ihi real domain cell index
    integer(ik) :: jlo = 0 !< jlo real domain cell index
    integer(ik) :: jhi = 0 !< jhi real domain cell index
  contains
    procedure, public :: set_time
    procedure, public :: get_time
    procedure, public :: set_indices
    procedure(apply_primitive_var_bc), public, deferred :: apply_primitive_var_bc
    procedure(apply_gradient_bc), public, deferred :: apply_gradient_bc
    procedure(apply_reconstructed_state_bc), public, deferred :: apply_reconstructed_state_bc
  end type boundary_condition_t

  abstract interface
    subroutine apply_primitive_var_bc(self, rho, u, v, p, lbounds)
      import :: boundary_condition_t
      import :: ik, rk
      class(boundary_condition_t), intent(inout) :: self
      integer(ik), dimension(2), intent(in) :: lbounds
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: rho
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: u
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: v
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: p
    end subroutine apply_primitive_var_bc

    subroutine apply_gradient_bc(self, grad_x, grad_y, lbounds)
      import :: boundary_condition_t
      import :: ik, rk
      class(boundary_condition_t), intent(inout) :: self
      integer(ik), dimension(2), intent(in) :: lbounds
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: grad_x
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: grad_y
    end subroutine apply_gradient_bc

    subroutine apply_reconstructed_state_bc(self, recon_rho, recon_u, recon_v, recon_p, lbounds)
      import :: boundary_condition_t
      import :: ik, rk
      class(boundary_condition_t), intent(inout) :: self
      integer(ik), dimension(2), intent(in) :: lbounds
      real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_rho
      real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_u
      real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_v
      real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_p
    end subroutine apply_reconstructed_state_bc

  end interface
contains
  pure subroutine set_time(self, time)
    class(boundary_condition_t), intent(inout) :: self
    real(rk), intent(in) :: time
    self%time = time
  end subroutine

  pure real(rk) function get_time(self) result(time)
    class(boundary_condition_t), intent(in) :: self
    time = self%time
  end function

  subroutine set_indices(self, grid)
    !< Save the indices for the cells that are tagged as ghost cells. These indices will be used
    !< by the other procedures that apply the boundary conditions

    class(boundary_condition_t), intent(inout) :: self
    class(grid_t), intent(in) :: grid
    integer(ik) :: i
    ! integer(ik), dimension(:, :), intent(in) :: ghost_layers
    !< (ilo_layers(n), ihi_layers(n), jlo_layers(n), jhi_layers(n)); indices to the ghost layers.
    !< The ilo_layers type var can be scalar or a vector, e.g. ilo_layers = [-1,0] or ilo_layers = 0

    ! Example:
    ! the cells are numbered as follows: [-1 0 | 1 2 3 4 | 5 6]
    ! the real cells are [1,2,3,4] and the ghost cells are [-1 0] and [5 6]
    ! or
    ! the cells are numbered as follows: [ 0 | 1 2 3 4 | 5 ]
    ! the real cells are [1,2,3,4] and the ghost cells are [0] and [5]

    self%n_ghost_layers = grid%ilo_cell - grid%ilo_bc_cell

    allocate(self%ilo_ghost(self%n_ghost_layers))
    allocate(self%ihi_ghost(self%n_ghost_layers))
    allocate(self%jlo_ghost(self%n_ghost_layers))
    allocate(self%jhi_ghost(self%n_ghost_layers))

    if(self%n_ghost_layers == 1) then
      self%ilo_ghost = grid%ilo_bc_cell
      self%ihi_ghost = grid%ihi_bc_cell
      self%jlo_ghost = grid%jlo_bc_cell
      self%jhi_ghost = grid%jhi_bc_cell
      self%ilo = grid%ilo_cell
      self%ihi = grid%ihi_cell
      self%jlo = grid%jlo_cell
      self%jhi = grid%jhi_cell
    else if(self%n_ghost_layers > 1) then
      do i = 1, self%n_ghost_layers
        self%ilo_ghost(i) = grid%ilo_bc_cell + i - 1
        self%ihi_ghost(i) = grid%ihi_bc_cell - i + 1

        self%jlo_ghost(i) = grid%jlo_bc_cell + i - 1
        self%jhi_ghost(i) = grid%jhi_bc_cell - i + 1
      end do
      self%ihi_ghost = .reverse.self%ihi_ghost
      self%jhi_ghost = .reverse.self%jhi_ghost

      self%ilo = grid%ilo_cell
      self%ihi = grid%ihi_cell
      self%jlo = grid%jlo_cell
      self%jhi = grid%jhi_cell
    else
      error stop "Error in boundary_condition_t%set_indicies(), n_ghost_layers <= 0"
    end if

  end subroutine set_indices
end module mod_boundary_conditions
