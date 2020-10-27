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
  use mod_field, only: field_2d_t
  use mod_input, only: input_t
  use mod_grid_block, only: grid_block_t

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
    integer(ik), dimension(:), allocatable :: ilo_ghost !< (innermost:outermost ghost_layer); ghost layer indices; note: innermost (closest to real domain) layer is always first
    integer(ik), dimension(:), allocatable :: ihi_ghost !< (innermost:outermost ghost_layer); ghost layer indices; note: innermost (closest to real domain) layer is always first
    integer(ik), dimension(:), allocatable :: jlo_ghost !< (innermost:outermost ghost_layer); ghost layer indices; note: innermost (closest to real domain) layer is always first
    integer(ik), dimension(:), allocatable :: jhi_ghost !< (innermost:outermost ghost_layer); ghost layer indices; note: innermost (closest to real domain) layer is always first
    integer(ik) :: ilo = 0 !< ilo real domain cell index
    integer(ik) :: ihi = 0 !< ihi real domain cell index
    integer(ik) :: jlo = 0 !< jlo real domain cell index
    integer(ik) :: jhi = 0 !< jhi real domain cell index
  contains
    procedure, public :: set_time
    procedure, public :: get_time
    procedure, public :: set_indices
    procedure(apply), public, deferred :: apply

  end type boundary_condition_t

  abstract interface
    subroutine apply(self, rho, u, v, p)
      import :: boundary_condition_t, field_2d_t
      import :: ik, rk
      class(boundary_condition_t), intent(inout) :: self
      class(field_2d_t), intent(inout) :: rho
      class(field_2d_t), intent(inout) :: u
      class(field_2d_t), intent(inout) :: v
      class(field_2d_t), intent(inout) :: p
    end subroutine apply
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
    class(grid_block_t), intent(in) :: grid
    integer(ik) :: i
    ! integer(ik), dimension(:, :), intent(in) :: ghost_layers
    !< (ilo_layers(n), ihi_layers(n), jlo_layers(n), jhi_layers(n)); indices to the ghost layers.
    !< The ilo_layers type var can be scalar or a vector, e.g. ilo_layers = [-1,0] or ilo_layers = 0

    ! Note: The indexing order starts with the inner-most index, e.g the index closest to the real domain

    ! Example:
    ! the cells are numbered as follows: [-1 0 | 1 2 3 4 | 5 6]
    ! the real cells are [1,2,3,4] and the ghost cells are [-1 0] and [5 6].
    ! Then ilo_ghost = [0, -1] and ihi_ghost = [5, 6]

    self%n_ghost_layers = grid%n_halo_cells

    allocate(self%ilo_ghost(self%n_ghost_layers))
    allocate(self%ihi_ghost(self%n_ghost_layers))
    allocate(self%jlo_ghost(self%n_ghost_layers))
    allocate(self%jhi_ghost(self%n_ghost_layers))

    if(self%n_ghost_layers == 1) then
      self%ilo_ghost = grid%cell_lbounds_halo(1)
      self%ihi_ghost = grid%cell_ubounds_halo(1)
      self%jlo_ghost = grid%cell_lbounds_halo(2)
      self%jhi_ghost = grid%cell_ubounds_halo(2)
      self%ilo = grid%cell_lbounds(1)
      self%ihi = grid%cell_ubounds(1)
      self%jlo = grid%cell_lbounds(2)
      self%jhi = grid%cell_ubounds(2)
    else if(self%n_ghost_layers > 1) then
      do i = 1, self%n_ghost_layers
        self%ilo_ghost(i) = grid%cell_lbounds_halo(1) + i - 1
        self%ihi_ghost(i) = grid%cell_ubounds_halo(1) - i + 1

        self%jlo_ghost(i) = grid%cell_lbounds_halo(2) + i - 1
        self%jhi_ghost(i) = grid%cell_ubounds_halo(2) - i + 1
      end do
      self%ilo_ghost = .reverse.self%ilo_ghost
      self%jlo_ghost = .reverse.self%jlo_ghost
      self%ihi_ghost = .reverse.self%ihi_ghost
      self%jhi_ghost = .reverse.self%jhi_ghost

      self%ilo = grid%cell_lbounds(1)
      self%ihi = grid%cell_ubounds(1)
      self%jlo = grid%cell_lbounds(2)
      self%jhi = grid%cell_ubounds(2)
    else
      error stop "Error in boundary_condition_t%set_indicies(), n_ghost_layers <= 0"
    end if

  end subroutine set_indices
end module mod_boundary_conditions
