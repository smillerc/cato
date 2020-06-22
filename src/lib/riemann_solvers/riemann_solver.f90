module mod_riemann_solver
  !> Summary: Provide a base Riemann solver class structure
  !> Date: 06/22/2020
  !> Author: Sam Miller
  !> Notes:
  !> References:
  !      [1]

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_grid, only: grid_t
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_bc_factory, only: bc_factory

  implicit none

  private
  public :: riemann_solver_t

  type, abstract :: riemann_solver_t
    type(abstract_reconstruction_t), allocatable :: reconstructor
    class(boundary_condition_t), allocatable :: bc_plus_x
    class(boundary_condition_t), allocatable :: bc_plus_y
    class(boundary_condition_t), allocatable :: bc_minus_x
    class(boundary_condition_t), allocatable :: bc_minus_y
  contains
    procedure(solve), deferred, public :: solve
    procedure, public :: initialize_bcs
  end type riemann_solver_t

  abstract interface
    subroutine solve(self, grid, lbounds, rho, u, v, p, rho_u, rho_v, rho_E)
      import :: rk
      import :: riemann_solver_t
      import :: grid_t
      class(riemann_solver_t), intent(in) :: self
      class(grid_t), intent(in) :: grid
      real(rk), dimension(2), intent(in) :: lbounds
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: rho
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(in) :: u
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(in) :: v
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(in) :: p
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: rho_u
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: rho_v
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: rho_E
    end subroutine
  end interface

contains

  subroutine initialize_bcs()
    class(riemann_solver_t), intent(inout) :: self
    class(boundary_condition_t), pointer :: bc => null()

    allocate(ghost_layers(4, 2))

    ! The ghost_layers array just tells the boundary condition which cell indices
    ! are tagged as boundaries. See boundary_condition_t%set_indices() for details
    ghost_layers(1, :) = [self%grid%ilo_bc_cell, self%grid%ilo_cell - 1] ! ilo
    ghost_layers(3, :) = [self%grid%jlo_bc_cell, self%grid%jlo_cell - 1] ! jlo

    ghost_layers(2, :) = [self%grid%ihi_cell + 1, self%grid%ihi_bc_cell] ! ihi
    ghost_layers(4, :) = [self%grid%jhi_cell + 1, self%grid%jhi_bc_cell] ! jhi

    ! Set boundary conditions
    bc => bc_factory(bc_type=input%plus_x_bc, location='+x', input=input, ghost_layers=ghost_layers)
    allocate(self%bc_plus_x, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%bc_plus_x"
    deallocate(bc)

    bc => bc_factory(bc_type=input%plus_y_bc, location='+y', input=input, ghost_layers=ghost_layers)
    allocate(self%bc_plus_y, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%bc_plus_y"
    deallocate(bc)

    bc => bc_factory(bc_type=input%minus_x_bc, location='-x', input=input, ghost_layers=ghost_layers)
    allocate(self%bc_minus_x, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%bc_minus_x"
    deallocate(bc)

    bc => bc_factory(bc_type=input%minus_y_bc, location='-y', input=input, ghost_layers=ghost_layers)
    allocate(self%bc_minus_y, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate finite_volume_scheme_t%bc_minus_y"
    deallocate(bc)
  end subroutine

end module mod_riemann_solver
