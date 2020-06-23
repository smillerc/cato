module mod_flux_solver
  !> Summary: Provide a base Riemann solver class structure
  !> Date: 06/22/2020
  !> Author: Sam Miller
  !> Notes:
  !> References:
  !      [1]

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_grid, only: grid_t
  use mod_input, only: input_t
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_bc_factory, only: bc_factory

  implicit none

  private
  public :: flux_solver_t

  type, abstract :: flux_solver_t
    class(abstract_reconstruction_t), allocatable :: reconstructor
    class(boundary_condition_t), allocatable :: bc_plus_x
    class(boundary_condition_t), allocatable :: bc_plus_y
    class(boundary_condition_t), allocatable :: bc_minus_x
    class(boundary_condition_t), allocatable :: bc_minus_y
    integer(ik) :: iteration = 0
  contains
    procedure(initialize), deferred, public :: initialize
    procedure(solve), deferred, public :: solve
    procedure, public :: init_boundary_conditions
  end type flux_solver_t

  abstract interface
    subroutine solve(self, time, grid, lbounds, rho, u, v, p, rho_u, rho_v, rho_E)
      import :: ik, rk
      import :: flux_solver_t
      import :: grid_t
      class(flux_solver_t), intent(inout) :: self
      class(grid_t), intent(in) :: grid
      integer(ik), dimension(2), intent(in) :: lbounds
      real(rk), intent(in) :: time
      real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: rho
      real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: u
      real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: v
      real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: p
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: rho_u
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: rho_v
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: rho_E
    end subroutine

    subroutine initialize(self, grid, input)
      import :: flux_solver_t
      import :: grid_t
      import :: input_t
      class(flux_solver_t), intent(inout) :: self
      class(grid_t), intent(in), target :: grid
      class(input_t), intent(in) :: input
    end subroutine
  end interface

contains

  subroutine init_boundary_conditions(self, input, grid)
    class(flux_solver_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in) :: grid
    class(boundary_condition_t), pointer :: bc => null()

    ! Locals
    integer(ik) :: alloc_status
    integer(ik), dimension(4, 2) :: ghost_layers

    ! The ghost_layers array just tells the boundary condition which cell indices
    ! are tagged as boundaries. See boundary_condition_t%set_indices() for details
    ghost_layers(1, :) = [grid%ilo_bc_cell, grid%ilo_cell - 1] ! ilo
    ghost_layers(3, :) = [grid%jlo_bc_cell, grid%jlo_cell - 1] ! jlo

    ghost_layers(2, :) = [grid%ihi_cell + 1, grid%ihi_bc_cell] ! ihi
    ghost_layers(4, :) = [grid%jhi_cell + 1, grid%jhi_bc_cell] ! jhi

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

    write(*, '(a)') "Boundary Conditions"
    write(*, '(a)') "==================="
    write(*, '(3(a),i0,a)') "+x: ", trim(self%bc_plus_x%name), ' (priority = ', self%bc_plus_x%priority, ')'
    write(*, '(3(a),i0,a)') "-x: ", trim(self%bc_minus_x%name), ' (priority = ', self%bc_minus_x%priority, ')'
    write(*, '(3(a),i0,a)') "+y: ", trim(self%bc_plus_y%name), ' (priority = ', self%bc_plus_y%priority, ')'
    write(*, '(3(a),i0,a)') "-y: ", trim(self%bc_minus_y%name), ' (priority = ', self%bc_minus_y%priority, ')'
    write(*, *)

  end subroutine

end module mod_flux_solver
