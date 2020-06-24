module mod_ausm_solver
  !< Summary: Provide a base AUSM solver class structure
  !< Date: 06/22/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !<     [1] M. S. Liou "A sequel to AUSM, Part II AUSM+-up for all speeds",
  !<         Journal of Computational Physics 214 (2006) 137–170

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_flux_solver, only: flux_solver_t
  use mod_grid, only: grid_t
  use mod_input, only: input_t

  implicit none

  private
  public :: ausm_solver_t

  type, extends(flux_solver_t) :: ausm_solver_t
  contains
    ! Public methods
    procedure, public :: initialize => initialize_ausm
    procedure, public :: solve => solve_ausm
    procedure, public, pass(lhs) :: copy => copy_ausm

    ! Private methods
    procedure, private, nopass :: flux_edges
    procedure, private :: apply_primitive_bc
    final :: finalize

    ! Operators
    generic :: assignment(=) => copy
  end type ausm_solver_t

contains

  subroutine initialize_ausm(self, grid, input)
    class(ausm_solver_t), intent(inout) :: self
    class(grid_t), intent(in), target :: grid
    class(input_t), intent(in) :: input

    call debug_print('Running ausm_solver_t%initialize_ausm()', __FILE__, __LINE__)
  end subroutine initialize_ausm

  subroutine copy_ausm(lhs, rhs)
    !< Implement LHS = RHS
    class(ausm_solver_t), intent(inout) :: lhs
    type(ausm_solver_t), intent(in) :: rhs

    call debug_print('Running ausm_solver_t%copy()', __FILE__, __LINE__)

    ! allocate(lhs%bc_plus_x, source=rhs%bc_plus_x)
    ! allocate(lhs%bc_plus_y, source=rhs%bc_plus_y)
    ! allocate(lhs%bc_minus_x, source=rhs%bc_minus_x)
    ! allocate(lhs%bc_minus_y, source=rhs%bc_minus_y)

    lhs%iteration = rhs%iteration
    lhs%time = rhs%time
    lhs%dt = rhs%dt

  end subroutine copy_ausm

  subroutine solve_ausm(self, dt, grid, lbounds, rho, u, v, p, d_rho_dt, d_rho_u_dt, d_rho_v_dt, d_rho_E_dt)
    !< Solve and flux the edges
    class(ausm_solver_t), intent(inout) :: self
    class(grid_t), intent(in) :: grid
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), intent(in) :: dt
    real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: rho !< density
    real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: u   !< x-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: v   !< y-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: p   !< pressure
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) ::   d_rho_dt    !< d/dt of the density field
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: d_rho_u_dt    !< d/dt of the rhou field
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: d_rho_v_dt    !< d/dt of the rhov field
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: d_rho_E_dt    !< d/dt of the rhoE field

    call debug_print('Running ausm_solver_t%solve_ausm()', __FILE__, __LINE__)

    if(dt < tiny(1.0_rk)) error stop "Error in ausm_solver_t%solve_ausm(), the timestep dt is < tiny(1.0_rk)"
    self%time = self%time + dt
    self%dt = dt
    self%iteration = self%iteration + 1

  end subroutine solve_ausm

  subroutine flux_edges()
  end subroutine flux_edges

  subroutine apply_primitive_bc(self, lbounds, rho, u, v, p, &
                                bc_plus_x, bc_minus_x, bc_plus_y, bc_minus_y)
    class(ausm_solver_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: rho
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: u
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: v
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: p
    class(boundary_condition_t), intent(inout):: bc_plus_x
    class(boundary_condition_t), intent(inout):: bc_plus_y
    class(boundary_condition_t), intent(inout):: bc_minus_x
    class(boundary_condition_t), intent(inout):: bc_minus_y

    integer(ik) :: priority
    integer(ik) :: max_priority_bc !< highest goes first

    call debug_print('Running fvleg_solver_t%apply_primitive_var_bc()', __FILE__, __LINE__)

    max_priority_bc = max(bc_plus_x%priority, bc_plus_y%priority, &
                          bc_minus_x%priority, bc_minus_y%priority)

    do priority = max_priority_bc, 0, -1

      if(bc_plus_x%priority == priority) then
        call bc_plus_x%apply_primitive_var_bc(rho=rho, u=u, v=v, p=p, lbounds=lbounds)
      end if

      if(bc_plus_y%priority == priority) then
        call bc_plus_y%apply_primitive_var_bc(rho=rho, u=u, v=v, p=p, lbounds=lbounds)
      end if

      if(bc_minus_x%priority == priority) then
        call bc_minus_x%apply_primitive_var_bc(rho=rho, u=u, v=v, p=p, lbounds=lbounds)
      end if

      if(bc_minus_y%priority == priority) then
        call bc_minus_y%apply_primitive_var_bc(rho=rho, u=u, v=v, p=p, lbounds=lbounds)
      end if

    end do

  end subroutine apply_primitive_bc

  subroutine finalize(self)
    !< Class finalizer
    type(ausm_solver_t), intent(inout) :: self

    call debug_print('Running ausm_solver_t%finalize()', __FILE__, __LINE__)
    ! if(allocated(self%reconstructor)) deallocate(self%reconstructor)
    ! if(allocated(self%bc_plus_x)) deallocate(self%bc_plus_x)
    ! if(allocated(self%bc_plus_y)) deallocate(self%bc_plus_y)
    ! if(allocated(self%bc_minus_x)) deallocate(self%bc_minus_x)
    ! if(allocated(self%bc_minus_y)) deallocate(self%bc_minus_y)
  end subroutine finalize
end module mod_ausm_solver
