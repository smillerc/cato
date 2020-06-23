module mod_ausm_solver
  !< Summary: Provide a base AUSM solver class structure
  !< Date: 06/22/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !<     [1] M. S. Liou "A sequel to AUSM, Part II AUSM+-up for all speeds",
  !<         Journal of Computational Physics 214 (2006) 137–170

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_flux_solver, only: flux_solver_t
  use mod_grid, only: grid_t
  use mod_input, only: input_t

  implicit none

  private
  public :: ausm_solver_t

  type, extends(flux_solver_t) :: ausm_solver_t
  contains
    procedure, public :: initialize => initialize_ausm
    procedure, public :: solve => solve_ausm
    final :: finalize
  end type ausm_solver_t

contains

  subroutine initialize_ausm(self, grid, input)
    class(ausm_solver_t), intent(inout) :: self
    class(grid_t), intent(in), target :: grid
    class(input_t), intent(in) :: input
  end subroutine initialize_ausm

  subroutine solve_ausm(self, time, grid, lbounds, rho, u, v, p, rho_u, rho_v, rho_E)
    class(ausm_solver_t), intent(inout) :: self
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
  end subroutine solve_ausm

  subroutine finalize(self)
    !< Class finalizer
    type(ausm_solver_t), intent(inout) :: self
    if(allocated(self%reconstructor)) deallocate(self%reconstructor)
    if(allocated(self%bc_plus_x)) deallocate(self%bc_plus_x)
    if(allocated(self%bc_plus_y)) deallocate(self%bc_plus_y)
    if(allocated(self%bc_minus_x)) deallocate(self%bc_minus_x)
    if(allocated(self%bc_minus_y)) deallocate(self%bc_minus_y)
  end subroutine finalize
end module mod_ausm_solver
