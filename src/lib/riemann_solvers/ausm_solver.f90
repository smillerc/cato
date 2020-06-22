module mod_ausm_solver
  !< Summary: Provide a base AUSM solver class structure
  !< Date: 06/22/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !<     [1] M. S. Liou "A sequel to AUSM, Part II AUSM+-up for all speeds",
  !<         Journal of Computational Physics 214 (2006) 137–170

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_riemann_solver, only: riemann_solver_t
  use mod_grid, only: grid_t

  implicit none

  private
  public :: riemann_solver_t

  type, extends(riemann_solver_t) :: ausm_solver_t
  contains
    procedure, public :: solve
    final :: finalize
  end type ausm_solver_t

contains
  subroutine solve(self, grid, lbounds, rho, u, v, p, rho_u, rho_v, rho_E)
    class(ausm_solver_t), intent(in) :: self
    class(grid_t), intent(in) :: grid
    real(rk), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: rho
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in) :: u
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in) :: v
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in) :: p
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: rho_u
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: rho_v
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: rho_E
  end subroutine solve

  subroutine finalize(self)
    !< Class finalizer
    type(ausm_solver_t), intent(inout) :: self
    if(allocated(self%reconstructor)) deallocate(self%reconstructor)
  end subroutine
end module mod_ausm_solver
