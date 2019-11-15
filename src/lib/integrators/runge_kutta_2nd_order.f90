module runge_kutta_2nd_module
  use iso_fortran_env, only: ik => int32, rk => real64

  use mod_surrogate, only: surrogate
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t

  implicit none
  private

  !< 2nd-order Runge-Kutta time integration
  type, extends(strategy), public :: runge_kutta_2nd
  contains
    procedure, nopass :: integrate ! integration procedure
  end type

contains

  subroutine integrate(self, dt)
    !< Time integrator
    class(surrogate), intent(inout) :: self
    real(rk), intent(in) :: dt
    class(integrand_t), allocatable :: self_half !< function evaluation at interval t+dt/2.

    select type(self)
    class is(integrand_t)
      allocate(self_half, source=self)
      self_half = self + self%t() * (0.5 * dt)
      self = self + self_half%t() * dt
    class default
      stop 'integrate: unsupported class'
    end select
  end subroutine
end module
