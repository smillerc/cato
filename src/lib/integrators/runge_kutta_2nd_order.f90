module mod_2nd_order_runge_kutta

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_surrogate, only: surrogate
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t

  implicit none
  private
  public :: runge_kutta_2nd

  type, extends(strategy) :: runge_kutta_2nd
    !< 2nd-order Runge-Kutta time integration
  contains
    procedure, nopass :: integrate ! integration procedure
  end type

contains

  subroutine integrate(self, dt)
    !< Time integrator implementation
    class(surrogate), intent(inout) :: self
    real(rk), intent(in) :: dt
    class(integrand_t), allocatable :: self_half !< function evaluation at interval t+dt/2

    call debug_print('Running runge_kutta_2nd%integrate()', __FILE__, __LINE__)

    select type(self)
    class is(integrand_t)
      allocate(self_half, source=self)

      print *, '1st step'
      self_half = self - self%t() * dt
      print *

      print *, '2nd step'
      self = 0.5_rk * self + 0.5_rk * (self_half - (self_half%t() * dt))
      print *

      deallocate(self_half)
    class default
      error stop 'Error in runge_kutta_2nd%integrate - unsupported class'
    end select

  end subroutine
end module mod_2nd_order_runge_kutta
