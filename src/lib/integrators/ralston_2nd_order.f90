module mod_2nd_order_ralston

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_surrogate, only: surrogate
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t
  use mod_finite_volume_schemes, only: finite_volume_scheme_t

  implicit none
  private
  public :: ralston_2nd

  type, extends(strategy) :: ralston_2nd
    !< Ralston's 2nd-order time integration
  contains
    procedure, nopass :: integrate ! integration procedure
  end type

contains

  subroutine integrate(U, finite_volume_scheme, dt)
    !< Time integrator implementation
    class(surrogate), intent(inout) :: U
    class(finite_volume_scheme_t), intent(in) :: finite_volume_scheme
    real(rk), intent(in) :: dt
    class(integrand_t), allocatable :: U_half !< function evaluation at interval t+dt/2

    call debug_print('Running ralston_2nd%integrate()', __FILE__, __LINE__)

    select type(U)
      ! class is(integrand_t)
      !   allocate(U_half, source=U)
      !   U_half = U + U%t(finite_volume_scheme) * dt
      !   U = 0.5_rk * U + 0.5_rk * (U_half + (U_half%t(finite_volume_scheme) * dt))
      !   deallocate(U_half)
    class default
      error stop 'Error in ralston_2nd%integrate - unsupported class'
    end select

  end subroutine
end module mod_2nd_order_ralston
