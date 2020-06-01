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
    class(finite_volume_scheme_t), intent(inout) :: finite_volume_scheme
    real(rk), intent(inout) :: dt
    class(integrand_t), allocatable :: U_1
    class(integrand_t), allocatable :: dU_dt

    call debug_print('Running ralston_2nd%integrate()', __FILE__, __LINE__)

    select type(U)
    class is(integrand_t)
      allocate(U_1, source=U)
      allocate(dU_dt, source=U)

      ! 1st stage
      call debug_print('Running ralston_2nd 1st stage', __FILE__, __LINE__)
      dU_dt = U%t(finite_volume_scheme, stage=1)
      U_1 = U + (2.0_rk * dt / 3.0_rk) * dU_dt
      call U_1%residual_smoother()

      ! Final stage
      call debug_print('Running ralston_2nd 2nd stage', __FILE__, __LINE__)
      U = (U + (dt / 4.0_rk) * dU_dt) + (3.0_rk * dt / 4.0_rk) * U_1%t(finite_volume_scheme, stage=2)
      call U%residual_smoother()
      deallocate(U_1)
      deallocate(dU_dt)
    class default
      error stop 'Error in ralston_2nd%integrate - unsupported class'
    end select

  end subroutine
end module mod_2nd_order_ralston
