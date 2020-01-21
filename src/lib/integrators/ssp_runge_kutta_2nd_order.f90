module mod_ssp_2nd_order_runge_kutta

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_surrogate, only: surrogate
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t
  use mod_finite_volume_schemes, only: finite_volume_scheme_t

  implicit none
  private
  public :: ssp_runge_kutta_2nd

  type, extends(strategy) :: ssp_runge_kutta_2nd
    !< 2nd-order Heun's Method time integration
  contains
    procedure, nopass :: integrate ! integration procedure
  end type

contains

  subroutine integrate(U, finite_volume_scheme, dt)
    !< Time integrator implementation
    class(surrogate), intent(inout) :: U
    class(finite_volume_scheme_t), intent(in) :: finite_volume_scheme
    real(rk), intent(in) :: dt
    class(integrand_t), allocatable :: U_1 !< first stage

    call debug_print('Running ssp_runge_kutta_2nd%integrate()', __FILE__, __LINE__)

    select type(U)
    class is(integrand_t)
      allocate(U_1, source=U)

      ! 1st stage
      U_1 = U + U%t(finite_volume_scheme) * dt

      ! Final stage
      U = 0.5_rk * U + &
          0.5_rk * U_1 + &
          0.5_rk * U_1%t(finite_volume_scheme) * dt
      ! // TODO: Do I need to have bc's applied at each stage?
      deallocate(U_1)
    class default
      error stop 'Error in ssp_runge_kutta_2nd%integrate - unsupported class'
    end select

  end subroutine
end module mod_ssp_2nd_order_runge_kutta
