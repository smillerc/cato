module mod_2nd_order_heun

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_surrogate, only: surrogate
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t
  use mod_finite_volume_schemes, only: finite_volume_scheme_t

  implicit none
  private
  public :: heun_2nd

  type, extends(strategy) :: heun_2nd
    !< 2nd-order Heun's Method time integration
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

    real(rk) :: min_dt = 1.0e-30_rk
    call debug_print('Running heun_2nd%integrate()', __FILE__, __LINE__)

    select type(U)
    class is(integrand_t)
      allocate(U_1, source=U)

      ! 1st stage
      call debug_print('Running heun_2nd 1st stage', __FILE__, __LINE__)
      U_1 = U + U%t(finite_volume_scheme) * dt

      ! Final stage
      call debug_print('Running heun_2nd 2nd stage', __FILE__, __LINE__)
      U = 0.5_rk * U + &
          0.5_rk * U_1 + &
          (0.5_rk * dt) * U_1%t(finite_volume_scheme)
      deallocate(U_1)
    class default
      error stop 'Error in heun_2nd%integrate - unsupported class'
    end select

  end subroutine
end module mod_2nd_order_heun

! subroutine integrate(U, finite_volume_scheme, dt)
!   !< Time integrator implementation
!   class(surrogate), intent(inout) :: U
!   class(finite_volume_scheme_t), intent(inout) :: finite_volume_scheme
!   real(rk), intent(inout) :: dt
!   class(integrand_t), allocatable :: U_1 !< first stage
!   class(integrand_t), allocatable :: U_0 !< first stage
!   class(integrand_t), allocatable :: dU_dt !< first stage

!   real(rk) :: min_dt = 1.0e-30_rk
!   call debug_print('Running heun_2nd%integrate()', __FILE__, __LINE__)

!   select type(U)
!   class is(integrand_t)
!     allocate(U_0, source=U)
!     allocate(U_1, source=U)
!     allocate(dU_dt, source=U)

!     dU_dt = U%t(finite_volume_scheme)
!     if (finite_volume_scheme%error_code /= 0) then
!       write(*,'(a, es12.4)') "Errors in dU/dt, exiting..."
!       error stop
!     end if

!     do
!       finite_volume_scheme%error_code = 0
!       ! 1st stage
!       call debug_print('Running heun_2nd 1st stage', __FILE__, __LINE__)
!       U_1 = U_0 + dU_dt * dt

!       ! Final stage
!       call debug_print('Running heun_2nd 2nd stage', __FILE__, __LINE__)

!       U = U_0 + 0.5_rk * dt * dU_dt + 0.5_rk * dt *  U_1%t(finite_volume_scheme)

!       if (finite_volume_scheme%error_code == 0) then
!         exit
!       else
!         if (dt > min_dt) then
!           dt = dt / 10.0_rk
!           write(*,'(a, es12.4)') "Errors in dU/dt, reducing the timestep to dt =", dt
!         else
!           write(*,'(a, es12.4)') 'Attempted to reduce dt, but now dt < min_dt, dt =', min_dt
!           error stop
!         end if
!       end if
!     end do

!     deallocate(U_0)
!     deallocate(U_1)
!     deallocate(dU_dt)
!   class default
!     error stop 'Error in heun_2nd%integrate - unsupported class'
!   end select

! end subroutine
