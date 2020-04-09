module mod_strategy
  !< Substitute for integrand_t (avoiding circular references)

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_surrogate, only: surrogate
  use mod_finite_volume_schemes, only: finite_volume_scheme_t

  implicit none
  private
  public :: strategy

  type, abstract :: strategy
    !< Abstract time integration strategy
  contains

    procedure(integrator_interface), nopass, deferred :: integrate !< integration procedure interface
    ! procedure, public :: write_residuals
  end type

  abstract interface
    subroutine integrator_interface(U, finite_volume_scheme, dt)
      import :: surrogate
      import :: rk
      import :: finite_volume_scheme_t
      class(surrogate), intent(inout) :: U
      class(finite_volume_scheme_t), intent(inout) :: finite_volume_scheme
      real(rk), intent(inout) :: dt !< time step size
    end subroutine
  end interface
contains

end module
