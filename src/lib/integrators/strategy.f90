module mod_strategy
  !< Substitute for integrand (avoiding circular references)

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_surrogate, only: surrogate

  implicit none
  private
  public :: strategy

  type, abstract :: strategy
    !< Abstract time integration strategy
  contains

    procedure(integrator_interface), nopass, deferred :: integrate !< integration procedure interface
  end type

  abstract interface
    subroutine integrator_interface(self, dt)
      import :: surrogate
      import :: rk
      class(surrogate), intent(inout) :: self
      real(rk), intent(in) :: dt !< time step size
    end subroutine
  end interface

end module
