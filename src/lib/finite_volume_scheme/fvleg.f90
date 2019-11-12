module mod_fvleg

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_surrogate, only: surrogate
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand
  use mod_input, only: input_t
  use mod_finite_volume_schemes, only: finite_volume_scheme_t

  implicit none
  private
  public :: fvleg_t

  type, extends(finite_volume_scheme_t) :: fvleg_t
  contains
    procedure, public :: reconstruct => reconstruct_fvleg
    procedure, public :: eval_fluxes => eval_fvleg_fluxes
    procedure, public :: t => time_derivative

  end type

contains

  type(fvleg_t) function constructor(input) result(self)
    !< Construct the finite volume local evolution Galerkin (fvleg) scheme
    class(input_t), intent(in) :: input
  end function

  subroutine reconstruct_fvleg(self)
    class(fvleg_t), intent(inout) :: self
  end subroutine

  subroutine eval_fvleg_fluxes(self)
    class(fvleg_t), intent(inout) :: self
  end subroutine

  function time_derivative(self) result(dState_dt)
    class(fvleg_t), intent(in)  :: self
    class(integrand), allocatable :: dState_dt
  end function time_derivative

end module mod_fvleg
