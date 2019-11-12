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

    procedure, public :: add => add_fvleg
    procedure, public :: multiply => multiply_fvleg
    procedure, public :: assign => assign_fvleg
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
    class(fvleg_t), intent(in) :: self
    class(integrand), allocatable :: dState_dt
  end function time_derivative

  function add_fvleg(lhs, rhs) result(operator_result)
    class(fvleg_t), intent(in) :: lhs
    class(integrand), intent(in) :: rhs
    class(integrand), allocatable :: operator_result
  end function add_fvleg

  function multiply_fvleg(lhs, rhs) result(operator_result)
    class(fvleg_t), intent(in) :: lhs
    class(integrand), allocatable :: operator_result
    real(rk), intent(in) :: rhs
  end function multiply_fvleg

  subroutine assign_fvleg(lhs, rhs)
    class(fvleg_t), intent(inout) :: lhs
    class(integrand), intent(in) :: rhs
  end subroutine assign_fvleg

end module mod_fvleg
