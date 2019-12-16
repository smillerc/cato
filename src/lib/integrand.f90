module mod_integrand
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_surrogate, only: surrogate
  use mod_globals, only: debug_print
  use mod_strategy, only: strategy

  implicit none
  private

  type, abstract, public, extends(surrogate) :: integrand_t
    class(strategy), allocatable :: time_integrator
    logical :: initiated = .false.
  contains
    procedure, non_overridable :: integrate   ! Time integrator
    procedure, non_overridable :: set_time_integrator
    procedure, non_overridable :: get_time_integrator
    procedure(time_derivative), deferred :: t ! Time derivative that evaluates evolution equations

    procedure(symmetric_operator), deferred :: type_plus_type
    procedure(symmetric_operator), deferred :: type_minus_type
    procedure(asymmetric_operator_rhs), pass(rhs), deferred :: real_mul_type
    procedure(asymmetric_operator_lhs), pass(lhs), deferred :: type_mul_real
    procedure(symmetric_assignment), deferred :: assign

    ! Map operators to corresponding procedures
    generic :: operator(+) => type_plus_type
    generic :: operator(-) => type_minus_type
    generic :: operator(*) => real_mul_type, type_mul_real
    generic :: assignment(=) => assign
  end type

  abstract interface
    function time_derivative(self) result(dState_dt)
      import :: integrand_t
      class(integrand_t), intent(in) :: self
      class(integrand_t), allocatable :: dState_dt
    end function time_derivative

    function symmetric_operator(lhs, rhs) result(operator_result)
      import :: integrand_t
      class(integrand_t), intent(in) :: lhs, rhs
      class(integrand_t), allocatable :: operator_result
    end function symmetric_operator

    function asymmetric_operator_lhs(lhs, rhs) result(operator_result)
      import :: integrand_t
      import :: rk
      real(rk), intent(in) :: rhs
      class(integrand_t), intent(in) :: lhs
      class(integrand_t), allocatable :: operator_result
    end function asymmetric_operator_lhs

    function asymmetric_operator_rhs(lhs, rhs) result(operator_result)
      import :: integrand_t
      import :: rk
      real(rk), intent(in) :: lhs
      class(integrand_t), intent(in) :: rhs
      class(integrand_t), allocatable :: operator_result
    end function asymmetric_operator_rhs

    subroutine symmetric_assignment(lhs, rhs)
      import :: integrand_t
      class(integrand_t), intent(in) :: rhs
      class(integrand_t), intent(inout) :: lhs
    end subroutine symmetric_assignment
  end interface

contains

  subroutine set_time_integrator(self, s)
    class(integrand_t), intent(inout) :: self
    class(strategy), intent(in) :: s

    if(allocated(self%time_integrator)) deallocate(self%time_integrator)

    allocate(self%time_integrator, source=s)
  end subroutine

  function get_time_integrator(self) result(self_strategy)
    class(integrand_t), intent(in) :: self
    class(strategy), allocatable :: self_strategy

    allocate(self_strategy, source=self%time_integrator)
  end function

  subroutine integrate(model, dt)
    class(integrand_t), intent(inout) :: model ! integrand_t
    real(rk), intent(in) :: dt ! time step size

    call debug_print('Running integrand_t%integrate', __FILE__, __LINE__)
    if(allocated(model%time_integrator)) then
      call model%time_integrator%integrate(model, dt)
    else
      stop 'integrate: no integration procedure available.'
    end if
  end subroutine
end module mod_integrand
