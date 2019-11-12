module mod_integrand
  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_surrogate, only: surrogate
  use mod_strategy, only: strategy

  implicit none
  private

  type, abstract, public, extends(surrogate) :: integrand
    private
    class(strategy), allocatable :: time_integrator
  contains
    procedure, non_overridable :: integrate   ! Time integrator
    procedure, non_overridable :: set_time_integrator
    procedure, non_overridable :: get_time_integrator
    procedure(time_derivative), deferred :: t ! Time derivative that evaluates evolution equations

    ! procedure(symmetric_operator), deferred :: add
    ! procedure(asymmetric_operator), deferred :: multiply
    ! procedure(symmetric_assignment), deferred :: assign

    ! ! Map operators to corresponding procedures
    ! generic :: operator(+) => add
    ! generic :: operator(*) => multiply
    ! generic :: assignment(=) => assign
  end type

  abstract interface
    function time_derivative(self) result(dState_dt)
      import :: integrand
      class(integrand), intent(in)  :: self
      class(integrand), allocatable :: dState_dt
    end function time_derivative
    ! function symmetric_operator(lhs, rhs) result(operator_result)
    !   import :: integrand
    !   class(integrand), intent(in)  :: lhs, rhs
    !   class(integrand), allocatable :: operator_result
    ! end function symmetric_operator
    ! function asymmetric_operator(lhs, rhs) result(operator_result)
    !   import :: integrand
    !   class(integrand), intent(in)  :: lhs
    !   class(integrand), allocatable :: operator_result
    !   real, intent(in)  :: rhs
    ! end function asymmetric_operator
    ! subroutine symmetric_assignment(lhs, rhs)
    !   import :: integrand
    !   class(integrand), intent(in)    :: rhs
    !   class(integrand), intent(inout) :: lhs
    ! end subroutine symmetric_assignment
  end interface

contains

  subroutine set_time_integrator(self, s)
    class(integrand), intent(inout) :: self
    class(strategy), intent(in) :: s

    if(allocated(self%time_integrator)) deallocate(self%time_integrator)

    allocate(self%time_integrator, source=s)
  end subroutine

  function get_time_integrator(self) result(self_strategy)
    class(integrand), intent(in) :: self
    class(strategy), allocatable :: self_strategy

    allocate(self_strategy, source=self%time_integrator)
  end function

  subroutine integrate(model, dt)
    class(integrand) :: model ! integrand
    real(rk), intent(in) :: dt ! time step size
    if(allocated(model%time_integrator)) then
      call model%time_integrator%integrate(model, dt)
    else
      stop 'integrate: no integration procedure available.'
    end if
  end subroutine
end module mod_integrand
