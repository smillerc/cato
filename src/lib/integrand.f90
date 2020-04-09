module mod_integrand
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mod_surrogate, only: surrogate
  use mod_globals, only: debug_print
  use mod_strategy, only: strategy
  use mod_local_field, only: local_field_t
  use mod_finite_volume_schemes, only: finite_volume_scheme_t

  implicit none
  private

  type, abstract, public, extends(surrogate) :: integrand_t
    class(strategy), allocatable :: time_integrator
    integer(ik) :: error_code = 0
  contains
    procedure, non_overridable :: integrate   ! Time integrator
    procedure, non_overridable :: set_time_integrator
    procedure, non_overridable :: get_time_integrator
    procedure(time_derivative), deferred :: t ! Time derivative that evaluates evolution equations
    procedure(sanity_check), deferred :: sanity_check
    procedure(global_op_local), deferred :: global_plus_local
    procedure(global_op_local), deferred :: global_minus_local
    ! procedure(local_op_global), pass(rhs), deferred :: local_plus_global
    ! procedure(local_op_global), pass(rhs), deferred :: local_minus_global
    ! procedure(symmetric_operator), deferred :: type_minus_type
    procedure(asymmetric_operator_rhs), pass(rhs), deferred :: real_mul_type
    procedure(asymmetric_operator_lhs), pass(lhs), deferred :: type_mul_real
    procedure(symmetric_assignment), deferred :: assign

    ! Map operators to corresponding procedures
    generic :: operator(+) => global_plus_local
    generic :: operator(-) => global_minus_local
    generic :: operator(*) => real_mul_type, type_mul_real
    generic :: assignment(=) => assign
  end type

  abstract interface
    function time_derivative(self, fv) result(dState_dt)
      !< Implementation of d/dt for the local_field_t
      import :: integrand_t, local_field_t
      import :: finite_volume_scheme_t
      import :: rk
      class(integrand_t), intent(in) :: self
      class(finite_volume_scheme_t), intent(inout) :: fv !< finite volume scheme
      class(local_field_t), allocatable :: dState_dt
    end function time_derivative

    subroutine residual_smoother(self)
      import :: integrand_t
      class(integrand_t), intent(inout) :: self
    end subroutine residual_smoother

    subroutine sanity_check(self, error_code)
      import :: ik, integrand_t, local_field_t
      class(integrand_t), intent(in) :: self
      integer(ik), intent(out) :: error_code
    end subroutine sanity_check

    function symmetric_operator(lhs, rhs) result(operator_result)
      !< LHS +-*/ RHS
      import :: integrand_t, local_field_t
      class(integrand_t), intent(in) :: lhs, rhs
      type(local_field_t) :: operator_result
    end function symmetric_operator

    function global_op_local(lhs, rhs) result(operator_result)
      !< LHS +-*/ RHS
      import :: integrand_t, local_field_t
      class(integrand_t), intent(in) :: lhs
      class(local_field_t), intent(in) :: rhs
      type(local_field_t) :: operator_result
    end function global_op_local

    function local_op_global(lhs, rhs) result(operator_result)
      !< LHS +-*/ RHS
      import :: integrand_t, local_field_t
      class(local_field_t), intent(in) :: lhs
      class(integrand_t), intent(in) :: rhs
      type(local_field_t) :: operator_result
    end function local_op_global

    function asymmetric_operator_lhs(lhs, rhs) result(operator_result)
      !< LHS +-*/ RHS (real64)
      import :: integrand_t, local_field_t
      import :: rk
      real(rk), intent(in) :: rhs
      class(integrand_t), intent(in) :: lhs
      type(local_field_t) :: operator_result
    end function asymmetric_operator_lhs

    function asymmetric_operator_rhs(lhs, rhs) result(operator_result)
      !< LHS (real64) +-*/ RHS
      import :: integrand_t, local_field_t
      import :: rk
      real(rk), intent(in) :: lhs
      class(integrand_t), intent(in) :: rhs
      type(local_field_t) :: operator_result
    end function asymmetric_operator_rhs

    subroutine symmetric_assignment(lhs, rhs)
      !< LHS = RHS
      import :: integrand_t, local_field_t
      class(integrand_t), intent(inout) :: lhs
      class(local_field_t), intent(in) :: rhs
    end subroutine symmetric_assignment
  end interface

contains

  subroutine set_time_integrator(self, s)
    !< Set the time integration scheme, e.g. RK-4, RK-2, etc...
    class(integrand_t), intent(inout) :: self
    class(strategy), intent(in) :: s

    if(allocated(self%time_integrator)) deallocate(self%time_integrator)
    allocate(self%time_integrator, source=s)
  end subroutine

  function get_time_integrator(self) result(self_strategy)
    !< Get the currently set time integration scheme
    class(integrand_t), intent(in) :: self
    class(strategy), allocatable :: self_strategy
    allocate(self_strategy, source=self%time_integrator)
  end function

  subroutine integrate(model, finite_volume_scheme, dt, error_code)
    !< Integration implementation
    class(integrand_t), intent(inout) :: model ! integrand_t
    class(finite_volume_scheme_t), intent(inout) :: finite_volume_scheme
    real(rk), intent(inout) :: dt ! time step size
    integer(ik), intent(out) :: error_code

    if(.not. ieee_is_finite(dt)) then
      error stop 'The timestep "dt" in integrand_t%integrate() is not a finite number'
    end if

    call debug_print('Running integrand_t%integrate', __FILE__, __LINE__)
    if(allocated(model%time_integrator)) then
      call model%time_integrator%integrate(model, finite_volume_scheme, dt)
      ! call model%residual_smoother()
      call model%sanity_check(error_code)
    else
      error stop 'Error: No integration procedure available in integrand_t%integrate()'
    end if
  end subroutine
end module mod_integrand
