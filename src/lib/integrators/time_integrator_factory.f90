module mod_time_integrator_factory

  use mod_strategy, only: strategy
  use mod_input, only: input_t
  use mod_2nd_order_runge_kutta, only: runge_kutta_2nd

  implicit none

  private
  public :: time_integrator_factory

contains
  function time_integrator_factory(input) result(integrator)
    !! Factory function to create a time integration strategy object
    class(input_t), intent(in) :: input
    class(strategy), allocatable :: integrator

    select case(input%time_integration_strategy)

    case('rk2')
      allocate(runge_kutta_2nd :: integrator)
    case default
      error stop 'Error in time_integrator_factory - unsupported time integration scheme'
    end select

  end function time_integrator_factory
end module mod_time_integrator_factory
