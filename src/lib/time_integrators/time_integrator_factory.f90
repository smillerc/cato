module mod_time_integrator_factory

  use mod_strategy, only: strategy
  use mod_globals, only: debug_print
  use mod_input, only: input_t
  use mod_ssp_rk2, only: ssp_rk2_t
  use mod_2nd_order_ralston, only: ralston_2nd
  use mod_ssp_rk3, only: ssp_rk3_t

  implicit none

  private
  public :: time_integrator_factory

contains
  function time_integrator_factory(input) result(integrator)
    !! Factory function to create a time integration strategy object
    class(input_t), intent(in) :: input
    class(strategy), pointer :: integrator

    call debug_print('Making an "'//input%time_integration_strategy//'" time integrator', &
                     __FILE__, __LINE__)

    select case(input%time_integration_strategy)
    case('heun', 'ssp_rk2')
      allocate(ssp_rk2_t :: integrator)
    case('ralston')
      allocate(ralston_2nd :: integrator)
    case('ssp_rk3')
      allocate(ssp_rk3_t :: integrator)
    case default
      error stop 'Error in time_integrator_factory - unsupported time integration scheme'
    end select
  end function time_integrator_factory
end module mod_time_integrator_factory
