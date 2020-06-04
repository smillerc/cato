module mod_bc_factory

  use mod_input, only: input_t
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_periodic_bc, only: periodic_bc_t, periodic_bc_constructor
  ! use mod_reflection_bc, only: reflection_bc_t, reflection_bc_constructor
  use mod_pressure_input_bc, only: pressure_input_bc_t, pressure_input_bc_constructor
  use mod_zero_gradient_bc, only: zero_gradient_bc_t, zero_gradient_bc_constructor
  ! use mod_vacuum_bc, only: vacuum_bc_t, vacuum_bc_constructor
  implicit none

  private
  public :: bc_factory
contains

  function bc_factory(bc_type, location, input) result(bc)
    !< Factory function to create a boundary condition object
    character(len=*), intent(in) :: bc_type
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(boundary_condition_t), pointer :: bc
    class(input_t), intent(in) :: input

    select case(trim(bc_type))
    case('periodic')
      bc => periodic_bc_constructor(location, input)
      bc%priority = 0
      ! case('reflection')
      !   bc => reflection_bc_constructor(location, input)
      !   bc%priority = 1
    case('pressure_input')
      bc => pressure_input_bc_constructor(location, input)
      bc%priority = 2
    case('zero_gradient')
      bc => zero_gradient_bc_constructor(location, input)
      bc%priority = 1
      ! case('vacuum')
      !   bc => vacuum_bc_constructor(location, input)
      !   bc%priority = 1
    case default
      write(*, '(3(a))') "Unsupported boundary condition type in bc_factory: '", trim(bc_type), "'"
      error stop "Unsupported boundary condition type in bc_factory"
    end select

  end function bc_factory

end module mod_bc_factory
