module mod_bc_factory

  use mod_input, only: input_t
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_periodic_bc, only: periodic_bc_t, periodic_bc_constructor
  implicit none

  private
  public :: bc_factory
contains

  function bc_factory(bc_type, location) result(bc)
    !< Factory function to create a boundary condition object
    character(len=*), intent(in) :: bc_type
    character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    class(boundary_condition_t), pointer :: bc

    select case(trim(bc_type))
    case('periodic')
      bc => periodic_bc_constructor(location)

    case default
      error stop "Unsupported boundary condition type in bc_factory"
    end select

  end function bc_factory

end module mod_bc_factory
