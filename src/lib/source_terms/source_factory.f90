module mod_source_factory

  use mod_input, only: input_t
  use mod_source, only: source_t
  use mod_pressure_source, only: pressure_source_t, new_pressure_source
  implicit none

  private
  public :: source_factory
contains

  function source_factory(input) result(source)
    !< Factory function to create a boundary condition object
    class(source_t), pointer :: source
    class(input_t), intent(in) :: input

    select case(trim(input%source_term_type))
    case('pressure')
      source => new_pressure_source(input)
    case default
      error stop "Unsupported source type in source_factory"
    end select

  end function source_factory

end module mod_source_factory
