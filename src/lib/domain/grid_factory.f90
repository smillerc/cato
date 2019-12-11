module mod_grid_factory

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t
  use mod_grid, only: grid_t
  use mod_regular_2d_grid, only: regular_2d_grid_t

  implicit none

  private
  public :: grid_factory

contains

  function grid_factory(input) result(grid)
    !! Factory function to create a grid object
    class(input_t), intent(in) :: input
    class(grid_t), allocatable :: grid

    select case(trim(input%grid_type))

    case('2d_regular')
      allocate(regular_2d_grid_t :: grid)
      call grid%initialize(input)

    case default
      error stop 'Unsupported grid type in grid_factory'
    end select

  end function grid_factory

end module mod_grid_factory
