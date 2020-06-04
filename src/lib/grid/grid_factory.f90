module mod_grid_factory

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print, set_domain_dimensionality
  use mod_input, only: input_t
  use mod_grid, only: grid_t
  use mod_regular_2d_grid, only: regular_2d_grid_t, new_regular_2d_grid

  implicit none

  private
  public :: grid_factory

  ! type :: grid_factory
  !   class(grid_t), pointer :: grid_ptr => null()
  ! contains
  !   procedure :: create_grid
  !   final :: finalize
  ! end type
contains

  function grid_factory(input) result(grid)
    !! Factory function to create a grid object
    ! class(grid_factory), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), pointer :: grid

    call debug_print('Creating a grid in grid_factory', __FILE__, __LINE__)

    select case(trim(input%grid_type))
    case('2d_regular')
      grid => new_regular_2d_grid(input)
      call set_domain_dimensionality(dimensionality='2D_XY', grid_orthogonality=.true.)
      ! call grid%initialize(input)

    case default
      error stop 'Unsupported grid type in grid_factory'
    end select
  end function grid_factory

  ! subroutine finalize()
  !   type(grid_factory), intent(inout) :: self
  !   nullify(self%grid_ptr)
  ! end subroutine
end module mod_grid_factory
