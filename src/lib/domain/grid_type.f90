module mod_grid

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t

  implicit none

  private
  public :: grid_t

  type, abstract :: grid_t
  contains
    procedure(initialize), deferred :: initialize
  end type grid_t

  abstract interface
    subroutine initialize(self, input)
      import :: grid_t
      import :: input_t
      class(grid_t), intent(inout) :: self
      class(input_t), intent(in) :: input
    end subroutine
  end interface

end module mod_grid
