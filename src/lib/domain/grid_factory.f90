module mod_grid_factory

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t
  use mod_grid, only: grid_t
  use mod_regular_2d_grid, only: regular_2d_grid_t

  implicit none

  private
  public :: grid_factory_t

  type grid_factory_t
    class(grid_t), pointer :: grid_ptr => null()
  contains
    private
    ! procedure, public :: init
    procedure, public :: create_grid
    final :: final_factory
  end type grid_factory_t

contains

  ! subroutine init(self, input)
  !   class(grid_factory_t), intent(inout) :: self
  !   class(input_t), intent(in) :: input
  !   self%grid_ptr => null()
  ! end subroutine init

  subroutine final_factory(self)
    type(grid_factory_t), intent(inout) :: self
    print *, 'Finalizing grid_factory_t'
    nullify(self%grid_ptr)
  end subroutine final_factory

  function create_grid(self, input) result(ptr)
    !! Factory function to create a pointer to a grid object
    class(grid_factory_t) :: self
    class(input_t), intent(in) :: input
    class(grid_t), pointer :: ptr

    select case(input%grid_type)

    case('2d_regular')
      allocate(regular_2d_grid_t :: self%grid_ptr)
      call self%grid_ptr%initialize(input)
      ptr => self%grid_ptr

    case default
      allocate(regular_2d_grid_t :: self%grid_ptr)
      call self%grid_ptr%initialize(input)
      ptr => self%grid_ptr
    end select

  end function create_grid

end module mod_grid_factory
