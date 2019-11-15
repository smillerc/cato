module mod_reconstruction_factory

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_first_order_reconstruction, only: first_order_reconstruction_t
  use mod_second_order_reconstruction, only: second_order_reconstruction_t

  implicit none

  private
  public :: reconstruction_factory_t

  type reconstruction_factory_t
    class(abstract_reconstruction_t), pointer :: reconstruction_ptr => null()
  contains
    private
    ! procedure, public :: init
    procedure, public :: create_reconstruction
    final :: final_factory
  end type reconstruction_factory_t

contains

  ! subroutine init(self, input)
  !   class(reconstruction_factory_t), intent(inout) :: self
  !   class(input_t), intent(in) :: input
  !   self%reconstruction_ptr => null()
  ! end subroutine init

  subroutine final_factory(self)
    type(reconstruction_factory_t), intent(inout) :: self
    nullify(self%reconstruction_ptr)
  end subroutine final_factory

  function create_reconstruction(self, input) result(ptr)
    !< Factory function to create a pointer to a reconstruction operator
    class(reconstruction_factory_t) :: self
    class(input_t), intent(in) :: input

    class(abstract_reconstruction_t), pointer :: ptr

    select case(input%reconstruction_type)
    case('0th_order')
      allocate(first_order_reconstruction_t :: self%reconstruction_ptr)
      call self%reconstruction_ptr%initialize()
      ptr => self%reconstruction_ptr

    case('piecewise_linear')
      allocate(second_order_reconstruction_t :: self%reconstruction_ptr)
      call self%reconstruction_ptr%initialize()
      ptr => self%reconstruction_ptr

    case default
      if(this_image() == 1) then
        error stop "Error in reconstruction_factory_t, unrecognizable reconstruction type"
      end if
    end select

  end function create_reconstruction

end module mod_reconstruction_factory
