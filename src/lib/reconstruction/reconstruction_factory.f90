module mod_reconstruction_factory

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_input, only: input_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  ! use mod_first_order_reconstruction, only: first_order_reconstruction_t
  use mod_second_order_reconstruction, only: second_order_reconstruction_t
  use mod_grid, only: grid_t

  implicit none

  private
  public :: reconstruction_factory

  ! type reconstruction_factory_t
  !   class(abstract_reconstruction_t), allocatable :: reconstruction_ptr
  !   class(input_t), pointer :: input
  !   character(:), allocatable :: reconstruction_type
  ! contains
  !   private
  !   procedure, public :: create_reconstruction
  !   final :: final_factory
  ! end type reconstruction_factory_t

  ! interface reconstruction_factory_t
  !   module procedure :: constructor
  ! end interface

contains

  ! type(reconstruction_factory_t) function constructor(input) result(factory)
  !   class(input_t), intent(in), target :: input
  !   ! factory%reconstruction_ptr => null()
  !   factory%reconstruction_type = trim(input%reconstruction_type)
  !   factory%input => input
  ! end function

  ! subroutine final_factory(self)
  !   type(reconstruction_factory_t), intent(inout) :: self
  !   print *, 'Finalizing reconstruction_factory_t'
  !   ! nullify(self%reconstruction_ptr)
  !   nullify(self%input)
  ! end subroutine final_factory

  ! function create_reconstruction(self, grid) result(ptr)
  !   !< Factory function to create a pointer to a reconstruction operator
  !   class(reconstruction_factory_t) :: self
  !   class(grid_t), intent(in), target :: grid
  !   class(abstract_reconstruction_t), allocatable :: ptr

  !   select case(self%reconstruction_type)
  !     ! case('cell_average')
  !     !   allocate(first_order_reconstruction_t :: self%reconstruction_ptr)
  !     !   call self%reconstruction_ptr%initialize(input)
  !     !   ptr => self%reconstruction_ptr

  !   case('piecewise_linear')
  !     allocate(second_order_reconstruction_t :: ptr)
  !     call ptr%initialize(input=self%input, grid=grid)

  !   case default
  !     ! if(this_image() == 1) then
  !     error stop "Error in reconstruction_factory_t, unrecognizable reconstruction type"
  !     ! end if
  !   end select

  ! end function create_reconstruction

  function reconstruction_factory(input, grid_target) result(operator)
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid_target
    ! integer(ik), dimension(3), intent(in) :: lbounds
    ! real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), &
    !   intent(in), target :: conserved_vars_target

    class(abstract_reconstruction_t), pointer :: operator

    character(len=:), allocatable :: recon_type

    recon_type = trim(input%reconstruction_type)
    call debug_print('Making a "'//recon_type//'" reconstruction operator', __FILE__, __LINE__)

    select case(recon_type)
    case('piecewise_linear')
      allocate(second_order_reconstruction_t :: operator)
      call operator%initialize(input=input, grid_target=grid_target)
    case default
      ! if(this_image() == 1) then
      error stop "Error in reconstruction_factory_t, unrecognizable reconstruction type"
      ! end if
    end select
    ! call debug_print('Done', __FILE__, __LINE__)
  end function

end module mod_reconstruction_factory
