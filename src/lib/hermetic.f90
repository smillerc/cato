module mod_hermetic
  use, intrinsic :: iso_fortran_env, only: std_out => output_unit, std_err => error_unit
  use mod_globals, only: enable_debug_print

  private
  public :: hermetic ! Expose type and type-bound procedures

  type, abstract :: hermetic
    private
    integer, pointer :: temporary => null() !Null marks non-temporary data
  contains
    procedure :: set_temp   ! Mark object as temporary
    procedure :: guard_temp ! Increment the depth count
    procedure, pass(self) :: clean_temp ! Decrement depth count/Free memory if 1
    procedure(final_interface), deferred :: force_finalization
  end type

  abstract interface
    subroutine final_interface(self)
      import :: hermetic
      class(hermetic), intent(inout) :: self
    end subroutine
  end interface
contains

  subroutine set_temp(self, calling_function, line)
    class(hermetic), intent(inout) :: self
    character(len=*), intent(in) :: calling_function
    integer, intent(in) :: line
    integer, target :: val
    val = 1

    if(.not. associated(self%temporary)) then
      allocate(self%temporary, source=1)
      ! self%temporary => val
    end if

    if(enable_debug_print) then
      write(std_out, '(a,l1)') calling_function//' -> set_temp(), associated(self%temporary) = ', associated(self%temporary)
    end if

  end subroutine

  subroutine guard_temp(self, calling_function, line)
    class(hermetic) :: self
    character(len=*), intent(in) :: calling_function
    integer, intent(in) :: line

    if(enable_debug_print) then
      write(std_out, '(a,l1)') calling_function//' -> guard_temp(), associated(self%temporary) = ', associated(self%temporary)
    end if

    if(associated(self%temporary)) then
      self%temporary = self%temporary + 1
      ! write(std_out, '(a,i0)') calling_function // ' -> guard_temp(), self%temporary = ', self%temporary
    end if
  end subroutine

  subroutine clean_temp(self, calling_function, line)
    class(hermetic) :: self
    character(len=*), intent(in) :: calling_function
    integer, intent(in) :: line

    if(enable_debug_print) then
      write(std_out, '(a,l1)') calling_function//' -> clean_temp(), associated(self%temporary) = ', associated(self%temporary)
    end if

    if(associated(self%temporary)) then

      if(self%temporary > 1) then
        self%temporary = self%temporary - 1
      end if

      if(self%temporary == 1) then

        if(enable_debug_print) then
          write(std_out, '(a)') calling_function//' -> cleaning up!'
        end if

        call self%force_finalization()
        deallocate(self%temporary)
      end if
    end if
  end subroutine
end module mod_hermetic
