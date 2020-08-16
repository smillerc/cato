! MIT License
! Copyright (c) 2019 Sam Miller
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

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

    if(.not. associated(self%temporary)) allocate(self%temporary)
    self%temporary = 1

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
      if(enable_debug_print) then
        write(std_out, '(a,i0)') calling_function//' -> guard_temp(), self%temporary = ', self%temporary
      end if
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
        if(enable_debug_print) then
          write(std_out, '(a,l1)') calling_function//' -> clean_temp(),  self%temporary = ', self%temporary
        end if
        self%temporary = self%temporary - 1

      end if

      if(self%temporary == 1) then

        if(enable_debug_print) then
          write(std_out, '(a)') calling_function//' self%temporary = 1 -> cleaning up!'
        end if

        call self%force_finalization()
        deallocate(self%temporary)
      end if
    end if
  end subroutine
end module mod_hermetic
