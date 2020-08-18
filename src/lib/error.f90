module mod_error
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit

  implicit none

  integer(ik), parameter :: ALL_OK = 0
  integer(ik), parameter :: NEG_DENSITY = 1
  integer(ik), parameter :: NEG_PRESSURE = 2
  integer(ik), parameter :: NANS_FOUND = 3

contains
  subroutine error_msg(message, module, class, procedure, file_name, line_number, error_stop)
    !< Standardized error reporting
    character(len=*), intent(in) :: message
    character(len=*), intent(in) :: module
    character(len=*), intent(in), optional :: class
    character(len=*), intent(in) :: procedure
    character(len=*), intent(in) :: file_name
    integer(ik), intent(in) :: line_number
    logical, intent(in), optional :: error_stop

    if(present(class)) then
      write(std_err, '(3(a, i0))') "Error: (image: ", this_image(), "/", num_images(), ") : " // trim(message) // "; in " // module // "::" // class // "%" // procedure // "() in " // file_name // ":", line_number
    else
      write(std_err, '(3(a, i0))') "Error: (image: ", this_image(), "/", num_images(), ") : " // trim(message) // "; in " // module // "::" // procedure // "() in " // file_name // ":", line_number
    end if

    ! Catch-all error code
    if(present(error_stop)) then
      if(error_stop) error stop 1
    else
      error stop 1
    end if
  end subroutine

end module mod_error
