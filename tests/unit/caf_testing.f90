module caf_testing
  use iso_fortran_env
  implicit none(type, external)

  real(real64), parameter :: eps_dp = epsilon(1.0_real64)
  real(real64), parameter :: eps_sp = epsilon(1.0_real32)

  interface assert_equal
    procedure :: assert_equal_real64_0d_0d
    procedure :: assert_equal_real64_0d_1d, assert_equal_real64_1d_1d
    procedure :: assert_equal_real64_0d_2d, assert_equal_real64_2d_2d
    procedure :: assert_equal_real64_0d_3d, assert_equal_real64_3d_3d

    procedure :: assert_equal_int32_0d_0d
  end interface

contains

  subroutine int_fail_message()
  end subroutine int_fail_message

  subroutine real_fail_message()
  end subroutine real_fail_message

  subroutine assert_equal_real64_0d_0d(desired, actual, tol)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual
    real(real64), intent(in), optional :: tol
    
    real(real64) :: eps
    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = eps_dp
    end if

    if(abs(desired - actual) < eps) correct = .true.
    if(.not. correct) then
      print*, "assert_equal (=) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if
  end subroutine assert_equal_real64_0d_0d

  subroutine assert_equal_real64_0d_1d(desired, actual, tol)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:)
    real(real64), intent(in), optional :: tol
    
    real(real64) :: eps
    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = eps_dp
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.
    if(.not. correct) then
      print*, "assert_equal (=) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if
  end subroutine assert_equal_real64_0d_1d

  subroutine assert_equal_real64_1d_1d(desired, actual, tol)
    real(real64), intent(in) :: desired(:)
    real(real64), intent(in) :: actual(:)
    real(real64), intent(in), optional :: tol
    
    real(real64) :: eps
    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = eps_dp
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.
    if(.not. correct) then
      print*, "assert_equal (=) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if
  end subroutine assert_equal_real64_1d_1d

  subroutine assert_equal_real64_0d_2d(desired, actual, tol)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:, :)
    real(real64), intent(in), optional :: tol
    
    real(real64) :: eps
    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = eps_dp
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.
    if(.not. correct) then
      print*, "assert_equal (=) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if
  end subroutine assert_equal_real64_0d_2d

  subroutine assert_equal_real64_2d_2d(desired, actual, tol)
    real(real64), intent(in) :: desired(:, :)
    real(real64), intent(in) :: actual(:, :)
    real(real64), intent(in), optional :: tol
    
    real(real64) :: eps
    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = eps_dp
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.
    if(.not. correct) then
      print*, "assert_equal (=) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if
  end subroutine assert_equal_real64_2d_2d

  subroutine assert_equal_real64_0d_3d(desired, actual, tol)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:, :, :)
    real(real64), intent(in), optional :: tol
    
    real(real64) :: eps
    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = eps_dp
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.
    if(.not. correct) then
      print*, "assert_equal (=) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if
  end subroutine assert_equal_real64_0d_3d

  subroutine assert_equal_real64_3d_3d(desired, actual, tol)
    real(real64), intent(in) :: desired(:, :, :)
    real(real64), intent(in) :: actual(:, :, :)
    real(real64), intent(in), optional :: tol
    
    real(real64) :: eps
    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = eps_dp
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.
    if(.not. correct) then
      print*, "assert_equal (=) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if
  end subroutine assert_equal_real64_3d_3d

  subroutine assert_equal_int32_0d_0d(desired, actual)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual
    
    logical :: correct
    correct = .false.

    if(desired == actual) correct = .true.
    if(.not. correct) then
      print*, "assert_equal (=) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if
  end subroutine assert_equal_int32_0d_0d

  subroutine assert_equal_int32_0d_1d(desired, actual)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:)
    
    logical :: correct
    correct = .false.

    if(all(desired == actual)) correct = .true.
    if(.not. correct) then
      print*, "assert_equal (=) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if
  end subroutine assert_equal_int32_0d_1d

  subroutine assert_equal_int32_1d_1d(desired, actual)
    integer(int32), intent(in) :: desired(:)
    integer(int32), intent(in) :: actual(:)
    
    logical :: correct
    correct = .false.

    if(all(desired == actual)) correct = .true.
    if(.not. correct) then
      print*, "assert_equal (=) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if
  end subroutine assert_equal_int32_1d_1d

  subroutine assert_equal_int32_0d_2d(desired, actual)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:, :)
    
    logical :: correct
    correct = .false.

    if(all(desired == actual)) correct = .true.
    if(.not. correct) then
      print*, "assert_equal (=) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if
  end subroutine assert_equal_int32_0d_2d
end module