module caf_testing
  !< Summary: Provide a simple coarray-aware assertion module to be used in unit tests
  !< Date: 11/09/2020
  !< Author: Sam Miller
  use iso_fortran_env
  implicit none 

  private
  public :: assert_equal, assert_less, assert_greater

  interface assert_equal
    module procedure assert_r0_equal_r0_real32
    module procedure assert_r0_equal_r1_real32
    module procedure assert_r1_equal_r1_real32
    module procedure assert_r0_equal_r2_real32
    module procedure assert_r2_equal_r2_real32
    module procedure assert_r0_equal_r3_real32
    module procedure assert_r3_equal_r3_real32
    module procedure assert_r0_equal_r4_real32
    module procedure assert_r4_equal_r4_real32
    module procedure assert_r0_equal_r5_real32
    module procedure assert_r5_equal_r5_real32
    module procedure assert_r0_equal_r6_real32
    module procedure assert_r6_equal_r6_real32
    module procedure assert_r0_equal_r0_real64
    module procedure assert_r0_equal_r1_real64
    module procedure assert_r1_equal_r1_real64
    module procedure assert_r0_equal_r2_real64
    module procedure assert_r2_equal_r2_real64
    module procedure assert_r0_equal_r3_real64
    module procedure assert_r3_equal_r3_real64
    module procedure assert_r0_equal_r4_real64
    module procedure assert_r4_equal_r4_real64
    module procedure assert_r0_equal_r5_real64
    module procedure assert_r5_equal_r5_real64
    module procedure assert_r0_equal_r6_real64
    module procedure assert_r6_equal_r6_real64

    module procedure assert_r0_equal_r0_int32
    module procedure assert_r0_equal_r1_int32
    module procedure assert_r1_equal_r1_int32
    module procedure assert_r0_equal_r2_int32
    module procedure assert_r2_equal_r2_int32
    module procedure assert_r0_equal_r3_int32
    module procedure assert_r3_equal_r3_int32
    module procedure assert_r0_equal_r4_int32
    module procedure assert_r4_equal_r4_int32
    module procedure assert_r0_equal_r5_int32
    module procedure assert_r5_equal_r5_int32
    module procedure assert_r0_equal_r6_int32
    module procedure assert_r6_equal_r6_int32
    module procedure assert_r0_equal_r0_int64
    module procedure assert_r0_equal_r1_int64
    module procedure assert_r1_equal_r1_int64
    module procedure assert_r0_equal_r2_int64
    module procedure assert_r2_equal_r2_int64
    module procedure assert_r0_equal_r3_int64
    module procedure assert_r3_equal_r3_int64
    module procedure assert_r0_equal_r4_int64
    module procedure assert_r4_equal_r4_int64
    module procedure assert_r0_equal_r5_int64
    module procedure assert_r5_equal_r5_int64
    module procedure assert_r0_equal_r6_int64
    module procedure assert_r6_equal_r6_int64
    
    module procedure assert_r0_equal_r0_logical
    module procedure assert_r0_equal_r1_logical
    module procedure assert_r1_equal_r1_logical
    module procedure assert_r0_equal_r2_logical
    module procedure assert_r2_equal_r2_logical
    module procedure assert_r0_equal_r3_logical
    module procedure assert_r3_equal_r3_logical
    module procedure assert_r0_equal_r4_logical
    module procedure assert_r4_equal_r4_logical
    module procedure assert_r0_equal_r5_logical
    module procedure assert_r5_equal_r5_logical
    module procedure assert_r0_equal_r6_logical
    module procedure assert_r6_equal_r6_logical
    module procedure assert_equal_str
  end interface assert_equal

  interface assert_less
    module procedure assert_r0_less_r0_real32
    module procedure assert_r0_less_r1_real32
    module procedure assert_r1_less_r1_real32
    module procedure assert_r0_less_r2_real32
    module procedure assert_r2_less_r2_real32
    module procedure assert_r0_less_r3_real32
    module procedure assert_r3_less_r3_real32
    module procedure assert_r0_less_r4_real32
    module procedure assert_r4_less_r4_real32
    module procedure assert_r0_less_r5_real32
    module procedure assert_r5_less_r5_real32
    module procedure assert_r0_less_r6_real32
    module procedure assert_r6_less_r6_real32
    module procedure assert_r0_less_r0_real64
    module procedure assert_r0_less_r1_real64
    module procedure assert_r1_less_r1_real64
    module procedure assert_r0_less_r2_real64
    module procedure assert_r2_less_r2_real64
    module procedure assert_r0_less_r3_real64
    module procedure assert_r3_less_r3_real64
    module procedure assert_r0_less_r4_real64
    module procedure assert_r4_less_r4_real64
    module procedure assert_r0_less_r5_real64
    module procedure assert_r5_less_r5_real64
    module procedure assert_r0_less_r6_real64
    module procedure assert_r6_less_r6_real64

    module procedure assert_r0_less_r0_int32
    module procedure assert_r0_less_r1_int32
    module procedure assert_r1_less_r1_int32
    module procedure assert_r0_less_r2_int32
    module procedure assert_r2_less_r2_int32
    module procedure assert_r0_less_r3_int32
    module procedure assert_r3_less_r3_int32
    module procedure assert_r0_less_r4_int32
    module procedure assert_r4_less_r4_int32
    module procedure assert_r0_less_r5_int32
    module procedure assert_r5_less_r5_int32
    module procedure assert_r0_less_r6_int32
    module procedure assert_r6_less_r6_int32
    module procedure assert_r0_less_r0_int64
    module procedure assert_r0_less_r1_int64
    module procedure assert_r1_less_r1_int64
    module procedure assert_r0_less_r2_int64
    module procedure assert_r2_less_r2_int64
    module procedure assert_r0_less_r3_int64
    module procedure assert_r3_less_r3_int64
    module procedure assert_r0_less_r4_int64
    module procedure assert_r4_less_r4_int64
    module procedure assert_r0_less_r5_int64
    module procedure assert_r5_less_r5_int64
    module procedure assert_r0_less_r6_int64
    module procedure assert_r6_less_r6_int64
    
  end interface assert_less

  interface assert_greater
    module procedure assert_r0_greater_r0_real32
    module procedure assert_r0_greater_r1_real32
    module procedure assert_r1_greater_r1_real32
    module procedure assert_r0_greater_r2_real32
    module procedure assert_r2_greater_r2_real32
    module procedure assert_r0_greater_r3_real32
    module procedure assert_r3_greater_r3_real32
    module procedure assert_r0_greater_r4_real32
    module procedure assert_r4_greater_r4_real32
    module procedure assert_r0_greater_r5_real32
    module procedure assert_r5_greater_r5_real32
    module procedure assert_r0_greater_r6_real32
    module procedure assert_r6_greater_r6_real32
    module procedure assert_r0_greater_r0_real64
    module procedure assert_r0_greater_r1_real64
    module procedure assert_r1_greater_r1_real64
    module procedure assert_r0_greater_r2_real64
    module procedure assert_r2_greater_r2_real64
    module procedure assert_r0_greater_r3_real64
    module procedure assert_r3_greater_r3_real64
    module procedure assert_r0_greater_r4_real64
    module procedure assert_r4_greater_r4_real64
    module procedure assert_r0_greater_r5_real64
    module procedure assert_r5_greater_r5_real64
    module procedure assert_r0_greater_r6_real64
    module procedure assert_r6_greater_r6_real64

    module procedure assert_r0_greater_r0_int32
    module procedure assert_r0_greater_r1_int32
    module procedure assert_r1_greater_r1_int32
    module procedure assert_r0_greater_r2_int32
    module procedure assert_r2_greater_r2_int32
    module procedure assert_r0_greater_r3_int32
    module procedure assert_r3_greater_r3_int32
    module procedure assert_r0_greater_r4_int32
    module procedure assert_r4_greater_r4_int32
    module procedure assert_r0_greater_r5_int32
    module procedure assert_r5_greater_r5_int32
    module procedure assert_r0_greater_r6_int32
    module procedure assert_r6_greater_r6_int32
    module procedure assert_r0_greater_r0_int64
    module procedure assert_r0_greater_r1_int64
    module procedure assert_r1_greater_r1_int64
    module procedure assert_r0_greater_r2_int64
    module procedure assert_r2_greater_r2_int64
    module procedure assert_r0_greater_r3_int64
    module procedure assert_r3_greater_r3_int64
    module procedure assert_r0_greater_r4_int64
    module procedure assert_r4_greater_r4_int64
    module procedure assert_r0_greater_r5_int64
    module procedure assert_r5_greater_r5_int64
    module procedure assert_r0_greater_r6_int64
    module procedure assert_r6_greater_r6_int64
    
  end interface assert_greater


contains

! -------------------------------------------------------
! Character assertions
! -------------------------------------------------------
subroutine assert_equal_str(desired, actual, compare_trimed, file, line)
  character(len=*), intent(in) :: desired
  character(len=*), intent(in) :: actual
  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  logical, optional :: compare_trimed  !< trim the strings before compare?

  logical :: do_trim 
  logical :: correct
  
  correct = .false.
  do_trim = .true. !< defaults to true
  if (present(compare_trimed)) do_trim = compare_trimed

  if (do_trim) then
    if (len_trim(desired) /= len_trim(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): len_trim(desired) /= len_trim(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (trim(desired) == trim(actual)) correct = .true.
  else
    if (len(desired) /= len(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): len(desired) /= len(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (desired == actual) correct = .true.
  endif

  if (.not. correct) then
    write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
    write(*,'(a,i0)') "Image: ", this_image()
    if (do_trim) then
      write(*,'(3(a))') "Desired Value: '", trim(desired), "'"
      write(*,'(3(a))') "Actual Value:  '", trim(actual), "'"
    else
      write(*,'(3(a))') "Desired Value: '", desired, "'"
      write(*,'(3(a))') "Actual Value:  '", actual, "'"
    endif
    error stop "Test failure"
  end if

end subroutine

! -------------------------------------------------------
! Real (=) assertions
! -------------------------------------------------------
  subroutine assert_r0_equal_r0_real32(desired, actual, tol, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual
    real(real32), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real32) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line

    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real32)
    end if

    if(abs(desired - actual) < eps) correct = .true.

    if(.not. correct) then
      write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
      write(*,'(a,i0)') "Image: ", this_image()
      write(*,'(a,es16.6)') "Desired Value: ", desired
      write(*,'(a,es16.6)') "Actual Value:  ", actual
      write(*,'(a,es16.6)') "Difference:    ", desired - actual
      write(*,'(a,es16.6)') "Epsilon:       ", eps
      error stop "Test failure"
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r1_real32(desired, actual, tol, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:)
    real(real32), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real32) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real32)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired - actual(i1)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1)
            write(*,'(a,es16.6)')      "Difference:    ", desired - actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
    end if
  end subroutine 
  
  subroutine assert_r1_equal_r1_real32(desired, actual, tol, file, line)
    real(real32), intent(in) :: desired(:)
    real(real32), intent(in) :: actual(:)
    real(real32), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real32) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.
    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real32)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired(i1) - actual(i1)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired(i1)
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1)
            write(*,'(a,es16.6)')      "Difference:    ", desired(i1) - actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r2_real32(desired, actual, tol, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:,:)
    real(real32), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real32) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real32)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired - actual(i1,i2)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2)
            write(*,'(a,es16.6)')      "Difference:    ", desired - actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r2_equal_r2_real32(desired, actual, tol, file, line)
    real(real32), intent(in) :: desired(:,:)
    real(real32), intent(in) :: actual(:,:)
    real(real32), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real32) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.
    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real32)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired(i1,i2) - actual(i1,i2)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired(i1,i2)
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2)
            write(*,'(a,es16.6)')      "Difference:    ", desired(i1,i2) - actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r3_real32(desired, actual, tol, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:,:,:)
    real(real32), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real32) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real32)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired - actual(i1,i2,i3)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3)
            write(*,'(a,es16.6)')      "Difference:    ", desired - actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r3_equal_r3_real32(desired, actual, tol, file, line)
    real(real32), intent(in) :: desired(:,:,:)
    real(real32), intent(in) :: actual(:,:,:)
    real(real32), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real32) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.
    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real32)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired(i1,i2,i3) - actual(i1,i2,i3)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired(i1,i2,i3)
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3)
            write(*,'(a,es16.6)')      "Difference:    ", desired(i1,i2,i3) - actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r4_real32(desired, actual, tol, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:,:,:,:)
    real(real32), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real32) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real32)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired - actual(i1,i2,i3,i4)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*,'(a,es16.6)')      "Difference:    ", desired - actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r4_equal_r4_real32(desired, actual, tol, file, line)
    real(real32), intent(in) :: desired(:,:,:,:)
    real(real32), intent(in) :: actual(:,:,:,:)
    real(real32), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real32) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.
    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real32)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired(i1,i2,i3,i4) - actual(i1,i2,i3,i4)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired(i1,i2,i3,i4)
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*,'(a,es16.6)')      "Difference:    ", desired(i1,i2,i3,i4) - actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r5_real32(desired, actual, tol, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:,:,:,:,:)
    real(real32), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real32) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real32)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired - actual(i1,i2,i3,i4,i5)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*,'(a,es16.6)')      "Difference:    ", desired - actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r5_equal_r5_real32(desired, actual, tol, file, line)
    real(real32), intent(in) :: desired(:,:,:,:,:)
    real(real32), intent(in) :: actual(:,:,:,:,:)
    real(real32), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real32) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.
    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real32)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired(i1,i2,i3,i4,i5) - actual(i1,i2,i3,i4,i5)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired(i1,i2,i3,i4,i5)
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*,'(a,es16.6)')      "Difference:    ", desired(i1,i2,i3,i4,i5) - actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r6_real32(desired, actual, tol, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:,:,:,:,:,:)
    real(real32), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real32) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real32)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired - actual(i1,i2,i3,i4,i5,i6)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*,'(a,es16.6)')      "Difference:    ", desired - actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r6_equal_r6_real32(desired, actual, tol, file, line)
    real(real32), intent(in) :: desired(:,:,:,:,:,:)
    real(real32), intent(in) :: actual(:,:,:,:,:,:)
    real(real32), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real32) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.
    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real32)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired(i1,i2,i3,i4,i5,i6) - actual(i1,i2,i3,i4,i5,i6)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired(i1,i2,i3,i4,i5,i6)
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*,'(a,es16.6)')      "Difference:    ", desired(i1,i2,i3,i4,i5,i6) - actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r0_real64(desired, actual, tol, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual
    real(real64), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real64) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line

    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real64)
    end if

    if(abs(desired - actual) < eps) correct = .true.

    if(.not. correct) then
      write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
      write(*,'(a,i0)') "Image: ", this_image()
      write(*,'(a,es16.6)') "Desired Value: ", desired
      write(*,'(a,es16.6)') "Actual Value:  ", actual
      write(*,'(a,es16.6)') "Difference:    ", desired - actual
      write(*,'(a,es16.6)') "Epsilon:       ", eps
      error stop "Test failure"
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r1_real64(desired, actual, tol, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:)
    real(real64), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real64) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real64)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired - actual(i1)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1)
            write(*,'(a,es16.6)')      "Difference:    ", desired - actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
    end if
  end subroutine 
  
  subroutine assert_r1_equal_r1_real64(desired, actual, tol, file, line)
    real(real64), intent(in) :: desired(:)
    real(real64), intent(in) :: actual(:)
    real(real64), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real64) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.
    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real64)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired(i1) - actual(i1)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired(i1)
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1)
            write(*,'(a,es16.6)')      "Difference:    ", desired(i1) - actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r2_real64(desired, actual, tol, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:,:)
    real(real64), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real64) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real64)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired - actual(i1,i2)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2)
            write(*,'(a,es16.6)')      "Difference:    ", desired - actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r2_equal_r2_real64(desired, actual, tol, file, line)
    real(real64), intent(in) :: desired(:,:)
    real(real64), intent(in) :: actual(:,:)
    real(real64), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real64) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.
    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real64)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired(i1,i2) - actual(i1,i2)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired(i1,i2)
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2)
            write(*,'(a,es16.6)')      "Difference:    ", desired(i1,i2) - actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r3_real64(desired, actual, tol, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:,:,:)
    real(real64), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real64) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real64)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired - actual(i1,i2,i3)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3)
            write(*,'(a,es16.6)')      "Difference:    ", desired - actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r3_equal_r3_real64(desired, actual, tol, file, line)
    real(real64), intent(in) :: desired(:,:,:)
    real(real64), intent(in) :: actual(:,:,:)
    real(real64), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real64) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.
    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real64)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired(i1,i2,i3) - actual(i1,i2,i3)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired(i1,i2,i3)
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3)
            write(*,'(a,es16.6)')      "Difference:    ", desired(i1,i2,i3) - actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r4_real64(desired, actual, tol, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:,:,:,:)
    real(real64), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real64) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real64)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired - actual(i1,i2,i3,i4)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*,'(a,es16.6)')      "Difference:    ", desired - actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r4_equal_r4_real64(desired, actual, tol, file, line)
    real(real64), intent(in) :: desired(:,:,:,:)
    real(real64), intent(in) :: actual(:,:,:,:)
    real(real64), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real64) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.
    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real64)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired(i1,i2,i3,i4) - actual(i1,i2,i3,i4)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired(i1,i2,i3,i4)
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*,'(a,es16.6)')      "Difference:    ", desired(i1,i2,i3,i4) - actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r5_real64(desired, actual, tol, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:,:,:,:,:)
    real(real64), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real64) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real64)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired - actual(i1,i2,i3,i4,i5)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*,'(a,es16.6)')      "Difference:    ", desired - actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r5_equal_r5_real64(desired, actual, tol, file, line)
    real(real64), intent(in) :: desired(:,:,:,:,:)
    real(real64), intent(in) :: actual(:,:,:,:,:)
    real(real64), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real64) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.
    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real64)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired(i1,i2,i3,i4,i5) - actual(i1,i2,i3,i4,i5)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired(i1,i2,i3,i4,i5)
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*,'(a,es16.6)')      "Difference:    ", desired(i1,i2,i3,i4,i5) - actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r6_real64(desired, actual, tol, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:,:,:,:,:,:)
    real(real64), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real64) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real64)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired - actual(i1,i2,i3,i4,i5,i6)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*,'(a,es16.6)')      "Difference:    ", desired - actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r6_equal_r6_real64(desired, actual, tol, file, line)
    real(real64), intent(in) :: desired(:,:,:,:,:,:)
    real(real64), intent(in) :: actual(:,:,:,:,:,:)
    real(real64), intent(in), optional :: tol !< Optional user-specified floating point tolerance
    real(real64) :: eps !< Floating point epsilon
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.
    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (present(tol)) then
      eps = tol
    else
      eps = epsilon(1.0_real64)
    end if

    if (all(abs(desired - actual) < eps)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(abs(desired(i1,i2,i3,i4,i5,i6) - actual(i1,i2,i3,i4,i5,i6)) > eps) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')          "Image:         ", this_image()
            write(*,'(a,es16.6)')      "Desired Value: ", desired(i1,i2,i3,i4,i5,i6)
            write(*,'(a,es16.6)')      "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*,'(a,es16.6)')      "Difference:    ", desired(i1,i2,i3,i4,i5,i6) - actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            write(*,'(a,es16.6)')      "Epsilon:       ", eps
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  

! -------------------------------------------------------
! Integer (=) assertions
! -------------------------------------------------------
  subroutine assert_r0_equal_r0_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual
    character(len=*), intent(in) :: file
    integer, intent(in) :: line

    logical :: correct
    correct = .false.

    if(desired == actual) correct = .true.

    if(.not. correct) then
      write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
      write(*,'(a,i0)') "Image: ", this_image()
      write(*,'(a,i0)') "Desired Value: ", desired
      write(*,'(a,i0)') "Actual Value:  ", actual
      write(*,'(a,i0)') "Difference:    ", desired - actual
      error stop "Test failure"
    end if
  end subroutine 

  subroutine assert_r0_equal_r1_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired /= actual(i1)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1)
            write(*,'(a,i0)')       "Difference:    ", desired - actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if
  end subroutine 

  subroutine assert_r1_equal_r1_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:)
    integer(int32), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if
    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1) /= actual(i1)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired(i1)
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1)
            write(*,'(a,i0)')       "Difference:    ", desired(i1) - actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if
  end subroutine 

  subroutine assert_r0_equal_r2_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired /= actual(i1,i2)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2)
            write(*,'(a,i0)')       "Difference:    ", desired - actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r2_equal_r2_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:,:)
    integer(int32), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if
    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1,i2) /= actual(i1,i2)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired(i1,i2)
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2)
            write(*,'(a,i0)')       "Difference:    ", desired(i1,i2) - actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r0_equal_r3_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired /= actual(i1,i2,i3)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3)
            write(*,'(a,i0)')       "Difference:    ", desired - actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r3_equal_r3_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:,:,:)
    integer(int32), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if
    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1,i2,i3) /= actual(i1,i2,i3)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired(i1,i2,i3)
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3)
            write(*,'(a,i0)')       "Difference:    ", desired(i1,i2,i3) - actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r0_equal_r4_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired /= actual(i1,i2,i3,i4)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*,'(a,i0)')       "Difference:    ", desired - actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r4_equal_r4_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:,:,:,:)
    integer(int32), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if
    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1,i2,i3,i4) /= actual(i1,i2,i3,i4)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired(i1,i2,i3,i4)
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*,'(a,i0)')       "Difference:    ", desired(i1,i2,i3,i4) - actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r0_equal_r5_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired /= actual(i1,i2,i3,i4,i5)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*,'(a,i0)')       "Difference:    ", desired - actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r5_equal_r5_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:,:,:,:,:)
    integer(int32), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if
    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1,i2,i3,i4,i5) /= actual(i1,i2,i3,i4,i5)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired(i1,i2,i3,i4,i5)
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*,'(a,i0)')       "Difference:    ", desired(i1,i2,i3,i4,i5) - actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r0_equal_r6_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5,i6 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired /= actual(i1,i2,i3,i4,i5,i6)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*,'(a,i0)')       "Difference:    ", desired - actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r6_equal_r6_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:,:,:,:,:,:)
    integer(int32), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5,i6 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if
    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1,i2,i3,i4,i5,i6) /= actual(i1,i2,i3,i4,i5,i6)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired(i1,i2,i3,i4,i5,i6)
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*,'(a,i0)')       "Difference:    ", desired(i1,i2,i3,i4,i5,i6) - actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r0_equal_r0_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual
    character(len=*), intent(in) :: file
    integer, intent(in) :: line

    logical :: correct
    correct = .false.

    if(desired == actual) correct = .true.

    if(.not. correct) then
      write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
      write(*,'(a,i0)') "Image: ", this_image()
      write(*,'(a,i0)') "Desired Value: ", desired
      write(*,'(a,i0)') "Actual Value:  ", actual
      write(*,'(a,i0)') "Difference:    ", desired - actual
      error stop "Test failure"
    end if
  end subroutine 

  subroutine assert_r0_equal_r1_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired /= actual(i1)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1)
            write(*,'(a,i0)')       "Difference:    ", desired - actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if
  end subroutine 

  subroutine assert_r1_equal_r1_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:)
    integer(int64), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if
    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1) /= actual(i1)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired(i1)
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1)
            write(*,'(a,i0)')       "Difference:    ", desired(i1) - actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if
  end subroutine 

  subroutine assert_r0_equal_r2_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired /= actual(i1,i2)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2)
            write(*,'(a,i0)')       "Difference:    ", desired - actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r2_equal_r2_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:,:)
    integer(int64), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if
    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1,i2) /= actual(i1,i2)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired(i1,i2)
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2)
            write(*,'(a,i0)')       "Difference:    ", desired(i1,i2) - actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r0_equal_r3_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired /= actual(i1,i2,i3)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3)
            write(*,'(a,i0)')       "Difference:    ", desired - actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r3_equal_r3_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:,:,:)
    integer(int64), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if
    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1,i2,i3) /= actual(i1,i2,i3)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired(i1,i2,i3)
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3)
            write(*,'(a,i0)')       "Difference:    ", desired(i1,i2,i3) - actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r0_equal_r4_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired /= actual(i1,i2,i3,i4)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*,'(a,i0)')       "Difference:    ", desired - actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r4_equal_r4_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:,:,:,:)
    integer(int64), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if
    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1,i2,i3,i4) /= actual(i1,i2,i3,i4)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired(i1,i2,i3,i4)
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*,'(a,i0)')       "Difference:    ", desired(i1,i2,i3,i4) - actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r0_equal_r5_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired /= actual(i1,i2,i3,i4,i5)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*,'(a,i0)')       "Difference:    ", desired - actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r5_equal_r5_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:,:,:,:,:)
    integer(int64), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if
    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1,i2,i3,i4,i5) /= actual(i1,i2,i3,i4,i5)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired(i1,i2,i3,i4,i5)
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*,'(a,i0)')       "Difference:    ", desired(i1,i2,i3,i4,i5) - actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r0_equal_r6_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5,i6 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired /= actual(i1,i2,i3,i4,i5,i6)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*,'(a,i0)')       "Difference:    ", desired - actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 

  subroutine assert_r6_equal_r6_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:,:,:,:,:,:)
    integer(int64), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5,i6 !< loop indices to check assertions elementwise

    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if
    if (all(desired == actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1,i2,i3,i4,i5,i6) /= actual(i1,i2,i3,i4,i5,i6)) then
            write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')       "Image:         ", this_image()
            write(*,'(a,i0)')       "Desired Value: ", desired(i1,i2,i3,i4,i5,i6)
            write(*,'(a,i0)')       "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*,'(a,i0)')       "Difference:    ", desired(i1,i2,i3,i4,i5,i6) - actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 


! -------------------------------------------------------
! Logical assertions
! -------------------------------------------------------
  subroutine assert_r0_equal_r0_logical(desired, actual, file, line)
    logical, intent(in) :: desired
    logical, intent(in) :: actual
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    correct = .false.

    if(desired .eqv. actual) correct = .true.

    if(.not. correct) then
      write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
      write(*,'(a,i0)') "Image: ", this_image()
      write(*,'(a,l2)') "Desired Value: ", desired
      write(*,'(a,l2)') "Actual Value:  ", actual
      error stop "Test failure"
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r1_logical(desired, actual, file, line)
    logical, intent(in) :: desired
    logical, intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1 !< loop indices to check assertions elementwise
    logical :: correct
    correct = .false.


    if (all(desired .eqv. actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired .neqv. actual(i1)) then
          write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')     "Image:         ", this_image()
            write(*,'(a,l2)')      "Desired Value: ", desired
            write(*,'(a,l2)')      "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if
  end subroutine 
  
  subroutine assert_r1_equal_r1_logical(desired, actual, file, line)
    logical, intent(in) :: desired(:)
    logical, intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1 !< loop indices to check assertions elementwise
    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired .eqv. actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1) .neqv. actual(i1)) then
          write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')     "Image:         ", this_image()
            write(*,'(a,l2)')      "Desired Value: ", desired(i1)
            write(*,'(a,l2)')      "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r2_logical(desired, actual, file, line)
    logical, intent(in) :: desired
    logical, intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2 !< loop indices to check assertions elementwise
    logical :: correct
    correct = .false.


    if (all(desired .eqv. actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired .neqv. actual(i1,i2)) then
          write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')     "Image:         ", this_image()
            write(*,'(a,l2)')      "Desired Value: ", desired
            write(*,'(a,l2)')      "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r2_equal_r2_logical(desired, actual, file, line)
    logical, intent(in) :: desired(:,:)
    logical, intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2 !< loop indices to check assertions elementwise
    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired .eqv. actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1,i2) .neqv. actual(i1,i2)) then
          write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')     "Image:         ", this_image()
            write(*,'(a,l2)')      "Desired Value: ", desired(i1,i2)
            write(*,'(a,l2)')      "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r3_logical(desired, actual, file, line)
    logical, intent(in) :: desired
    logical, intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3 !< loop indices to check assertions elementwise
    logical :: correct
    correct = .false.


    if (all(desired .eqv. actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired .neqv. actual(i1,i2,i3)) then
          write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')     "Image:         ", this_image()
            write(*,'(a,l2)')      "Desired Value: ", desired
            write(*,'(a,l2)')      "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r3_equal_r3_logical(desired, actual, file, line)
    logical, intent(in) :: desired(:,:,:)
    logical, intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3 !< loop indices to check assertions elementwise
    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired .eqv. actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1,i2,i3) .neqv. actual(i1,i2,i3)) then
          write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')     "Image:         ", this_image()
            write(*,'(a,l2)')      "Desired Value: ", desired(i1,i2,i3)
            write(*,'(a,l2)')      "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r4_logical(desired, actual, file, line)
    logical, intent(in) :: desired
    logical, intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4 !< loop indices to check assertions elementwise
    logical :: correct
    correct = .false.


    if (all(desired .eqv. actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired .neqv. actual(i1,i2,i3,i4)) then
          write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')     "Image:         ", this_image()
            write(*,'(a,l2)')      "Desired Value: ", desired
            write(*,'(a,l2)')      "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r4_equal_r4_logical(desired, actual, file, line)
    logical, intent(in) :: desired(:,:,:,:)
    logical, intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4 !< loop indices to check assertions elementwise
    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired .eqv. actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1,i2,i3,i4) .neqv. actual(i1,i2,i3,i4)) then
          write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')     "Image:         ", this_image()
            write(*,'(a,l2)')      "Desired Value: ", desired(i1,i2,i3,i4)
            write(*,'(a,l2)')      "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r5_logical(desired, actual, file, line)
    logical, intent(in) :: desired
    logical, intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5 !< loop indices to check assertions elementwise
    logical :: correct
    correct = .false.


    if (all(desired .eqv. actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired .neqv. actual(i1,i2,i3,i4,i5)) then
          write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')     "Image:         ", this_image()
            write(*,'(a,l2)')      "Desired Value: ", desired
            write(*,'(a,l2)')      "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r5_equal_r5_logical(desired, actual, file, line)
    logical, intent(in) :: desired(:,:,:,:,:)
    logical, intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5 !< loop indices to check assertions elementwise
    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired .eqv. actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1,i2,i3,i4,i5) .neqv. actual(i1,i2,i3,i4,i5)) then
          write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')     "Image:         ", this_image()
            write(*,'(a,l2)')      "Desired Value: ", desired(i1,i2,i3,i4,i5)
            write(*,'(a,l2)')      "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r0_equal_r6_logical(desired, actual, file, line)
    logical, intent(in) :: desired
    logical, intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5,i6 !< loop indices to check assertions elementwise
    logical :: correct
    correct = .false.


    if (all(desired .eqv. actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired .neqv. actual(i1,i2,i3,i4,i5,i6)) then
          write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')     "Image:         ", this_image()
            write(*,'(a,l2)')      "Desired Value: ", desired
            write(*,'(a,l2)')      "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
  subroutine assert_r6_equal_r6_logical(desired, actual, file, line)
    logical, intent(in) :: desired(:,:,:,:,:,:)
    logical, intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    integer :: i1,i2,i3,i4,i5,i6 !< loop indices to check assertions elementwise
    logical :: correct
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired .eqv. actual)) correct = .true.

    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(desired(i1,i2,i3,i4,i5,i6) .neqv. actual(i1,i2,i3,i4,i5,i6)) then
          write(*,'((3(a), i0))') "Assertion Failure (=) in ", trim(file), ":", line
            write(*,'(a,i0)')     "Image:         ", this_image()
            write(*,'(a,l2)')      "Desired Value: ", desired(i1,i2,i3,i4,i5,i6)
            write(*,'(a,l2)')      "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if
  end subroutine 
  
! -------------------------------------------------------
! < and > assertions for Reals and Integers
! -------------------------------------------------------
  subroutine assert_r0_less_r0_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    correct = .false.


    if(desired < actual) correct = .true.

    if(.not. correct) then
      write(*,'((3(a), i0))') "Assertion Failure (<) in ", trim(file), ":", line
      write(*,'(a,i0)') "Image: ", this_image()
      write(*,'(a,i0)') "Desired Value: ", desired
      write(*,'(a,i0)') "Actual Value:  ", actual
      error stop "Test failure"
    end if

    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r1_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r1_less_r1_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired(:)
    real(real32), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1) < actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1)
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r2_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r2_less_r2_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired(:,:)
    real(real32), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2) < actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2)
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r3_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r3_less_r3_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired(:,:,:)
    real(real32), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3) < actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r4_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r4_less_r4_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired(:,:,:,:)
    real(real32), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4) < actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r5_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r5_less_r5_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired(:,:,:,:,:)
    real(real32), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5) < actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r6_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r6_less_r6_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired(:,:,:,:,:,:)
    real(real32), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5,i6) < actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5,i6)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r0_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    correct = .false.


    if(desired < actual) correct = .true.

    if(.not. correct) then
      write(*,'((3(a), i0))') "Assertion Failure (<) in ", trim(file), ":", line
      write(*,'(a,i0)') "Image: ", this_image()
      write(*,'(a,i0)') "Desired Value: ", desired
      write(*,'(a,i0)') "Actual Value:  ", actual
      error stop "Test failure"
    end if

    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r1_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r1_less_r1_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired(:)
    real(real64), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1) < actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1)
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r2_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r2_less_r2_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired(:,:)
    real(real64), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2) < actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2)
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r3_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r3_less_r3_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired(:,:,:)
    real(real64), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3) < actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r4_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r4_less_r4_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired(:,:,:,:)
    real(real64), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4) < actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r5_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r5_less_r5_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired(:,:,:,:,:)
    real(real64), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5) < actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r6_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r6_less_r6_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired(:,:,:,:,:,:)
    real(real64), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5,i6) < actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5,i6)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r0_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    correct = .false.


    if(desired < actual) correct = .true.

    if(.not. correct) then
      write(*,'((3(a), i0))') "Assertion Failure (<) in ", trim(file), ":", line
      write(*,'(a,i0)') "Image: ", this_image()
      write(*,'(a,i0)') "Desired Value: ", desired
      write(*,'(a,i0)') "Actual Value:  ", actual
      error stop "Test failure"
    end if

    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r1_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r1_less_r1_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:)
    integer(int32), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1) < actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1)
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r2_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r2_less_r2_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:,:)
    integer(int32), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2) < actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2)
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r3_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r3_less_r3_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:,:,:)
    integer(int32), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3) < actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r4_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r4_less_r4_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:,:,:,:)
    integer(int32), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4) < actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r5_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r5_less_r5_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:,:,:,:,:)
    integer(int32), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5) < actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r6_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r6_less_r6_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:,:,:,:,:,:)
    integer(int32), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5,i6) < actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5,i6)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r0_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    correct = .false.


    if(desired < actual) correct = .true.

    if(.not. correct) then
      write(*,'((3(a), i0))') "Assertion Failure (<) in ", trim(file), ":", line
      write(*,'(a,i0)') "Image: ", this_image()
      write(*,'(a,i0)') "Desired Value: ", desired
      write(*,'(a,i0)') "Actual Value:  ", actual
      error stop "Test failure"
    end if

    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r1_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r1_less_r1_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:)
    integer(int64), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1) < actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1)
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r2_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r2_less_r2_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:,:)
    integer(int64), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2) < actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2)
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r3_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r3_less_r3_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:,:,:)
    integer(int64), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3) < actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r4_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r4_less_r4_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:,:,:,:)
    integer(int64), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4) < actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r5_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r5_less_r5_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:,:,:,:,:)
    integer(int64), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5) < actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_less_r6_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired < actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r6_less_r6_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:,:,:,:,:,:)
    integer(int64), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired < actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5,i6) < actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (<) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5,i6)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_less (<) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r0_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    correct = .false.


    if(desired > actual) correct = .true.

    if(.not. correct) then
      write(*,'((3(a), i0))') "Assertion Failure (>) in ", trim(file), ":", line
      write(*,'(a,i0)') "Image: ", this_image()
      write(*,'(a,i0)') "Desired Value: ", desired
      write(*,'(a,i0)') "Actual Value:  ", actual
      error stop "Test failure"
    end if

    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r1_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r1_greater_r1_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired(:)
    real(real32), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1) > actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1)
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r2_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r2_greater_r2_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired(:,:)
    real(real32), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2) > actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2)
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r3_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r3_greater_r3_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired(:,:,:)
    real(real32), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3) > actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r4_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r4_greater_r4_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired(:,:,:,:)
    real(real32), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4) > actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r5_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r5_greater_r5_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired(:,:,:,:,:)
    real(real32), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5) > actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r6_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired
    real(real32), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r6_greater_r6_real32(desired, actual, file, line)
    real(real32), intent(in) :: desired(:,:,:,:,:,:)
    real(real32), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5,i6) > actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5,i6)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r0_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    correct = .false.


    if(desired > actual) correct = .true.

    if(.not. correct) then
      write(*,'((3(a), i0))') "Assertion Failure (>) in ", trim(file), ":", line
      write(*,'(a,i0)') "Image: ", this_image()
      write(*,'(a,i0)') "Desired Value: ", desired
      write(*,'(a,i0)') "Actual Value:  ", actual
      error stop "Test failure"
    end if

    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r1_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r1_greater_r1_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired(:)
    real(real64), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1) > actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1)
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r2_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r2_greater_r2_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired(:,:)
    real(real64), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2) > actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2)
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r3_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r3_greater_r3_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired(:,:,:)
    real(real64), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3) > actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r4_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r4_greater_r4_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired(:,:,:,:)
    real(real64), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4) > actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r5_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r5_greater_r5_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired(:,:,:,:,:)
    real(real64), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5) > actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r6_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired
    real(real64), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r6_greater_r6_real64(desired, actual, file, line)
    real(real64), intent(in) :: desired(:,:,:,:,:,:)
    real(real64), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5,i6) > actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5,i6)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r0_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    correct = .false.


    if(desired > actual) correct = .true.

    if(.not. correct) then
      write(*,'((3(a), i0))') "Assertion Failure (>) in ", trim(file), ":", line
      write(*,'(a,i0)') "Image: ", this_image()
      write(*,'(a,i0)') "Desired Value: ", desired
      write(*,'(a,i0)') "Actual Value:  ", actual
      error stop "Test failure"
    end if

    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r1_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r1_greater_r1_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:)
    integer(int32), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1) > actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1)
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r2_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r2_greater_r2_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:,:)
    integer(int32), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2) > actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2)
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r3_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r3_greater_r3_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:,:,:)
    integer(int32), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3) > actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r4_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r4_greater_r4_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:,:,:,:)
    integer(int32), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4) > actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r5_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r5_greater_r5_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:,:,:,:,:)
    integer(int32), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5) > actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r6_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired
    integer(int32), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r6_greater_r6_int32(desired, actual, file, line)
    integer(int32), intent(in) :: desired(:,:,:,:,:,:)
    integer(int32), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5,i6) > actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5,i6)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r0_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    correct = .false.


    if(desired > actual) correct = .true.

    if(.not. correct) then
      write(*,'((3(a), i0))') "Assertion Failure (>) in ", trim(file), ":", line
      write(*,'(a,i0)') "Image: ", this_image()
      write(*,'(a,i0)') "Desired Value: ", desired
      write(*,'(a,i0)') "Actual Value:  ", actual
      error stop "Test failure"
    end if

    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r1_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r1_greater_r1_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:)
    integer(int64), intent(in) :: actual(:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1) > actual(i1))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1)
            write(*,*)                 "Actual Value:  ", actual(i1)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1
            error stop "Test failure"
        end if
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r2_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r2_greater_r2_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:,:)
    integer(int64), intent(in) :: actual(:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2) > actual(i1,i2))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2)
            write(*,*)                 "Actual Value:  ", actual(i1,i2)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2
            error stop "Test failure"
        end if
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r3_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r3_greater_r3_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:,:,:)
    integer(int64), intent(in) :: actual(:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3) > actual(i1,i2,i3))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3
            error stop "Test failure"
        end if
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r4_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r4_greater_r4_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:,:,:,:)
    integer(int64), intent(in) :: actual(:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4) > actual(i1,i2,i3,i4))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r5_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r5_greater_r5_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:,:,:,:,:)
    integer(int64), intent(in) :: actual(:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5) > actual(i1,i2,i3,i4,i5))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r0_greater_r6_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired
    integer(int64), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.


    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired > actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

  subroutine assert_r6_greater_r6_int64(desired, actual, file, line)
    integer(int64), intent(in) :: desired(:,:,:,:,:,:)
    integer(int64), intent(in) :: actual(:,:,:,:,:,:)
    character(len=*), intent(in) :: file
    integer, intent(in) :: line
    logical :: correct
    integer :: i1,i2,i3,i4,i5,i6  !< loop indices to check assertions elementwise
    correct = .false.

    if (size(desired) /= size(actual)) then
      write(*,'((3(a), i0))') "Assertion Failure (=): size(desired) /= size(actual) in ", trim(file), ":", line
      error stop "Test failure"
    end if

    if (all(desired > actual)) correct = .true.
    ! If all() is not correct, then loop through and figure out where
    if(.not. correct) then
      do i6 = lbound(actual, dim=6), ubound(actual, dim=6)
      do i5 = lbound(actual, dim=5), ubound(actual, dim=5)
      do i4 = lbound(actual, dim=4), ubound(actual, dim=4)
      do i3 = lbound(actual, dim=3), ubound(actual, dim=3)
      do i2 = lbound(actual, dim=2), ubound(actual, dim=2)
      do i1 = lbound(actual, dim=1), ubound(actual, dim=1)
        if(.not. (desired(i1,i2,i3,i4,i5,i6) > actual(i1,i2,i3,i4,i5,i6))) then
            write(*,'((3(a), i0))')    "Assertion Failure (>) in ", trim(file), ":", line
            write(*,*)                 "Image:         ", this_image()
            write(*,*)                 "Desired Value: ", desired(i1,i2,i3,i4,i5,i6)
            write(*,*)                 "Actual Value:  ", actual(i1,i2,i3,i4,i5,i6)
            write(*, '(a,10(i0, 1x))') "At index:      ", i1,i2,i3,i4,i5,i6
            error stop "Test failure"
        end if
      end do
      end do
      end do
      end do
      end do
      end do
    end if


    if(.not. correct) then
      write(*,*) "assert_greater (>) on image ", this_image(), " failed. desired= ", desired, " actual= ", actual
      error stop "Test failure"
    end if

  end subroutine

end module caf_testing