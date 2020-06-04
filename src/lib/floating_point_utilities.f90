module mod_floating_point_utils
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  private
  public :: nearly_equal, near_zero, equal, EPS, neumaier_sum, neumaier_sum_2, neumaier_sum_3, neumaier_sum_4

  real(rk), parameter :: TINY_RK = tiny(1.0_rk)
  real(rk), parameter :: EPS = epsilon(1.0_rk)

  interface near_zero
    module procedure :: near_zero_base, near_zero_custom
  end interface

  interface equal
    module procedure :: equal_base, equal_custom
  end interface

contains

  elemental function near_zero_base(number) result(return_value)
    !< Test if a floading point is near 0
    real(rk), intent(in) :: number
    logical :: return_value

    return_value = abs(number) < EPS
  end function near_zero_base

  elemental function near_zero_custom(number, epsilon) result(return_value)
    !< Test if a floading point is near 0. Use this instead of doing a
    !< `if (number == 0.0_rk)` check. Taken from "Modern Fortran: Style and Usage"
    !<  by N. Clerman and W. Spector in Chapter 13, Section 2, Rule 171
    real(rk), intent(in) :: number
    real(rk), intent(in) :: epsilon
    logical :: return_value
    real(rk) :: local_epsilon

    if(abs(epsilon) >= EPS) then
      local_epsilon = abs(epsilon)
    else
      local_epsilon = EPS
    end if
    return_value = abs(number) < local_epsilon

  end function near_zero_custom

  elemental function equal_base(lhs, rhs) result(return_value)
    !< Test if two floating point numbers are equal
    real(rk), intent(in)  :: lhs
    real(rk), intent(in)  :: rhs
    logical :: return_value

    return_value = abs(lhs - rhs) < EPS

  end function equal_base

  logical function nearly_equal(a, b) result(equal)
    real(rk), intent(in) :: a, b
    real(rk) :: diff, rel_error

    equal = .false.
    diff = abs(a - b)

    if(a == b) then ! infinities
      equal = .true.
    else if(a == 0.0_rk .or. b == 0.0_rk .or. &
            (abs(a) + abs(b)) < tiny(1.0_rk)) then
      equal = diff < (epsilon(1.0_rk) * tiny(1.0_rk))
    else
      equal = (diff / min((abs(a) + abs(b)), huge(1.0_rk))) < epsilon(1.0_rk)
      print *, 'rel error: ', diff / min((abs(a) + abs(b)), huge(1.0_rk))
    end if

  end function

  elemental function equal_custom(lhs, rhs, epsilon) result(return_value)
    !< Test if two floating point numbers are equal

    real(rk), intent(in)  :: lhs
    real(rk), intent(in)  :: rhs
    real(rk), intent(in) :: epsilon
    logical :: return_value
    real(rk) :: local_epsilon

    if(abs(epsilon) >= EPS) then
      local_epsilon = abs(epsilon)
    else
      local_epsilon = EPS
    end if
    return_value = abs(lhs - rhs) < local_epsilon

  end function equal_custom

  pure real(rk) function neumaier_sum(vector) result(sum)
    !< Implementation of an accurate floating point summation algorithm. See
    !< https://en.wikipedia.org/wiki/Kahan_summation_algorithm. This significantly reduces round off error.
    !< In order for the compiler not to aggressively optimize the benifits of this away, the ifort
    !< compiler flag 'fp-model source` needs to be included.

    real(rk), dimension(:), intent(in) :: vector

    integer(ik) :: i
    real(rk) :: c ! accumulator
    real(rk) :: t ! total

    c = 0.0_rk
    sum = 0.0_rk

    do i = 1, size(vector)
      t = sum + vector(i)
      if(abs(sum) >= abs(vector(i))) then
        c = c + ((sum - t) + vector(i)) ! If sum is bigger, low-order digits of vector(i) are lost.
      else
        c = c + ((vector(i) - t) + sum) ! Else low-order digits of sum are lost.
      endif
      sum = t
    end do
    sum = sum + c
  end function neumaier_sum

  pure real(rk) function neumaier_sum_2(vector) result(sum)
    !< Implementation of an accurate floating point summation algorithm. See
    !< https://en.wikipedia.org/wiki/Kahan_summation_algorithm. This significantly reduces round off error.
    !< In order for the compiler not to aggressively optimize the benifits of this away, the ifort
    !< compiler flag 'fp-model source` needs to be included.

    real(rk), dimension(2), intent(in) :: vector

    integer(ik) :: i
    real(rk) :: c ! accumulator
    real(rk) :: t ! total

    c = 0.0_rk
    sum = 0.0_rk

    do i = 1, 2
      t = sum + vector(i)
      if(abs(sum) >= abs(vector(i))) then
        c = c + ((sum - t) + vector(i)) ! If sum is bigger, low-order digits of vector(i) are lost.
      else
        c = c + ((vector(i) - t) + sum) ! Else low-order digits of sum are lost.
      endif
      sum = t
    end do
    sum = sum + c
  end function neumaier_sum_2

  pure real(rk) function neumaier_sum_3(vector) result(sum)
    !< Implementation of an accurate floating point summation algorithm. See
    !< https://en.wikipedia.org/wiki/Kahan_summation_algorithm. This significantly reduces round off error.
    !< In order for the compiler not to aggressively optimize the benifits of this away, the ifort
    !< compiler flag 'fp-model source` needs to be included.

    real(rk), dimension(3), intent(in) :: vector

    integer(ik) :: i
    real(rk) :: c ! accumulator
    real(rk) :: t ! total

    c = 0.0_rk
    sum = 0.0_rk

    do i = 1, 3
      t = sum + vector(i)
      if(abs(sum) >= abs(vector(i))) then
        c = c + ((sum - t) + vector(i)) ! If sum is bigger, low-order digits of vector(i) are lost.
      else
        c = c + ((vector(i) - t) + sum) ! Else low-order digits of sum are lost.
      endif
      sum = t
    end do
    sum = sum + c
  end function neumaier_sum_3

  pure real(rk) function neumaier_sum_4(vector) result(sum)
    !< Implementation of an accurate floating point summation algorithm. See
    !< https://en.wikipedia.org/wiki/Kahan_summation_algorithm. This significantly reduces round off error.
    !< In order for the compiler not to aggressively optimize the benifits of this away, the ifort
    !< compiler flag 'fp-model source` needs to be included.

    real(rk), dimension(4), intent(in) :: vector

    integer(ik) :: i
    real(rk) :: c ! accumulator
    real(rk) :: t ! total

    c = 0.0_rk
    sum = 0.0_rk

    do i = 1, 4
      t = sum + vector(i)
      if(abs(sum) >= abs(vector(i))) then
        c = c + ((sum - t) + vector(i)) ! If sum is bigger, low-order digits of vector(i) are lost.
      else
        c = c + ((vector(i) - t) + sum) ! Else low-order digits of sum are lost.
      endif
      sum = t
    end do
    sum = sum + c
  end function neumaier_sum_4

end module mod_floating_point_utils
