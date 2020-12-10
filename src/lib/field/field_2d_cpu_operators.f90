! MIT License
! Copyright (c) 2020 Sam Miller
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

#ifdef __OPENMP_THREADS__
#define  
#else
#define  pure
#endif

submodule(mod_field) field_2d_cpu_operators
  !< Summary: This is the submodule to the field_2d_t type that contains
  !<          all of the CPU arithmetic operators

contains

! --------------------------------------------------------------------
! Reduction operators (min/max/sum val and loc)
! --------------------------------------------------------------------
 module function field_maxval_cpu(f) result(res)
  class(field_2d_t), intent(in) :: f
  real(rk) :: res

  !$omp parallel workshare
  res = maxval(f%data(f%ilo:f%ihi, f%jlo:f%jhi))
  !$omp end parallel workshare
end function field_maxval_cpu

 module function field_maxloc_cpu(f) result(res)
  class(field_2d_t), intent(in) :: f
  integer(ik), dimension(2) :: res

  !$omp parallel workshare
  res = maxloc(f%data(f%ilo:f%ihi, f%jlo:f%jhi))
  !$omp end parallel workshare
end function field_maxloc_cpu

 module function field_minval_cpu(f) result(res)
  class(field_2d_t), intent(in) :: f
  real(rk) :: res

  !$omp parallel workshare
  res = minval(f%data(f%ilo:f%ihi, f%jlo:f%jhi))
  !$omp end parallel workshare
end function field_minval_cpu

 module function field_minloc_cpu(f) result(res)
  class(field_2d_t), intent(in) :: f
  integer(ik), dimension(2) :: res

  !$omp parallel workshare
  res = minloc(f%data(f%ilo:f%ihi, f%jlo:f%jhi))
  !$omp end parallel workshare
end function field_minloc_cpu

 module function field_sum_cpu(f) result(res)
  class(field_2d_t), intent(in) :: f
  real(rk) :: res

  !$omp parallel workshare
  res = sum(f%data(f%ilo:f%ihi, f%jlo:f%jhi))
  !$omp end parallel workshare
end function field_sum_cpu

! --------------------------------------------------------------------
! Arithmetic operators
! --------------------------------------------------------------------
 module subroutine field_add_field_cpu(lhs, f, res)
  !< Implementation of the field_2d_t + field_2d_t operation
  class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
  class(field_2d_t), intent(in) :: f
  type(field_2d_t), intent(inout) :: res

  integer(ik) :: i, j

  !$omp parallel default(none), &
  !$omp private(i, j) shared(lhs, res, f)
  !$omp do
  do j = lhs%jlo, lhs%jhi
    !$omp simd
    do i = lhs%ilo, lhs%ihi
      res%data(i, j) = lhs%data(i, j) + f%data(i, j)
    end do
    !$omp end simd
  end do
  !$omp end do
  !$omp end parallel
end subroutine field_add_field_cpu

 module subroutine field_sub_field_cpu(lhs, f, res)
  !< Implementation of the field_2d_t + field_2d_t operation
  class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
  class(field_2d_t), intent(in) :: f
  type(field_2d_t), intent(inout) :: res

  integer(ik) :: i, j

  !$omp parallel default(none), &
  !$omp private(i, j) shared(lhs, res, f)
  !$omp do
  do j = lhs%jlo, lhs%jhi
    !$omp simd
    do i = lhs%ilo, lhs%ihi
      res%data(i, j) = lhs%data(i, j) - f%data(i, j)
    end do
    !$omp end simd
  end do
  !$omp end do
  !$omp end parallel
end subroutine field_sub_field_cpu

 module subroutine field_mul_field_cpu(lhs, f, res)
  !< Implementation of the field_2d_t * field_2d_t operation
  class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
  class(field_2d_t), intent(in) :: f
  type(field_2d_t), intent(inout) :: res

  integer(ik) :: i, j

  !$omp parallel default(none), &
  !$omp private(i, j) shared(lhs, res, f)
  !$omp do
  do j = lhs%jlo, lhs%jhi
    !$omp simd
    do i = lhs%ilo, lhs%ihi
      res%data(i, j) = lhs%data(i, j) * f%data(i, j)
    end do
    !$omp end simd
  end do
  !$omp end do
  !$omp end parallel
end subroutine field_mul_field_cpu

 module subroutine field_div_field_cpu(lhs, f, res)
  !< Implementation of the field_2d_t * field_2d_t operation
  class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
  class(field_2d_t), intent(in) :: f
  type(field_2d_t), intent(inout) :: res

  integer(ik) :: i, j

  !$omp parallel default(none), &
  !$omp private(i, j) shared(lhs, res, f)
  !$omp do
  do j = lhs%jlo, lhs%jhi
    !$omp simd
    do i = lhs%ilo, lhs%ihi
      res%data(i, j) = lhs%data(i, j) / f%data(i, j)
    end do
    !$omp end simd
  end do
  !$omp end do
  !$omp end parallel
end subroutine field_div_field_cpu

 module subroutine field_add_real_1d_cpu(lhs, x, res)
  !< Implementation of the field_2d_t + real64 operation
  class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
  real(rk), intent(in) :: x
  type(field_2d_t), intent(inout) :: res

  integer(ik) :: i, j

  !$omp parallel default(none), &
  !$omp private(i, j) shared(lhs, res, x)
  !$omp do
  do j = lhs%jlo, lhs%jhi
    !$omp simd
    do i = lhs%ilo, lhs%ihi
      res%data(i, j) = lhs%data(i, j) + x
    end do
    !$omp end simd
  end do
  !$omp end do
  !$omp end parallel
end subroutine field_add_real_1d_cpu

 module subroutine field_add_real_2d_cpu(lhs, x, res)
  !< Implementation of the field_2d_t + real64 operation
  class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
  real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
  type(field_2d_t), intent(inout) :: res

  integer(ik) :: i, j

  !$omp parallel default(none), &
  !$omp private(i, j) shared(lhs, res, x)
  !$omp do
  do j = lhs%jlo, lhs%jhi
    !$omp simd
    do i = lhs%ilo, lhs%ihi
      res%data(i, j) = lhs%data(i, j) + x(i, j)
    end do
    !$omp end simd
  end do
  !$omp end do
  !$omp end parallel

end subroutine field_add_real_2d_cpu

 module subroutine field_sub_real_1d_cpu(lhs, x, res)
  !< Implementation of the field_2d_t + real64 operation
  class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
  real(rk), intent(in) :: x
  type(field_2d_t), intent(inout) :: res

  integer(ik) :: i, j

  !$omp parallel default(none), &
  !$omp private(i, j) shared(lhs, res, x)
  !$omp do
  do j = lhs%jlo, lhs%jhi
    !$omp simd
    do i = lhs%ilo, lhs%ihi
      res%data(i, j) = lhs%data(i, j) - x
    end do
    !$omp end simd
  end do
  !$omp end do
  !$omp end parallel

end subroutine field_sub_real_1d_cpu

 module subroutine field_sub_real_2d_cpu(lhs, x, res)
  !< Implementation of the field_2d_t + real64 operation
  class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
  real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
  type(field_2d_t), intent(inout) :: res

  integer(ik) :: i, j

  !$omp parallel default(none), &
  !$omp private(i, j) shared(lhs, res, x)
  !$omp do
  do j = lhs%jlo, lhs%jhi
    !$omp simd
    do i = lhs%ilo, lhs%ihi
      res%data(i, j) = lhs%data(i, j) - x(i, j)
    end do
    !$omp end simd
  end do
  !$omp end do
  !$omp end parallel

end subroutine field_sub_real_2d_cpu

 module subroutine real_1d_sub_field_cpu(x, rhs, res)
  !< Implementation of the field_2d_t + real64 operation
  class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
  real(rk), intent(in) :: x
  type(field_2d_t), intent(inout) :: res

  res%data = rhs%data - x
end subroutine real_1d_sub_field_cpu

 module subroutine real_2d_sub_field_cpu(x, rhs, res)
  !< Implementation of the field_2d_t + real64 operation
  class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
  real(rk), dimension(rhs%ilo:, rhs%jlo:), intent(in) :: x
  type(field_2d_t), intent(inout) :: res

  integer(ik) :: i, j

  !$omp parallel default(none), &
  !$omp private(i, j) shared(rhs, res, x)
  !$omp do
  do j = rhs%jlo, rhs%jhi
    !$omp simd
    do i = rhs%ilo, rhs%ihi
      res%data(i, j) = rhs%data(i, j) - x(i, j)
    end do
    !$omp end simd
  end do
  !$omp end do
  !$omp end parallel
end subroutine real_2d_sub_field_cpu

 module subroutine field_div_real_1d_cpu(lhs, x, res)
  !< Implementation of the field_2d_t / real64 operation
  class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
  real(rk), intent(in) :: x
  type(field_2d_t), intent(inout) :: res

  res%data = lhs%data / x
end subroutine field_div_real_1d_cpu

 module subroutine field_div_real_2d_cpu(lhs, x, res)
  !< Implementation of the field_2d_t / real64 operation
  class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
  real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
  type(field_2d_t), intent(inout) :: res

  integer(ik) :: i, j

  !$omp parallel default(none), &
  !$omp private(i, j) shared(lhs, res, x)
  !$omp do
  do j = lhs%jlo, lhs%jhi
    !$omp simd
    do i = lhs%ilo, lhs%ihi
      res%data(i, j) = lhs%data(i, j) / x(i, j)
    end do
    !$omp end simd
  end do
  !$omp end do
  !$omp end parallel
end subroutine field_div_real_2d_cpu

 module subroutine real_1d_div_field_cpu(x, rhs, res)
  !< Implementation of the field_2d_t / real64 operation
  class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
  real(rk), intent(in) :: x
  type(field_2d_t), intent(inout) :: res

  integer(ik) :: i, j

  !$omp parallel default(none), &
  !$omp private(i, j) shared(rhs, res, x)
  !$omp do
  do j = rhs%jlo, rhs%jhi
    !$omp simd
    do i = rhs%ilo, rhs%ihi
      res%data(i, j) = x / rhs%data(i, j)
    end do
    !$omp end simd
  end do
  !$omp end do
  !$omp end parallel

end subroutine real_1d_div_field_cpu

 module subroutine real_2d_div_field_cpu(x, rhs, res)
  !< Implementation of the field_2d_t / real64 operation
  class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
  real(rk), dimension(rhs%ilo:, rhs%jlo:), intent(in) :: x
  type(field_2d_t), intent(inout) :: res

  integer(ik) :: i, j

  !$omp parallel default(none), &
  !$omp private(i, j) shared(rhs, res, x)
  !$omp do
  do j = rhs%jlo, rhs%jhi
    !$omp simd
    do i = rhs%ilo, rhs%ihi
      res%data(i, j) = x(i, j) / rhs%data(i, j)
    end do
    !$omp end simd
  end do
  !$omp end do
  !$omp end parallel
end subroutine real_2d_div_field_cpu

 module subroutine field_mul_real_2d_cpu(lhs, x, res)
  !< Implementation of the field_2d_t * array operation
  class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
  real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
  type(field_2d_t), intent(inout) :: res

  integer(ik) :: i, j

  !$omp parallel default(none), &
  !$omp private(i, j) shared(lhs, res, x)
  !$omp do
  do j = lhs%jlo, lhs%jhi
    !$omp simd
    do i = lhs%ilo, lhs%ihi
      res%data(i, j) = lhs%data(i, j) * x(i, j)
    end do
    !$omp end simd
  end do
  !$omp end do
  !$omp end parallel
end subroutine field_mul_real_2d_cpu

 module subroutine field_mul_real_1d_cpu(lhs, x, res)
  !< Implementation of the field_2d_t * real64 operation
  class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
  real(rk), intent(in) :: x
  type(field_2d_t), intent(inout) :: res

  integer(ik) :: i, j

  !$omp parallel default(none), &
  !$omp private(i, j) shared(lhs, res, x)
  !$omp do
  do j = lhs%jlo, lhs%jhi
    !$omp simd
    do i = lhs%ilo, lhs%ihi
      res%data(i, j) = lhs%data(i, j) * x
    end do
    !$omp end simd
  end do
  !$omp end do
  !$omp end parallel
end subroutine field_mul_real_1d_cpu
end submodule
