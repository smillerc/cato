submodule(mod_field) field_2d_gpu_operators
contains

  ! --------------------------------------------------------------------
  ! Reduction operators (min/max/sum val and loc)
  ! --------------------------------------------------------------------
  module function field_maxval_gpu(f) result(res)
    class(field_2d_t), intent(in) :: f
    real(rk) :: res
  end function field_maxval_gpu

  module function field_maxloc_gpu(f) result(res)
    class(field_2d_t), intent(in) :: f
    integer(ik), dimension(2) :: res
  end function field_maxloc_gpu

  module function field_minval_gpu(f) result(res)
    class(field_2d_t), intent(in) :: f
    real(rk) :: res
  end function field_minval_gpu

  module function field_minloc_gpu(f) result(res)
    class(field_2d_t), intent(in) :: f
    integer(ik), dimension(2) :: res
  end function field_minloc_gpu

  module function field_sum_gpu(f) result(res)
    class(field_2d_t), intent(in) :: f
    real(rk) :: res
  end function field_sum_gpu

  ! --------------------------------------------------------------------
  ! Arithmetic operators
  ! --------------------------------------------------------------------
  module subroutine field_add_field_gpu(lhs, f, res)
    !< Implementation of the field_2d_t + field_2d_t operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    class(field_2d_t), intent(in) :: f
    type(field_2d_t), intent(inout) :: res
  end subroutine field_add_field_gpu

  module subroutine field_sub_field_gpu(lhs, f, res)
    !< Implementation of the field_2d_t + field_2d_t operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    class(field_2d_t), intent(in) :: f
    type(field_2d_t), intent(inout) :: res
  end subroutine field_sub_field_gpu

  module subroutine field_mul_field_gpu(lhs, f, res)
    !< Implementation of the field_2d_t * field_2d_t operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    class(field_2d_t), intent(in) :: f
    type(field_2d_t), intent(inout) :: res
  end subroutine field_mul_field_gpu

  module subroutine field_div_field_gpu(lhs, f, res)
    !< Implementation of the field_2d_t * field_2d_t operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    class(field_2d_t), intent(in) :: f
    type(field_2d_t), intent(inout) :: res
  end subroutine field_div_field_gpu

  module subroutine field_add_real_1d_gpu(lhs, x, res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), intent(in) :: x
    type(field_2d_t), intent(inout) :: res
  end subroutine field_add_real_1d_gpu

  module subroutine field_add_real_2d_gpu(lhs, x, res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
    type(field_2d_t), intent(inout) :: res
  end subroutine field_add_real_2d_gpu

  module subroutine field_sub_real_1d_gpu(lhs, x, res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), intent(in) :: x
    type(field_2d_t), intent(inout) :: res
  end subroutine field_sub_real_1d_gpu

  module subroutine field_sub_real_2d_gpu(lhs, x, res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
    type(field_2d_t), intent(inout) :: res
  end subroutine field_sub_real_2d_gpu

  module subroutine real_1d_sub_field_gpu(x, rhs, res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), intent(in) :: x
    type(field_2d_t), intent(inout) :: res
  end subroutine real_1d_sub_field_gpu

  module subroutine real_2d_sub_field_gpu(x, rhs, res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), dimension(rhs%ilo:, rhs%jlo:), intent(in) :: x
    type(field_2d_t), intent(inout) :: res
  end subroutine real_2d_sub_field_gpu

  module subroutine field_div_real_1d_gpu(lhs, x, res)
    !< Implementation of the field_2d_t / real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), intent(in) :: x
    type(field_2d_t), intent(inout) :: res
  end subroutine field_div_real_1d_gpu

  module subroutine field_div_real_2d_gpu(lhs, x, res)
    !< Implementation of the field_2d_t / real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
    type(field_2d_t), intent(inout) :: res
  end subroutine field_div_real_2d_gpu

  module subroutine real_1d_div_field_gpu(x, rhs, res)
    !< Implementation of the field_2d_t / real64 operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), intent(in) :: x
    type(field_2d_t), intent(inout) :: res
  end subroutine real_1d_div_field_gpu

  module subroutine real_2d_div_field_gpu(x, rhs, res)
    !< Implementation of the field_2d_t / real64 operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), dimension(rhs%ilo:, rhs%jlo:), intent(in) :: x
    type(field_2d_t), intent(inout) :: res
  end subroutine real_2d_div_field_gpu

  module subroutine field_mul_real_2d_gpu(lhs, x, res)
    !< Implementation of the field_2d_t * array operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
    type(field_2d_t), intent(inout) :: res
  end subroutine field_mul_real_2d_gpu

  module subroutine field_mul_real_1d_gpu(lhs, x, res)
    !< Implementation of the field_2d_t * real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), intent(in) :: x
    type(field_2d_t), intent(inout) :: res
  end subroutine field_mul_real_1d_gpu
end submodule
