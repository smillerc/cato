module mod_local_field
    !< Summary: Provide a local field class that handles abstract calculus operations. Because the 
    !<          global fluid has a coarray component, it cannot be returned by any functions (coarray restriction),
    !<          so a local class facilitates this on a per image basis

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use, intrinsic :: ieee_arithmetic
  use mod_globals, only: debug_print
  use mod_floating_point_utils, only: near_zero
  use mod_hermetic, only: hermetic

  implicit none

  private
  public :: local_field_t, new_local_field

  type, extends(hermetic) :: local_field_t
    real(rk), dimension(:, :, :), allocatable :: conserved_vars !< ((rho, rho*u, rho*v, rho*E), i, j); Conserved quantities
    integer(ik) :: ilo = 0 !< local min i index (cell-based w/ ghost included)
    integer(ik) :: ihi = 0 !< local max i index (cell-based w/ ghost included)
    integer(ik) :: jlo = 0 !< local min j index (cell-based w/ ghost included)
    integer(ik) :: jhi = 0 !< local max j index (cell-based w/ ghost included)
    ! logical :: on_ilo_bc = .false. !< is this image on the boundary for ilo?
    ! logical :: on_ihi_bc = .false. !< is this image on the boundary for ihi?
    ! logical :: on_jlo_bc = .false. !< is this image on the boundary for jlo?
    ! logical :: on_jhi_bc = .false. !< is this image on the boundary for jhi?
  contains
    procedure, pass(lhs), public :: add_field
    procedure, pass(lhs), public :: subtract_field
    procedure, pass(lhs), public :: field_mul_real
    procedure, pass(rhs), public :: real_mul_field
    procedure, pass(lhs), public :: assign_field
    procedure, public :: force_finalization
    procedure, private, nopass :: add_fields
    procedure, private, nopass :: subtract_fields
    procedure, private, nopass :: mult_fields
    final :: finalize

    ! Map operators to corresponding procedures
    generic :: operator(+) => add_field
    generic :: operator(-) => subtract_field
    generic :: operator(*) => field_mul_real, real_mul_field
    generic :: assignment(=) => assign_field
  end type local_field_t
contains

  function new_local_field(field_data, lbounds, ubounds) result(field)
    !< Constructor for local_field_t
    real(rk), dimension(:, :, :), intent(in) :: field_data !< ((rho, rho*u, rho*v, rho*E), i, j); field data, aka conserved variables
    integer(ik), dimension(2), intent(in) :: lbounds !< (ilo, jlo); lower bound
    integer(ik), dimension(2), intent(in) :: ubounds !< (ihi, jhi); upper bound
    class(local_field_t), allocatable :: field

    associate(ilo=>lbounds(1), ihi=>ubounds(2), &
              jlo=>lbounds(1), jhi=>ubounds(2))
      field%ilo = ilo
      field%ihi = ihi
      field%jlo = jlo
      field%jhi = jhi
      allocate(field%conserved_vars(4, ilo:ihi, jlo:jhi))
    end associate

    field%conserved_vars = field_data
  end function

  subroutine force_finalization(self)
    class(local_field_t), intent(inout) :: self

    call debug_print('Running local_field_t%force_finalization()', __FILE__, __LINE__)
    if(allocated(self%conserved_vars)) deallocate(self%conserved_vars)
  end subroutine force_finalization

  subroutine finalize(self)
    type(local_field_t), intent(inout) :: self

    call debug_print('Running local_field_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%conserved_vars)) deallocate(self%conserved_vars)
  end subroutine finalize

  function subtract_field(lhs, rhs) result(difference)
    !< Implementation of the (-) operator for the field type
    class(local_field_t), intent(in) :: lhs, rhs
    type(local_field_t), allocatable :: difference
    type(local_field_t), allocatable :: local_difference

    call debug_print('Running local_field_t%subtract_field()', __FILE__, __LINE__)

    select type(rhs)
    class is(local_field_t)
      allocate(local_difference, source=lhs)
      call add_fields(a=lhs%conserved_vars, b=rhs%conserved_vars, c=local_difference%conserved_vars) ! c=a-b
    class default
      error stop 'local_field_t%subtract_field: unsupported rhs class'
    end select

    call move_alloc(local_difference, difference)
    call difference%set_temp(calling_function='subtract_field (difference)', line=__LINE__)
  end function subtract_field

  function add_field(lhs, rhs) result(sum)
    !< Implementation of the (+) operator for the field type
    class(local_field_t), intent(in) :: lhs
    class(local_field_t), intent(in) :: rhs

    class(local_field_t), allocatable :: sum
    type(local_field_t), allocatable :: local_sum

    call debug_print('Running local_field_t%add_field()', __FILE__, __LINE__)

    select type(rhs)
    class is(local_field_t)
      allocate(local_sum, source=lhs)
      ! local_sum%conserved_vars = lhs%conserved_vars + rhs%conserved_vars
      call add_fields(a=lhs%conserved_vars, b=rhs%conserved_vars, c=local_sum%conserved_vars) ! c=a+b
    class default
      error stop 'local_field_t%add_field: unsupported rhs class'
    end select

    call local_sum%set_temp(calling_function='local_field_t%add_field(local_sum)', line=__LINE__)
    call move_alloc(local_sum, sum)
    call sum%set_temp(calling_function='local_field_t%add_field(sum)', line=__LINE__)
  end function add_field

  pure subroutine add_fields(a, b, c)
    !< Dumb routine for a vectorized version of c = a + b
    real(rk), dimension(:, :, :), intent(in), contiguous :: a
    real(rk), dimension(:, :, :), intent(in), contiguous :: b
    real(rk), dimension(:, :, :), intent(inout), contiguous :: c
    integer(ik) :: i, j, k

    do k = lbound(a, dim=3), ubound(a, dim=3)
      do j = lbound(a, dim=2), ubound(a, dim=2)
        do i = lbound(a, dim=1), ubound(a, dim=1)
          c(i, j, k) = a(i, j, k) + b(i, j, k)
        end do
      end do
    end do
  end subroutine add_fields

  pure subroutine subtract_fields(a, b, c)
    !< Dumb routine for a vectorized version of c = a - b
    real(rk), dimension(:, :, :), intent(in), contiguous :: a
    real(rk), dimension(:, :, :), intent(in), contiguous :: b
    real(rk), dimension(:, :, :), intent(inout), contiguous :: c
    integer(ik) :: i, j, k

    do k = lbound(a, dim=3), ubound(a, dim=3)
      do j = lbound(a, dim=2), ubound(a, dim=2)
        do i = lbound(a, dim=1), ubound(a, dim=1)
          c(i, j, k) = a(i, j, k) - b(i, j, k)
        end do
      end do
    end do
  end subroutine subtract_fields

  pure subroutine mult_fields(a, b, c)
    !< Dumb routine for a vectorized version of c = a * b
    real(rk), dimension(:, :, :), intent(in), contiguous :: a
    real(rk), dimension(:, :, :), intent(in), contiguous :: b
    real(rk), dimension(:, :, :), intent(inout), contiguous :: c
    integer(ik) :: i, j, k

    do k = lbound(a, dim=3), ubound(a, dim=3)
      do j = lbound(a, dim=2), ubound(a, dim=2)
        do i = lbound(a, dim=1), ubound(a, dim=1)
          c(i, j, k) = a(i, j, k) * b(i, j, k)
        end do
      end do
    end do
  end subroutine mult_fields

  pure subroutine mult_field_by_real(a, b, c)
    !< Dumb routine for a vectorized version of c = a * b
    real(rk), dimension(:, :, :), intent(in), contiguous :: a
    real(rk), intent(in) :: b
    real(rk), dimension(:, :, :), intent(inout), contiguous :: c
    integer(ik) :: i, j, k

    do k = lbound(a, dim=3), ubound(a, dim=3)
      do j = lbound(a, dim=2), ubound(a, dim=2)
        do i = lbound(a, dim=1), ubound(a, dim=1)
          c(i, j, k) = a(i, j, k) * b
        end do
      end do
    end do
  end subroutine mult_field_by_real

  function field_mul_real(lhs, rhs) result(product)
    !< Implementation of the field * real operation
    class(local_field_t), intent(in) :: lhs
    real(rk), intent(in) :: rhs
    class(local_field_t), allocatable :: product

    class(local_field_t), allocatable :: local_product
    integer(ik) :: alloc_status

    call debug_print('Running local_field_t%field_mul_real()', __FILE__, __LINE__)

    allocate(local_product, source=lhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate local_product in local_field_t%field_mul_real"

    ! local_product%conserved_vars = lhs%conserved_vars * rhs
    call mult_field_by_real(a=lhs%conserved_vars, b=rhs, c=local_product%conserved_vars)

    call move_alloc(local_product, product)
    call product%set_temp(calling_function='field_mul_real (product)', line=__LINE__)
  end function field_mul_real

  function real_mul_field(lhs, rhs) result(product)
    !< Implementation of the real * field operation
    class(local_field_t), intent(in) :: rhs
    real(rk), intent(in) :: lhs
    class(local_field_t), allocatable :: product

    type(local_field_t), allocatable :: local_product
    integer(ik) :: alloc_status

    call debug_print('Running local_field_t%real_mul_field()', __FILE__, __LINE__)

    allocate(local_product, source=rhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate local_product in local_field_t%real_mul_field"

    call mult_field_by_real(a=rhs%conserved_vars, b=lhs, c=local_product%conserved_vars)

    call move_alloc(local_product, product)
    call product%set_temp(calling_function='real_mul_field (product)', line=__LINE__)
  end function real_mul_field

  subroutine assign_field(lhs, rhs)
    !< Implementation of the (=) operator for the field type. e.g. lhs = rhs
    class(local_field_t), intent(inout) :: lhs
    class(local_field_t), intent(in) :: rhs
    integer(ik) :: alloc_status

    alloc_status = 0
    call debug_print('Running local_field_t%assign_field', __FILE__, __LINE__)

    call rhs%guard_temp(calling_function='assign_field (rhs)', line=__LINE__)
    select type(rhs)
    class is(local_field_t)
      lhs%conserved_vars = rhs%conserved_vars
    class default
      error stop 'Error in local_field_t%assign_field: unsupported class'
    end select

    call rhs%clean_temp(calling_function='assign_field (rhs)', line=__LINE__)
  end subroutine assign_field
end module mod_local_field
