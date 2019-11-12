module mod_conserved_vars
  use iso_fortran_env, only: ik => int32, rk => real64

  implicit none

  private
  public :: conserved_vars_t

  type :: conserved_vars_t
    integer(ik), dimension(2) :: shape
    real(rk), allocatable, dimension(:, :) :: density
    real(rk), allocatable, dimension(:, :) :: x_velocity
    real(rk), allocatable, dimension(:, :) :: y_velocity
    real(rk), allocatable, dimension(:, :) :: pressure
  contains
    procedure :: add_conserved_vars
    procedure :: multiply_conserved_vars
    procedure :: assign_conserved_vars
    final :: finalize_conserved_vars

    generic :: operator(+) => add_conserved_vars
    generic :: operator(*) => multiply_conserved_vars
    generic :: assignment(=) => assign_conserved_vars
  end type conserved_vars_t

  interface conserved_vars_t
    module procedure :: constructor
  end interface

contains

  pure function constructor(i, j) result(U)
    !< Implementation of the conserved_vars_t, aka U, vector type

    integer(ik), intent(in) :: i, j !< Dimensions of the type
    type(conserved_vars_t), allocatable :: U

    integer(ik) :: alloc_status

    allocate(conserved_vars_t :: U)

    allocate(U%density(i, j), stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate conserved_vars_t%density"
    U%density = 0.0_rk

    allocate(U%x_velocity(i, j), stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate conserved_vars_t%x_velocity"
    U%x_velocity = 0.0_rk

    allocate(U%y_velocity(i, j), stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate conserved_vars_t%y_velocity"
    U%y_velocity = 0.0_rk

    allocate(U%pressure(i, j), stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate conserved_vars_t%pressure"
    U%pressure = 0.0_rk

    U%shape = [i, j]

  end function

  subroutine finalize_conserved_vars(self)
    !< Cleanup of the conserved_vars_t type

    type(conserved_vars_t), intent(inout) :: self
    integer(ik) :: alloc_status

    if(allocated(self%density)) deallocate(self%density, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate conserved_vars_t%density"

    if(allocated(self%x_velocity)) deallocate(self%x_velocity, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate conserved_vars_t%x_velocity"

    if(allocated(self%y_velocity)) deallocate(self%y_velocity, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate conserved_vars_t%y_velocity"

    if(allocated(self%pressure)) deallocate(self%pressure, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate conserved_vars_t%pressure"
  end subroutine

  function add_conserved_vars(lhs, rhs) result(sum)
    !< Implementation of the + operator for the conserved_vars_t type

    class(conserved_vars_t), intent(in) :: lhs, rhs
    class(conserved_vars_t), allocatable :: sum
    class(conserved_vars_t), allocatable :: local_sum
    ! integer(ik) :: alloc_status

    local_sum = conserved_vars_t(lhs%shape(1), lhs%shape(2))

    ! allocate (conserved_vars_t :: local_sum)
    ! allocate (local_sum%density, mold=rhs%density, stat=alloc_status)
    ! if(alloc_status /= 0) error stop "Unable to allocate local_sum%density in conserved_vars_t%add()"

    ! allocate (local_sum%x_velocity, mold=rhs%x_velocity, stat=alloc_status)
    ! if(alloc_status /= 0) error stop "Unable to allocate local_sum%x_velocity in conserved_vars_t%add()"

    ! allocate (local_sum%y_velocity, mold=rhs%y_velocity, stat=alloc_status)
    ! if(alloc_status /= 0) error stop "Unable to allocate local_sum%y_velocity in conserved_vars_t%add()"

    ! allocate (local_sum%pressure, mold=rhs%pressure, stat=alloc_status)
    ! if(alloc_status /= 0) error stop "Unable to allocate local_sum%pressure in conserved_vars_t%add()"

    local_sum%density = lhs%density + rhs%density
    local_sum%x_velocity = lhs%x_velocity + rhs%x_velocity
    local_sum%y_velocity = lhs%y_velocity + rhs%y_velocity
    local_sum%pressure = lhs%pressure + rhs%pressure

    call move_alloc(local_sum, sum)

  end function add_conserved_vars

  function multiply_conserved_vars(lhs, factor) result(product)
    !< Implementation of the * operator for the conserved_vars_t type

    class(conserved_vars_t), intent(in) :: lhs
    real(rk), intent(in) :: factor
    class(conserved_vars_t), allocatable :: product

    class(conserved_vars_t), allocatable :: local_product

    local_product = conserved_vars_t(lhs%shape(1), lhs%shape(2))

    local_product%density = lhs%density * factor
    local_product%x_velocity = lhs%x_velocity * factor
    local_product%y_velocity = lhs%y_velocity * factor
    local_product%pressure = lhs%pressure * factor

    call move_alloc(local_product, product)
  end function multiply_conserved_vars

  subroutine assign_conserved_vars(lhs, rhs)
    class(conserved_vars_t), intent(in) :: rhs
    class(conserved_vars_t), intent(inout) :: lhs

    lhs%density = rhs%density
    lhs%x_velocity = rhs%x_velocity
    lhs%y_velocity = rhs%y_velocity
    lhs%pressure = rhs%pressure

  end subroutine assign_conserved_vars

end module mod_conserved_vars
