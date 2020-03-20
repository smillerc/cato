module local_fluid
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use, intrinsic :: ieee_arithmetic
  use mod_globals, only: debug_print, print_evolved_cell_data, print_recon_data
  use mod_floating_point_utils, only: near_zero
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t
  use mod_grid, only: grid_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use hdf5_interface, only: hdf5_file
  use mod_finite_volume_schemes, only: finite_volume_scheme_t
  use mod_flux_tensor, only: operator(.dot.), H => flux_tensor_t
  use mod_time_integrator_factory, only: time_integrator_factory

  implicit none

  private
  public :: local_fluid_t, new_fluid

  type :: local_fluid_t
    real(rk), dimension(:, :, :), allocatable :: conserved_vars !< ((rho, rho*u, rho*v, rho*E), i, j); Conserved quantities
  contains
    procedure, public :: initialize
    procedure, private :: initialize_from_ini
    procedure, private :: initialize_from_hdf5
    procedure, public :: t => time_derivative
    procedure, public :: get_sound_speed
    procedure, public :: get_max_sound_speed
    procedure, public :: get_primitive_vars
    procedure, public :: sanity_check
    procedure, nopass, private :: flux_edges
    procedure, pass(lhs), public :: type_plus_type => add_fluid
    procedure, pass(lhs), public :: type_minus_type => subtract_fluid
    procedure, pass(lhs), public :: type_mul_real => fluid_mul_real
    procedure, pass(rhs), public :: real_mul_type => real_mul_fluid
    procedure, pass(lhs), public :: assign => assign_fluid
    procedure, public :: force_finalization
    procedure, private, nopass :: add_fields
    final :: finalize
  end type local_fluid_t
contains

  function subtract_fluid(lhs, rhs) result(difference)
    !< Implementation of the (-) operator for the fluid type
    class(local_fluid_t), intent(in) :: lhs, rhs
    type(local_fluid_t), allocatable :: difference

    call debug_print('Running local_fluid_t%subtract_fluid()', __FILE__, __LINE__)

    select type(rhs)
    class is(local_fluid_t)
      allocate(local_difference, source=lhs)
      call add_fields(a=lhs%conserved_vars, b=rhs%conserved_vars, c=local_difference%conserved_vars) ! c=a-b
    class default
      error stop 'local_fluid_t%subtract_fluid: unsupported rhs class'
    end select

    call move_alloc(local_difference, difference)
    call difference%set_temp(calling_function='subtract_fluid (difference)', line=__LINE__)
  end function subtract_fluid

  function add_fluid(lhs, rhs) result(sum)
    !< Implementation of the (+) operator for the fluid type
    class(local_fluid_t), intent(in) :: lhs
    class(integrand_t), intent(in) :: rhs

    class(integrand_t), allocatable :: sum
    type(local_fluid_t), allocatable :: local_sum

    call debug_print('Running local_fluid_t%add_fluid()', __FILE__, __LINE__)

    select type(rhs)
    class is(local_fluid_t)
      allocate(local_sum, source=lhs)
      ! local_sum%conserved_vars = lhs%conserved_vars + rhs%conserved_vars
      call add_fields(a=lhs%conserved_vars, b=rhs%conserved_vars, c=local_sum%conserved_vars) ! c=a+b
    class default
      error stop 'local_fluid_t%add_fluid: unsupported rhs class'
    end select

    call local_sum%set_temp(calling_function='local_fluid_t%add_fluid(local_sum)', line=__LINE__)
    call move_alloc(local_sum, sum)
    call sum%set_temp(calling_function='local_fluid_t%add_fluid(sum)', line=__LINE__)
  end function add_fluid

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

  function fluid_mul_real(lhs, rhs) result(product)
    !< Implementation of the fluid * real operation
    class(local_fluid_t), intent(in) :: lhs
    real(rk), intent(in) :: rhs
    class(integrand_t), allocatable :: product

    class(local_fluid_t), allocatable :: local_product
    integer(ik) :: alloc_status

    call debug_print('Running local_fluid_t%fluid_mul_real()', __FILE__, __LINE__)

    allocate(local_product, source=lhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate local_product in local_fluid_t%fluid_mul_real"

    ! local_product%conserved_vars = lhs%conserved_vars * rhs
    call mult_field_by_real(a=lhs%conserved_vars, b=rhs, c=local_product%conserved_vars)

    call move_alloc(local_product, product)
    call product%set_temp(calling_function='fluid_mul_real (product)', line=__LINE__)
  end function fluid_mul_real

  function real_mul_fluid(lhs, rhs) result(product)
    !< Implementation of the real * fluid operation
    class(local_fluid_t), intent(in) :: rhs
    real(rk), intent(in) :: lhs
    class(integrand_t), allocatable :: product

    type(local_fluid_t), allocatable :: local_product
    integer(ik) :: alloc_status

    call debug_print('Running local_fluid_t%real_mul_fluid()', __FILE__, __LINE__)

    allocate(local_product, source=rhs, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate local_product in local_fluid_t%real_mul_fluid"

    local_product%time_integrator = rhs%time_integrator
    ! local_product%conserved_vars = rhs%conserved_vars * lhs
    call mult_field_by_real(a=rhs%conserved_vars, b=lhs, c=local_product%conserved_vars)

    call move_alloc(local_product, product)
    call product%set_temp(calling_function='real_mul_fluid (product)', line=__LINE__)
  end function real_mul_fluid

  subroutine assign_fluid(lhs, rhs)
    !< Implementation of the (=) operator for the fluid type. e.g. lhs = rhs
    class(local_fluid_t), intent(inout) :: lhs
    class(integrand_t), intent(in) :: rhs
    integer(ik) :: alloc_status

    alloc_status = 0
    call debug_print('Running local_fluid_t%assign_fluid', __FILE__, __LINE__)

    call rhs%guard_temp(calling_function='assign_fluid (rhs)', line=__LINE__)
    select type(rhs)
    class is(local_fluid_t)
      lhs%time_integrator = rhs%time_integrator
      lhs%conserved_vars = rhs%conserved_vars
    class default
      error stop 'Error in local_fluid_t%assign_fluid: unsupported class'
    end select

    call rhs%clean_temp(calling_function='assign_fluid (rhs)', line=__LINE__)
  end subroutine assign_fluid
end module
