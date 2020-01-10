module mod_boundary_conditions

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t

  implicit none

  private
  public :: boundary_condition_t

  type, abstract :: boundary_condition_t
    character(len=32) :: name
    character(len=2) :: location  !< Location (+x, -x, +y, or -y)
    real(rk), private :: time = 0.0_rk ! Solution time (for time dependent bc's)
    integer(ik) :: priority !< Certain b.c.'s take priority over others in terms of when they are applied (periodic is last)
  contains
    ! procedure(initialize), public, deferred :: initialize
    procedure, public :: set_time
    procedure, public :: get_time
    procedure(apply_conserved_var_bc), public, deferred :: apply_conserved_var_bc
    procedure(apply_reconstructed_state_bc), public, deferred :: apply_reconstructed_state_bc
    procedure(apply_cell_gradient_bc), public, deferred :: apply_cell_gradient_bc
    procedure(copy_bc), public, deferred :: copy
    generic :: assignment(=) => copy

  end type boundary_condition_t

  abstract interface
    ! subroutine initialize(self, location)
    !   import :: boundary_condition_t
    !   ! import :: input_t
    !   class(boundary_condition_t), intent(inout) :: self
    !   ! class(input_t), intent(in) :: input
    !   character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    ! end subroutine initialize

    subroutine copy_bc(out_bc, in_bc)
      import :: boundary_condition_t
      class(boundary_condition_t), intent(in) :: in_bc
      class(boundary_condition_t), intent(inout) :: out_bc
    end subroutine

    subroutine apply_conserved_var_bc(self, conserved_vars)
      import :: boundary_condition_t
      import :: rk
      class(boundary_condition_t), intent(inout) :: self
      real(rk), dimension(:, 0:, 0:), intent(inout) :: conserved_vars
    end subroutine apply_conserved_var_bc

    subroutine apply_reconstructed_state_bc(self, reconstructed_state)
      import :: boundary_condition_t
      import :: rk
      class(boundary_condition_t), intent(in) :: self
      real(rk), dimension(:, :, :, 0:, 0:), intent(inout) :: reconstructed_state
    end subroutine apply_reconstructed_state_bc

    subroutine apply_cell_gradient_bc(self, cell_gradient)
      import :: boundary_condition_t
      import :: rk
      class(boundary_condition_t), intent(in) :: self
      real(rk), dimension(:, :, 0:, 0:), intent(inout) :: cell_gradient
    end subroutine apply_cell_gradient_bc
  end interface
contains
  subroutine set_time(self, time)
    class(boundary_condition_t), intent(inout) :: self
    real(rk), intent(in) :: time
    self%time = time
  end subroutine

  function get_time(self) result(time)
    class(boundary_condition_t), intent(in) :: self
    real(rk) :: time
    time = self%time
  end function
end module mod_boundary_conditions
