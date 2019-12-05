module mod_boundary_conditions

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t

  implicit none

  private
  public :: boundary_condition_t

  type, abstract :: boundary_condition_t
    character(:), allocatable :: name
    character(len=2) :: location  !< Location (+x, -x, +y, or -y)
  contains
    procedure(initialize), public, deferred :: initialize
    procedure(apply_bc), public, deferred :: apply_bc
  end type boundary_condition_t

  abstract interface
    subroutine initialize(self, location)
      import :: boundary_condition_t
      ! import :: input_t
      class(boundary_condition_t), intent(inout) :: self
      ! class(input_t), intent(in) :: input
      character(len=2), intent(in) :: location !< Location (+x, -x, +y, or -y)
    end subroutine initialize

    subroutine apply_bc(self, conserved_vars, reconstructed_state)
      import :: boundary_condition_t
      import :: rk
      class(boundary_condition_t), intent(in) :: self
      real(rk), dimension(:, :, :), intent(inout) :: conserved_vars
      real(rk), dimension(:, :, :, :, :), intent(inout) :: reconstructed_state
    end subroutine apply_bc
  end interface
contains

end module mod_boundary_conditions
