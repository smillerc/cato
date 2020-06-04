module mod_boundary_conditions
  !< Define the based boundary condition class

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t

  implicit none

  private
  public :: boundary_condition_t

  type, abstract :: boundary_condition_t
    character(len=32) :: name
    character(len=2) :: location  !< Location (+x, -x, +y, or -y)
    real(rk), private :: time = 0.0_rk ! Solution time (for time dependent bc's)
    real(rk) :: max_time = 0.0_rk !< Max time in source (e.g. stop after this)
    integer(ik) :: priority !< Certain b.c.'s take priority over others in terms of when they are applied (periodic is last)
    integer(ik) :: io_unit
  contains
    procedure, public :: set_time
    procedure, public :: get_time
    procedure(apply_primitive_var_bc), public, deferred :: apply_primitive_var_bc
    procedure(apply_reconstructed_state_bc), public, deferred :: apply_reconstructed_state_bc
    procedure(copy_bc), public, deferred :: copy
    generic :: assignment(=) => copy

  end type boundary_condition_t

  abstract interface
    subroutine copy_bc(out_bc, in_bc)
      import :: boundary_condition_t
      class(boundary_condition_t), intent(in) :: in_bc
      class(boundary_condition_t), intent(inout) :: out_bc
    end subroutine

    subroutine apply_primitive_var_bc(self, rho, u, v, p, lbounds)
      import :: boundary_condition_t
      import :: ik, rk
      class(boundary_condition_t), intent(inout) :: self
      integer(ik), dimension(2), intent(in) :: lbounds
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: rho
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: u
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: v
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: p
    end subroutine apply_primitive_var_bc

    subroutine apply_reconstructed_state_bc(self, recon_rho, recon_u, recon_v, recon_p, lbounds)
      import :: boundary_condition_t
      import :: ik, rk
      class(boundary_condition_t), intent(inout) :: self
      integer(ik), dimension(2), intent(in) :: lbounds
      real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_rho
      real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_u
      real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_v
      real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_p
    end subroutine apply_reconstructed_state_bc

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
