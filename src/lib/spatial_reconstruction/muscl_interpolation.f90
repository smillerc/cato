module mod_muscl_interpolation
  !> Summary: Provide the implementation and base class for standard MUSCL edge interpolation
  !> Date: 08/03/2020
  !> Author: Sam Miller
  !> Notes:
  !> References:
  !      [1]

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use, intrinsic :: ieee_arithmetic
  use mod_globals, only: n_ghost_layers, debug_print
  use mod_error, only: error_msg

  implicit none
  private
  public :: muscl_interpolation_t

  type, abstract :: muscl_interpolation_t
    private
    integer(ik), public :: order = 2 !< interpolation order
    character(len=32), public :: name = ''
    character(len=32), public :: limiter_name = ''
  contains
    procedure(init), deferred, public ::  initialize
    procedure(basic_interface), deferred, public :: interpolate_edge_values
  end type

  abstract interface
    subroutine init(self, limiter)
      import :: muscl_interpolation_t
      class(muscl_interpolation_t), intent(inout) :: self
      character(len=*), intent(in) :: limiter
    end subroutine init

    subroutine basic_interface(self, q, lbounds, i_edges, j_edges)
      import :: muscl_interpolation_t, ik, rk
      class(muscl_interpolation_t), intent(in) :: self
      integer(ik), dimension(2), intent(in) :: lbounds

      real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
      !< (i,j); primitive variable to reconstruct at the edge

      real(rk), dimension(:, :, :), allocatable, intent(out) :: i_edges
      !<((L,R), i, j); L/R state for each edge
      real(rk), dimension(:, :, :), allocatable, intent(out) :: j_edges
    end subroutine
  end interface

contains

end module mod_muscl_interpolation
