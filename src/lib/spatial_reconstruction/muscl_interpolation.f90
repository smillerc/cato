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
    character(len=32), public :: name = '' !< name of the interpolation scheme, e.g. 'TVD2'
    character(len=32), public :: limiter_name = '' !< Name of the TVD limiter, e.g. 'minmod'
  contains
    ! procedure(initialize), deferred, public ::  initialize
    procedure(distinguish_continuous_regions), deferred, public :: distinguish_continuous_regions
    procedure(interpolate_edge_values), deferred, public :: interpolate_i_edge_values
    procedure(interpolate_edge_values), deferred, public :: interpolate_j_edge_values
  end type

  abstract interface
    subroutine initialize(self, limiter)
      import :: muscl_interpolation_t
      class(muscl_interpolation_t), intent(inout) :: self
      character(len=*), intent(in) :: limiter
    end subroutine initialize

    subroutine interpolate_edge_values(self, q, q_lbounds, edges, e_lbounds)
      import :: muscl_interpolation_t, ik, rk
      class(muscl_interpolation_t), intent(in) :: self
      integer(ik), dimension(2), intent(in) :: q_lbounds
      integer(ik), dimension(3), intent(in) :: e_lbounds
      real(rk), dimension(q_lbounds(1):, q_lbounds(2):), codimension[*], contiguous, intent(in) :: q !< (i,j); primitive variable to reconstruct at the edge
      real(rk), dimension(e_lbounds(1):, e_lbounds(2):, e_lbounds(3):), codimension[*], intent(inout) :: edges !<((L,R), i, j); L/R state for each i/j edge
    end subroutine interpolate_edge_values

    subroutine distinguish_continuous_regions(self, rho, u, v, p, lbounds)
      import :: muscl_interpolation_t, ik, rk
      class(muscl_interpolation_t), intent(inout) :: self
      integer(ik), dimension(2), intent(in) :: lbounds
      real(rk), dimension(lbounds(1):, lbounds(2):), codimension[*], contiguous, intent(in) :: rho !< (i,j); density
      real(rk), dimension(lbounds(1):, lbounds(2):), codimension[*], contiguous, intent(in) :: u !< (i,j); x-velocity
      real(rk), dimension(lbounds(1):, lbounds(2):), codimension[*], contiguous, intent(in) :: v !< (i,j); y-velocity
      real(rk), dimension(lbounds(1):, lbounds(2):), codimension[*], contiguous, intent(in) :: p !< (i,j); pressure

    end subroutine distinguish_continuous_regions
  end interface

contains
end module mod_muscl_interpolation
