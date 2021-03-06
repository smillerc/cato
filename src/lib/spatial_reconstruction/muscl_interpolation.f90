module mod_muscl_interpolation
  !< Summary: Provide the implementation and base class for standard MUSCL edge interpolation
  !< Date: 08/03/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !      [1]

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use, intrinsic :: ieee_arithmetic
  use mod_field, only: field_2d_t
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
    procedure(interpolate_edge_values), deferred, public :: interpolate_edge_values
  endtype

  abstract interface
    subroutine initialize(self, limiter)
      import :: muscl_interpolation_t
      class(muscl_interpolation_t), intent(inout) :: self
      character(len=*), intent(in) :: limiter
    endsubroutine initialize

    subroutine interpolate_edge_values(self, q, i_edges, j_edges)
      import :: muscl_interpolation_t, field_2d_t, ik, rk
      class(muscl_interpolation_t), intent(in) :: self
      class(field_2d_t), intent(in) :: q !< (i,j); primitive variable to reconstruct at the edge
      real(rk), dimension(:, :, :), allocatable, intent(out) :: i_edges !<((L,R), i, j); L/R state for each edge
      real(rk), dimension(:, :, :), allocatable, intent(out) :: j_edges !<((L,R), i, j); L/R state for each edge
    endsubroutine interpolate_edge_values

    subroutine distinguish_continuous_regions(self, rho, u, v, p)
      import :: muscl_interpolation_t, field_2d_t, ik, rk
      class(muscl_interpolation_t), intent(inout) :: self
      class(field_2d_t), intent(in) :: rho !< density
      class(field_2d_t), intent(in) :: u   !< x-velocity
      class(field_2d_t), intent(in) :: v   !< y-velocity
      class(field_2d_t), intent(in) :: p   !< pressure
    endsubroutine distinguish_continuous_regions
  endinterface

contains
endmodule mod_muscl_interpolation
