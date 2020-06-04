module mod_mlp
  !> Summary: Module to provide the Multidimensional Limiting Process spatial reconstruction
  !> Date: 06/04/2020
  !> Author: Sam Miller
  !> Notes:
  !> References:
  !>    [1] Kyu Hong Kim, Chongam Kim, "Accurate, efficient and monotonic numerical methods for multi-dimensional compressible flows",
  !>         https://doi.org/10.1016/j.jcp.2005.02.022

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_globals, only: debug_print, plot_limiters
  use mod_flux_limiter, only: flux_limiter_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_globals, only: MACHINE_EPS, n_ghost_layers

  implicit none

  private
  public :: reconstruct_edge_values_mlp, reconstruct_edge_values_mlp3, reconstruct_edge_values_mlp5

contains
  subroutine reconstruct_edge_values_mlp(q, lbounds, limiter, edge_values)
    !< Plain-jane MLP reconstruction

    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
    !< (i,j); primitive variable to reconstruct at the edge

    real(rk), dimension(:, :, :), allocatable, intent(out) :: edge_values
    !<((bottom, right, top, left), i, j); reconstructed edge values

    type(slope_limiter_t), intent(in) :: limiter !< slope limiter used to reconstruct the edge interface

  end subroutine reconstruct_edge_values_mlp

  subroutine reconstruct_edge_values_mlp3(q, lbounds, limiter, edge_values)
    !< 3rd order MLP reconstruction

    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
    !< (i,j); primitive variable to reconstruct at the edge

    real(rk), dimension(:, :, :), allocatable, intent(out) :: edge_values
    !<((bottom, right, top, left), i, j); reconstructed edge values

    type(slope_limiter_t), intent(in) :: limiter !< slope limiter used to reconstruct the edge interface

  end subroutine reconstruct_edge_values_mlp3

  subroutine reconstruct_edge_values_mlp5(q, lbounds, limiter, edge_values)
    !< 5th order MLP reconstruction

    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
    !< (i,j); primitive variable to reconstruct at the edge

    real(rk), dimension(:, :, :), allocatable, intent(out) :: edge_values
    !<((bottom, right, top, left), i, j); reconstructed edge values

    type(slope_limiter_t), intent(in) :: limiter !< slope limiter used to reconstruct the edge interface

  end subroutine reconstruct_edge_values_mlp5
end module mod_mlp
