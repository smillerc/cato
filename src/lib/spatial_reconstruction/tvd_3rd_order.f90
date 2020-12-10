! MIT License
! Copyright (c) 2019 Sam Miller
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module mod_tvd_3rd_order
  !< Summary: Provide class for 3rd order TVD edge interpolation
  !< Date: 06/08/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !<   [1] K.H. Kim, C. Kim, "Accurate, efficient and monotonic numerical methods for multi-dimensional compressible flows Part II: Multi-dimensional limiting process",
  !<       Journal of Computational Physics 208 (2005) 570–615, https://doi.org/10.1016/j.jcp.2005.02.022

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_flux_limiter, only: flux_limiter_t, smoothness, delta
  use mod_edge_interpolator, only: edge_iterpolator_t
  use mod_globals, only: n_ghost_layers, debug_print

  implicit none
  private
  public :: tvd_3rd_order_t, new_tvd_3rd_order_t

  type, extends(edge_iterpolator_t) :: tvd_3rd_order_t
    !< 3rd order edge interpolation with TVD filtering
    type(flux_limiter_t) :: limiter
  contains
    procedure, public :: initialize
    ! procedure, nopass, public :: get_beta => beta_3rd_order
    procedure, public :: interpolate_edge_values
  endtype tvd_3rd_order_t

contains

  function new_tvd_3rd_order_t(limiter) result(interpolator)
    type(tvd_3rd_order_t), pointer :: interpolator
    character(len=*), intent(in) :: limiter

    allocate(interpolator)
    interpolator%limiter_name = trim(limiter)
    interpolator%order = 3
    interpolator%limiter = flux_limiter_t(trim(limiter))
  endfunction

  subroutine initialize(self, limiter)
    !< Constructor for tvd_3rd_order_t
    class(tvd_3rd_order_t), intent(inout) :: self
    character(len=*), intent(in) :: limiter
    self%limiter_name = trim(limiter)
    self%order = 3
    self%limiter = flux_limiter_t(trim(limiter))
  endsubroutine initialize

  subroutine interpolate_edge_values(self, q, lbounds, edge_values)
    !< Reconstruct the edge values
    class(tvd_3rd_order_t), intent(in) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
    !< (i,j); primitive variable to reconstruct at the edge
    real(rk), dimension(:, :, :), allocatable, intent(out) :: edge_values
    !<((bottom, right, top, left), i, j); reconstructed edge values

    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: ilo_bc, ihi_bc, jlo_bc, jhi_bc
    real(rk), dimension(:, :), allocatable :: r_L_i  !< r_L,i in Ref[1]
    real(rk), dimension(:, :), allocatable :: r_R_i  !< r_R,i in Ref[1]
    real(rk), dimension(:, :), allocatable :: r_L_j  !< r_L,j in Ref[1]
    real(rk), dimension(:, :), allocatable :: r_R_j  !< r_R,j in Ref[1]

    real(rk), dimension(:, :), allocatable :: beta_L_i !< beta_L,i in Ref[1]
    real(rk), dimension(:, :), allocatable :: beta_R_i !< beta_R,i in Ref[1]
    real(rk), dimension(:, :), allocatable :: beta_L_j !< beta_L,j in Ref[1]
    real(rk), dimension(:, :), allocatable :: beta_R_j !< beta_R,j in Ref[1]

    real(rk) :: phi_top    !< limiter for the top edge, see Eq. 32 in Ref [1]
    real(rk) :: phi_bottom !< limiter for the bottom edge, see Eq. 32 in Ref [1]
    real(rk) :: phi_left   !< limiter for the left edge, see Eq. 32 in Ref [1]
    real(rk) :: phi_right  !< limiter for the right edge, see Eq. 32 in Ref [1]

    real(rk) :: delta_i_plus, delta_i_minus, delta_j_plus, delta_j_minus

    call debug_print('Running tvd_3rd_order_t%interpolate_edge_values()', __FILE__, __LINE__)

    ilo_bc = lbound(q, dim=1)
    ihi_bc = ubound(q, dim=1)
    jlo_bc = lbound(q, dim=2)
    jhi_bc = ubound(q, dim=2)

    ! Index limits for the real domain
    ilo = ilo_bc + n_ghost_layers
    ihi = ihi_bc - n_ghost_layers
    jlo = jlo_bc + n_ghost_layers
    jhi = jhi_bc - n_ghost_layers

    allocate(edge_values(4, ilo_bc:ihi_bc, jlo_bc:jhi_bc))

    allocate(r_L_i(ilo - 1:ihi + 1, jlo - 1:jhi + 1))
    allocate(r_R_i(ilo - 1:ihi + 1, jlo - 1:jhi + 1))
    allocate(r_L_j(ilo - 1:ihi + 1, jlo - 1:jhi + 1))
    allocate(r_R_j(ilo - 1:ihi + 1, jlo - 1:jhi + 1))

    allocate(beta_L_i(ilo - 1:ihi + 1, jlo - 1:jhi + 1))
    allocate(beta_R_i(ilo - 1:ihi + 1, jlo - 1:jhi + 1))
    allocate(beta_L_j(ilo - 1:ihi + 1, jlo - 1:jhi + 1))
    allocate(beta_R_j(ilo - 1:ihi + 1, jlo - 1:jhi + 1))

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp private(phi_bottom, phi_top, phi_left, phi_right) &
    !$omp private(delta_i_plus, delta_i_minus, delta_j_plus, delta_j_minus) &
    !$omp shared(r_L_i, r_R_i, r_L_j, r_R_j) &
    !$omp shared(beta_L_i, beta_R_i, beta_L_j, beta_R_j) &
    !$omp shared(q, self, edge_values)

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        r_L_i(i, j) = smoothness(q(i - 1, j), q(i, j), q(i + 1, j))
        r_R_i(i, j) = 1.0_rk / r_L_i(i, j)
        r_L_j(i, j) = smoothness(q(i, j - 1), q(i, j), q(i, j + 1))
        r_R_j(i, j) = 1.0_rk / r_L_j(i, j)
      enddo
    enddo
    !$omp end do
    !$omp barrier

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        ! 3rd order interpolation function
        beta_L_i(i, j) = (1.0_rk + 2.0_rk * r_L_i(i, j)) / 3.0_rk
        beta_R_i(i, j) = (1.0_rk + 2.0_rk * r_R_i(i, j)) / 3.0_rk
        beta_L_j(i, j) = (1.0_rk + 2.0_rk * r_L_j(i, j)) / 3.0_rk
        beta_R_j(i, j) = (1.0_rk + 2.0_rk * r_R_j(i, j)) / 3.0_rk
      enddo
    enddo
    !$omp end do
    !$omp barrier

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        ! (i+1/2, j), cell "right" edge -> corresponds to the "L" side of the interface, thus the "L" terms
        phi_right = max(0.0_rk, min(2.0_rk, 2.0_rk * r_L_i(i, j), beta_L_i(i, j)))
        delta_i_minus = delta(q(i, j), q(i - 1, j)) ! q(i,j) - q(i-1,j)
        edge_values(2, i, j) = q(i, j) + 0.5_rk * phi_right * delta_i_minus

        ! (i-1/2, j), cell "left" edge -> corresponds to the "R" side of the interface, thus the "R" terms
        phi_left = max(0.0_rk, min(2.0_rk, 2.0_rk * r_R_i(i, j), beta_R_i(i, j)))
        delta_i_plus = delta(q(i + 1, j), q(i, j)) ! q(i+1,j) - q(i,j)
        edge_values(4, i, j) = q(i, j) - 0.5_rk * phi_left * delta_i_plus

        ! (i, j+1/2), cell "top" edge -> corresponds to the "L" side of the interface, thus the "L" terms
        phi_top = max(0.0_rk, min(2.0_rk, 2.0_rk * r_L_j(i, j), beta_L_j(i, j)))
        delta_j_minus = delta(q(i, j), q(i, j - 1)) ! q(i,j) - q(i,j-1)
        edge_values(3, i, j) = q(i, j) + 0.5_rk * phi_top * delta_j_minus

        ! (i, j-1/2), cell "bottom" edge -> corresponds to the "R" side of the interface, thus the "R" terms
        phi_bottom = max(0.0_rk, min(2.0_rk, 2.0_rk * r_R_j(i, j), beta_R_j(i, j)))
        delta_j_plus = delta(q(i, j + 1), q(i, j)) ! q(i,j+1) - q(i,j)
        edge_values(1, i, j) = q(i, j) - 0.5_rk * phi_bottom * delta_j_plus
      enddo
    enddo
    !$omp end do
    !$omp end parallel

    deallocate(r_L_i)
    deallocate(r_R_i)
    deallocate(r_L_j)
    deallocate(r_R_j)

    deallocate(beta_L_i)
    deallocate(beta_R_i)
    deallocate(beta_L_j)
    deallocate(beta_R_j)

  endsubroutine interpolate_edge_values

endmodule mod_tvd_3rd_order
