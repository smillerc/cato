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

#ifdef __SIMD_ALIGN_OMP__
#define __EDGE_ALIGN__ aligned(q, edge_values, phi_r_L_i, phi_r_R_i, phi_r_L_j, phi_r_R_j, delta_i_minus, delta_i_plus, delta_j_minus, delta_j_plus)
#else
#define __EDGE_ALIGN__
#endif

module mod_tvd_2nd_order
  !< Summary: Provide class for 2nd order TVD edge interpolation
  !< Date: 06/08/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !<   [1] M. Berger, M. Aftosmis, S. Muman, "Analysis of Slope Limiters on Irregular Grids",
  !<       43rd AIAA Aerospace Sciences Meeting and Exhibit (2005), https://doi.org/10.2514/6.2005-490
  !<
  !<   [2] K.H. Kim, C. Kim, "Accurate, efficient and monotonic numerical methods for multi-dimensional compressible flows Part II: Multi-dimensional limiting process",
  !<       Journal of Computational Physics 208 (2005) 570–615, https://doi.org/10.1016/j.jcp.2005.02.022

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_flux_limiter, only: flux_limiter_t, smoothness, delta, van_leer, superbee, minmod
  use mod_edge_interpolator, only: edge_iterpolator_t
  use mod_globals, only: n_ghost_layers, debug_print

  implicit none
  private
  public :: tvd_2nd_order_t, new_tvd_2nd_order_t

  type, extends(edge_iterpolator_t) :: tvd_2nd_order_t
    !< 2nd order edge interpolation with TVD filtering
    type(flux_limiter_t) :: limiter
  contains
    procedure, public :: initialize
    procedure, public :: interpolate_edge_values
  end type tvd_2nd_order_t

contains

  function new_tvd_2nd_order_t(limiter) result(interpolator)
    type(tvd_2nd_order_t), pointer :: interpolator
    character(len=*), intent(in) :: limiter

    allocate(interpolator)
    interpolator%limiter_name = trim(limiter)
    interpolator%order = 2
    interpolator%limiter = flux_limiter_t(trim(limiter))
  end function

  subroutine initialize(self, limiter)
    class(tvd_2nd_order_t), intent(inout) :: self
    character(len=*), intent(in) :: limiter
    self%limiter_name = trim(limiter)
    self%order = 2
    self%limiter = flux_limiter_t(trim(limiter))
  end subroutine initialize

  subroutine interpolate_edge_values(self, q, lbounds, edge_values)
    !< Reconstruct the cell interface values, e.g. q_i-1/2, q_i+1/2. This assumes a cartesian
    !< structured square grid

    class(tvd_2nd_order_t), intent(in) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
    !< (i,j); primitive variable to reconstruct at the edge
    real(rk), dimension(:, :, :), allocatable, intent(out) :: edge_values
    !<((bottom, right, top, left), i, j); reconstructed edge values

    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: ilo_real, ihi_real, jlo_real, jhi_real
    real(rk), dimension(:, :), allocatable :: r_L_i  !< r_L,i in Ref[1]
    real(rk), dimension(:, :), allocatable :: r_R_i  !< r_R,i in Ref[1]
    real(rk), dimension(:, :), allocatable :: r_L_j  !< r_L,j in Ref[1]
    real(rk), dimension(:, :), allocatable :: r_R_j  !< r_R,j in Ref[1]

    real(rk), dimension(:, :), allocatable :: delta_i_plus   !<(i,j);  difference operator -> q(i+1,j) - q(i,j)
    real(rk), dimension(:, :), allocatable :: delta_i_minus  !<(i,j);  difference operator -> q(i,j) - q(i-1,j)
    real(rk), dimension(:, :), allocatable :: delta_j_plus   !<(i,j);  difference operator -> q(i,j+1) - q(i,j)
    real(rk), dimension(:, :), allocatable :: delta_j_minus  !<(i,j);  difference operator -> q(i,j) - q(i,j-1)

    real(rk), dimension(:, :), allocatable  :: phi_r_L_j !< limiter for the top edge, see Eq. 32 in Ref [1]
    real(rk), dimension(:, :), allocatable  :: phi_r_R_j !< limiter for the bottom edge, see Eq. 32 in Ref [1]
    real(rk), dimension(:, :), allocatable  :: phi_r_R_i !< limiter for the left edge, see Eq. 32 in Ref [1]
    real(rk), dimension(:, :), allocatable  :: phi_r_L_i !< limiter for the right edge, see Eq. 32 in Ref [1]

    call debug_print('Running tvd_2nd_order_t%interpolate_edge_values()', __FILE__, __LINE__)

    ilo = lbound(q, dim=1)
    ihi = ubound(q, dim=1)
    jlo = lbound(q, dim=2)
    jhi = ubound(q, dim=2)

    ! Index limits for the real domain
    ilo_real = ilo + n_ghost_layers !< first real i cell
    ihi_real = ihi - n_ghost_layers !< last  real i cell
    jlo_real = jlo + n_ghost_layers !< first real j cell
    jhi_real = jhi - n_ghost_layers !< last  real j cell

    allocate(edge_values(4, ilo:ihi, jlo:jhi))
    edge_values = 0.0_rk

    allocate(phi_r_L_i(ilo:ihi, jlo:jhi)) !dir$ attributes align:__ALIGNBYTES__ :: phi_r_L_i
    allocate(phi_r_R_i(ilo:ihi, jlo:jhi)) !dir$ attributes align:__ALIGNBYTES__ :: phi_r_R_i
    allocate(phi_r_L_j(ilo:ihi, jlo:jhi)) !dir$ attributes align:__ALIGNBYTES__ :: phi_r_L_j
    allocate(phi_r_R_j(ilo:ihi, jlo:jhi)) !dir$ attributes align:__ALIGNBYTES__ :: phi_r_R_j
    phi_r_L_i = 0.0_rk
    phi_r_R_j = 0.0_rk
    phi_r_L_j = 0.0_rk
    phi_r_R_i = 0.0_rk

    call self%get_solution_smoothness(q, lbounds=lbound(q), &
                                      delta_i_plus=delta_i_plus, delta_i_minus=delta_i_minus, &
                                      delta_j_plus=delta_j_plus, delta_j_minus=delta_j_minus, &
                                      r_L_i=r_L_i, r_R_i=r_R_i, r_L_j=r_L_j, r_R_j=r_R_j)

    call self%limit(R=r_L_i, psi=phi_r_L_i, name=self%limiter_name)
    call self%limit(R=r_R_i, psi=phi_r_R_i, name=self%limiter_name)
    call self%limit(R=r_L_j, psi=phi_r_L_j, name=self%limiter_name)
    call self%limit(R=r_R_j, psi=phi_r_R_j, name=self%limiter_name)

    !$omp parallel default(none), &
    !$omp firstprivate(ilo_real, ihi_real, jlo_real, jhi_real) &
    !$omp private(i, j) &
    !$omp shared(delta_i_plus, delta_i_minus, delta_j_plus, delta_j_minus) &
    !$omp shared(phi_r_L_i, phi_r_R_i, phi_r_L_j, phi_r_R_j) &
    !$omp shared(q, self, edge_values)

    !$omp do
    do j = jlo_real, jhi_real
      !$omp simd __EDGE_ALIGN__
      !dir$ vector aligned
      do i = ilo_real, ihi_real
        ! (i+1/2, j), cell "right" edge -> corresponds to the "L" side of the interface, thus the "L" terms
        ! phi_L = q + .5 * psi(r_L) * delta-
        edge_values(2, i, j) = q(i, j) + 0.5_rk * phi_r_L_i(i, j) * delta_i_minus(i, j)

        ! (i-1/2, j), cell "left" edge -> corresponds to the "R" side of the interface, thus the "R" terms
        edge_values(4, i, j) = q(i, j) - 0.5_rk * phi_r_R_i(i, j) * delta_i_plus(i, j)

        ! (i, j+1/2), cell "top" edge -> corresponds to the "L" side of the interface, thus the "L" terms
        edge_values(3, i, j) = q(i, j) + 0.5_rk * phi_r_L_j(i, j) * delta_j_minus(i, j)

        ! (i, j-1/2), cell "bottom" edge -> corresponds to the "R" side of the interface, thus the "R" terms
        edge_values(1, i, j) = q(i, j) - 0.5_rk * phi_r_R_j(i, j) * delta_j_plus(i, j)

      end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(r_L_i)
    deallocate(r_R_i)
    deallocate(r_L_j)
    deallocate(r_R_j)

    deallocate(phi_r_L_i)
    deallocate(phi_r_R_i)
    deallocate(phi_r_L_j)
    deallocate(phi_r_R_j)

    deallocate(delta_i_plus)
    deallocate(delta_i_minus)
    deallocate(delta_j_plus)
    deallocate(delta_j_minus)
  end subroutine interpolate_edge_values

end module mod_tvd_2nd_order
