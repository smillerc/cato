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
#define __PSI_ALIGN__ aligned(psi, R:__ALIGNBYTES__)
#define __DELTAQ_ALIGN__ aligned(delta_i_minus, delta_i_plus, delta_j_plus, delta_j_minus, q:__ALIGNBYTES__)
#define __DELTAR_ALIGN__ aligned(delta_plus, delta_minus, R, R_inv:__ALIGNBYTES__)
#else
#define __PSI_ALIGN__
#define __DELTAQ_ALIGN__
#define __DELTAR_ALIGN__
#endif

module mod_edge_interpolator
  !< Summary: Provide baseline class for edge interpolation schemes
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
  use mod_globals, only: n_ghost_layers, debug_print
  use mod_flux_limiter, only: flux_limiter_t, minmod, superbee

  implicit none

  type, abstract :: edge_iterpolator_t
    private
    integer(ik), public :: order = 2 !< interpolation order
    character(len=32), public :: name = ''
    character(len=32), public :: limiter_name = ''
  contains
    procedure(init), deferred, public ::  initialize
    procedure(basic_interface), deferred, public :: interpolate_edge_values
    procedure, public, nopass :: get_deltas
    procedure, public, nopass :: limit
    procedure, public, nopass :: get_smoothness_R
  end type edge_iterpolator_t

  abstract interface
    subroutine init(self, limiter)
      import :: edge_iterpolator_t
      class(edge_iterpolator_t), intent(inout) :: self
      character(len=*), intent(in) :: limiter
    end subroutine init

    subroutine basic_interface(self, q, lbounds, edge_values)
      import :: edge_iterpolator_t, ik, rk
      class(edge_iterpolator_t), intent(in) :: self
      integer(ik), dimension(2), intent(in) :: lbounds
      real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
      !< (i,j); primitive variable to reconstruct at the edge
      real(rk), dimension(:, :, :), allocatable, intent(out) :: edge_values
      !<((bottom, right, top, left), i, j); reconstructed edge values
    end subroutine
  end interface

contains

  subroutine get_deltas(q, q_lbounds, delta_i_plus, delta_i_minus, delta_j_plus, delta_j_minus, delta_lbounds)
    !< Get the solution smoothness at each cell. This is sent to the limiters and MUSCL interpolation. The lbound
    !< arrays are needed b/c q includes ghost regions, whereas the delta arrays do not, since they are only for
    !< non-ghost cells

    integer(ik), dimension(2), intent(in) :: q_lbounds     !< lower bounds for the q array
    integer(ik), dimension(2), intent(in) :: delta_lbounds !< lower bounds for the delta arrays; indexing is different that q
    real(rk), dimension(q_lbounds(1):, q_lbounds(2):), contiguous, intent(in) :: q
    real(rk), dimension(delta_lbounds(1):, delta_lbounds(2):), contiguous, intent(inout) :: delta_i_plus  !< (i,j); q(i+1,j) - q(i,j)
    real(rk), dimension(delta_lbounds(1):, delta_lbounds(2):), contiguous, intent(inout) :: delta_i_minus !< (i,j); q(i,j) - q(i-1,j)
    real(rk), dimension(delta_lbounds(1):, delta_lbounds(2):), contiguous, intent(inout) :: delta_j_plus  !< (i,j); q(i,j+1) - q(i,j)
    real(rk), dimension(delta_lbounds(1):, delta_lbounds(2):), contiguous, intent(inout) :: delta_j_minus !< (i,j); q(i,j) - q(i,j-1)

    ! Locals
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk), parameter :: REL_TOL = 1e-12_rk     !< relative error tolerance
    real(rk), parameter :: ABS_TOL = tiny(1.0_rk) !< absolute error tolerance
    real(rk) :: abs_err !< absolute error

    !dir$ assume_aligned q: __ALIGNBYTES__
    !dir$ assume_aligned delta_i_plus: __ALIGNBYTES__
    !dir$ assume_aligned delta_j_plus: __ALIGNBYTES__
    !dir$ assume_aligned delta_i_minus: __ALIGNBYTES__
    !dir$ assume_aligned delta_j_minus: __ALIGNBYTES__

    ilo = lbound(delta_i_plus, dim=1)  !< first real i cell
    ihi = ubound(delta_i_plus, dim=1)  !< last  real i cell
    jlo = lbound(delta_i_plus, dim=2)  !< first real j cell
    jhi = ubound(delta_i_plus, dim=2)  !< last  real j cell

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j, abs_err) &
    !$omp shared(delta_i_plus, delta_i_minus, delta_j_plus, delta_j_minus, q)
    !$omp do
    do j = jlo, jhi
      !$omp simd __DELTAQ_ALIGN__
      !dir$ vector aligned
      do i = ilo, ihi
        delta_i_plus(i, j) = q(i + 1, j) - q(i, j)
        abs_err = ABS_TOL + REL_TOL * max(abs(q(i + 1, j)), abs(q(i, j)))

        ! Error checks
        if(abs(delta_i_plus(i, j)) < epsilon(1.0_rk)) then
          delta_i_plus(i, j) = 0.0_rk
        else if(abs(q(i + 1, j)) < tiny(1.0_rk) .and. abs(q(i, j)) < tiny(1.0_rk)) then
          delta_i_plus(i, j) = 0.0_rk
        else if(abs(delta_i_plus(i, j)) < abs_err) then
          delta_i_plus(i, j) = 0.0_rk
        end if
      end do
    end do
    !$omp end do

    !$omp do
    do j = jlo, jhi
      !$omp simd __DELTAQ_ALIGN__
      !dir$ vector aligned
      do i = ilo, ihi
        delta_j_plus(i, j) = q(i, j + 1) - q(i, j)
        abs_err = ABS_TOL + REL_TOL * max(abs(q(i, j + 1)), abs(q(i, j)))

        ! Error checks
        if(abs(delta_j_plus(i, j)) < epsilon(1.0_rk)) then
          delta_j_plus(i, j) = 0.0_rk
        else if(abs(q(i, j + 1)) < tiny(1.0_rk) .and. abs(q(i, j)) < tiny(1.0_rk)) then
          delta_j_plus(i, j) = 0.0_rk
        else if(abs(delta_j_plus(i, j)) < abs_err) then
          delta_j_plus(i, j) = 0.0_rk
        end if
      end do
    end do
    !$omp end do

    !$omp do
    do j = jlo, jhi
      !$omp simd __DELTAQ_ALIGN__
      !dir$ vector aligned
      do i = ilo, ihi
        delta_i_minus(i, j) = delta_i_plus(i - 1, j)
        delta_j_minus(i, j) = delta_j_plus(i, j - 1)
      end do
    end do
    !$omp end do

    !$omp end parallel
  end subroutine get_deltas

  subroutine get_smoothness_R(delta_plus, delta_minus, R, R_inv)
    real(rk), dimension(:, :), contiguous, intent(in) :: delta_plus
    real(rk), dimension(:, :), contiguous, intent(in) :: delta_minus
    real(rk), dimension(:, :), contiguous, intent(inout) :: R     !<
    real(rk), dimension(:, :), contiguous, intent(inout) :: R_inv !< 1/R

    ! Locals
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk), parameter :: EPS = 1e-30_rk

    ! Intel compiler hints for alignment
    !dir$ assume_aligned delta_plus: __ALIGNBYTES__
    !dir$ assume_aligned delta_minus: __ALIGNBYTES__

    ilo = lbound(delta_plus, dim=1)  !< first real i cell
    ihi = ubound(delta_plus, dim=1)  !< last  real i cell
    jlo = lbound(delta_plus, dim=2)  !< first real j cell
    jhi = ubound(delta_plus, dim=2)  !< last  real j cell

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp shared(delta_plus, delta_minus, R, R_inv)
    !$omp do
    do j = jlo, jhi
      !$omp simd __DELTAR_ALIGN__

      !dir$ vector aligned
      do i = ilo, ihi
        R(i, j) = (delta_minus(i, j) + EPS) / (delta_plus(i, j) + EPS)
        R_inv(i, j) = (delta_plus(i, j) + EPS) / (delta_minus(i, j) + EPS)
      end do
    end do
    !$omp end do

    !$omp do
    do j = jlo, jhi
      !$omp simd __DELTAR_ALIGN__
      !dir$ vector aligned
      do i = ilo, ihi
        if(abs(R(i, j)) < EPS) R(i, j) = 0.0_rk
        if(abs(R_inv(i, j)) < EPS) R_inv(i, j) = 0.0_rk
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine get_smoothness_R

  subroutine limit(R, psi, name)
    !< Apply the flux limiter based on the smoothness
    real(rk), dimension(:, :), contiguous, intent(in) :: R !< smoothness
    real(rk), dimension(:, :), contiguous, intent(inout) :: psi !< flux limiter value
    character(len=*), intent(in) :: name !< name of the limiter

    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(RK), parameter :: PSI_EPS = 1e-5_rk
    !dir$ assume_aligned R: __ALIGNBYTES__
    !dir$ assume_aligned psi: __ALIGNBYTES__

    ilo = lbound(R, dim=1)
    ihi = ubound(R, dim=1)
    jlo = lbound(R, dim=2)
    jhi = ubound(R, dim=2)
    select case(trim(name))
    case('minmod')
      !$omp parallel default(none), &
      !$omp firstprivate(ilo, ihi, jlo, jhi) &
      !$omp private(i, j) shared(psi, R)
      !$omp do
      do j = jlo, jhi
        !dir$ vector aligned
        !$omp simd __PSI_ALIGN__
        do i = ilo, ihi
          psi(i, j) = max(0.0_rk, min(R(i, j), 1.0_rk))
        end do
      end do

      !$omp do
      do j = jlo, jhi
        !dir$ vector aligned
        !$omp simd __PSI_ALIGN__
        do i = ilo, ihi
          if(R(i, j) < PSI_EPS) then
            psi(i, j) = 0.0_rk
          else if(R(i, j) > 1.0_rk) then
            psi(i, j) = 1.0_rk
          endif
        end do
      end do
      !$omp end do
      !$omp end parallel

    case('superbee')
      !$omp parallel default(none), &
      !$omp firstprivate(ilo, ihi, jlo, jhi) &
      !$omp private(i, j) shared(psi, R)
      !$omp do
      do j = jlo, jhi
        !dir$ vector aligned
        !$omp simd __PSI_ALIGN__
        do i = ilo, ihi
          psi(i, j) = max(0.0_rk, min(2.0_rk * R(i, j), 1.0_rk), min(R(i, j), 2.0_rk))
        end do
      end do

      !$omp do
      do j = jlo, jhi
        !dir$ vector aligned
        !$omp simd __PSI_ALIGN__
        do i = ilo, ihi
          if(R(i, j) < PSI_EPS) then
            psi(i, j) = 0.0_rk
          else if(R(i, j) > 2.0_rk) then
            psi(i, j) = 2.0_rk
          endif
        end do
      end do
      !$omp end do
      !$omp end parallel
    case('van_leer')
      !$omp parallel default(none), &
      !$omp firstprivate(ilo, ihi, jlo, jhi) &
      !$omp private(i, j) shared(psi, R)
      !$omp do
      do j = jlo, jhi
        !dir$ vector aligned
        !$omp simd __PSI_ALIGN__
        do i = ilo, ihi
          psi(i, j) = (R(i, j) + abs(R(i, j))) / (1.0_rk + abs(R(i, j)))
        end do
      end do

      !$omp do
      do j = jlo, jhi
        !dir$ vector aligned
        !$omp simd __PSI_ALIGN__
        do i = ilo, ihi
          if(R(i, j) < PSI_EPS) psi(i, j) = 0.0_rk
        end do
      end do
      !$omp end do
      !$omp end parallel
    end select
  end subroutine limit

end module mod_edge_interpolator
