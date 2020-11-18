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
#define __PSI_ALIGN__ aligned(psi, r:__ALIGNBYTES__)
#define __DELTAQ_ALIGN__ aligned(delta_i_minus, delta_i_plus, delta_j_plus, delta_j_minus, q:__ALIGNBYTES__)
#define __DELTAR_ALIGN__ aligned(delta_plus, delta_minus, r, r_inv:__ALIGNBYTES__)
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

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
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
    procedure, public, nopass :: get_r_smoothness
    procedure, public, nopass :: get_solution_smoothness
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

  subroutine get_solution_smoothness(q, lbounds, delta_i_plus, delta_i_minus, delta_j_plus, delta_j_minus, &
                                     r_L_i, r_R_i, r_L_j, r_R_j)
    !< Find the smoothness of the solution based on the primitive variable q. This smoothness is dependent on the
    !< nearest neighbor

    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q !<(i,j); primitive variable

    real(rk), dimension(:, :), allocatable, intent(out) :: r_L_i  !<(i,j); r_L,i in Ref[1]
    real(rk), dimension(:, :), allocatable, intent(out) :: r_R_i  !<(i,j); r_R,i in Ref[1]
    real(rk), dimension(:, :), allocatable, intent(out) :: r_L_j  !<(i,j); r_L,j in Ref[1]
    real(rk), dimension(:, :), allocatable, intent(out) :: r_R_j  !<(i,j); r_R,j in Ref[1]

    real(rk), dimension(:, :), allocatable, intent(out) :: delta_i_plus   !<(i,j);  difference operator -> q(i+1,j) - q(i,j)
    real(rk), dimension(:, :), allocatable, intent(out) :: delta_i_minus  !<(i,j);  difference operator -> q(i,j) - q(i-1,j)
    real(rk), dimension(:, :), allocatable, intent(out) :: delta_j_plus   !<(i,j);  difference operator -> q(i,j+1) - q(i,j)
    real(rk), dimension(:, :), allocatable, intent(out) :: delta_j_minus  !<(i,j);  difference operator -> q(i,j) - q(i,j-1)

    integer(ik) :: ilo, ihi, jlo, jhi

    call debug_print('Running edge_iterpolator_t%get_solution_smoothness()', __FILE__, __LINE__)

    ilo = lbound(q, dim=1)
    ihi = ubound(q, dim=1)
    jlo = lbound(q, dim=2)
    jhi = ubound(q, dim=2)

    allocate(delta_i_plus(ilo:ihi, jlo:jhi))  !dir$ attributes align:__ALIGNBYTES__ :: delta_i_plus
    allocate(delta_i_minus(ilo:ihi, jlo:jhi)) !dir$ attributes align:__ALIGNBYTES__ :: delta_i_minus
    allocate(delta_j_plus(ilo:ihi, jlo:jhi))  !dir$ attributes align:__ALIGNBYTES__ :: delta_j_plus
    allocate(delta_j_minus(ilo:ihi, jlo:jhi)) !dir$ attributes align:__ALIGNBYTES__ :: delta_j_minus
    delta_i_plus = 0.0_rk
    delta_i_minus = 0.0_rk
    delta_j_plus = 0.0_rk
    delta_j_minus = 0.0_rk

    allocate(r_L_i(ilo:ihi, jlo:jhi)) !dir$ attributes align:__ALIGNBYTES__ :: r_L_i
    allocate(r_L_j(ilo:ihi, jlo:jhi)) !dir$ attributes align:__ALIGNBYTES__ :: r_L_j
    allocate(r_R_i(ilo:ihi, jlo:jhi)) !dir$ attributes align:__ALIGNBYTES__ :: r_R_i
    allocate(r_R_j(ilo:ihi, jlo:jhi)) !dir$ attributes align:__ALIGNBYTES__ :: r_R_j
    r_L_i = 0.0_rk
    r_R_i = 0.0_rk
    r_L_j = 0.0_rk
    r_R_j = 0.0_rk

    call get_deltas(q, &
                    delta_i_plus=delta_i_plus, &
                    delta_i_minus=delta_i_minus, &
                    delta_j_plus=delta_j_plus, &
                    delta_j_minus=delta_j_minus, &
                    lbounds=lbound(q))

    call get_r_smoothness(delta_plus=delta_i_plus, delta_minus=delta_i_minus, R=r_L_i, R_inv=r_R_i)
    call get_r_smoothness(delta_plus=delta_j_plus, delta_minus=delta_j_minus, R=r_L_j, R_inv=r_R_j)

  end subroutine get_solution_smoothness

  subroutine get_deltas(q, delta_i_plus, delta_i_minus, delta_j_plus, delta_j_minus, lbounds)
    !< Get the solution smoothness at each cell. This is sent to the limiters and MUSCL interpolation. The lbound
    !< arrays are needed b/c q includes ghost regions, whereas the delta arrays do not, since they are only for
    !< non-ghost cells

    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(inout) :: delta_i_plus  !< (i,j); q(i+1,j) - q(i,j)
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(inout) :: delta_i_minus !< (i,j); q(i,j) - q(i-1,j)
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(inout) :: delta_j_plus  !< (i,j); q(i,j+1) - q(i,j)
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(inout) :: delta_j_minus !< (i,j); q(i,j) - q(i,j-1)

    ! Locals
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk), parameter :: REL_TOL = 1e-12_rk     !< relative error tolerance
    real(rk), parameter :: ABS_TOL = tiny(1.0_rk) !< absolute error tolerance
    real(rk) :: abs_err !< absolute error

    !dir$ assume_aligned delta_i_plus: __ALIGNBYTES__
    !dir$ assume_aligned delta_j_plus: __ALIGNBYTES__
    !dir$ assume_aligned delta_i_minus: __ALIGNBYTES__
    !dir$ assume_aligned delta_j_minus: __ALIGNBYTES__

    ! The delta arrays have the same size and indexing as the q array does. Each cell
    ! has a delta+ and a delta- value. The very edges of the q array will not have a
    ! delta b/c they are on the edge of the total domain. With multiple ghost layers,
    ! the outermost ghost cells will not have a delta, but the inner ones will
    ilo = lbound(q, dim=1) + n_ghost_layers  !< first real i cell
    ihi = ubound(q, dim=1) - n_ghost_layers  !< last  real i cell
    jlo = lbound(q, dim=2) + n_ghost_layers  !< first real j cell
    jhi = ubound(q, dim=2) - n_ghost_layers  !< last  real j cell

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j, abs_err) &
    !$omp shared(delta_i_plus, delta_i_minus, delta_j_plus, delta_j_minus, q)

    ! First do the i-direction
    !$omp do
    do j = jlo, jhi
      !$omp simd __DELTAQ_ALIGN__
      !dir$ vector aligned
      do i = ilo, ihi
        delta_i_plus(i, j) = q(i + 1, j) - q(i, j)
        ! abs_err = ABS_TOL + REL_TOL * max(abs(q(i + 1, j)), abs(q(i, j)))

        ! ! Error checks
        ! if(abs(delta_i_plus(i, j)) < epsilon(1.0_rk)) then
        !   delta_i_plus(i, j) = 0.0_rk
        ! else if(abs(q(i + 1, j)) < tiny(1.0_rk) .and. abs(q(i, j)) < tiny(1.0_rk)) then
        !   delta_i_plus(i, j) = 0.0_rk
        ! else if(abs(delta_i_plus(i, j)) < abs_err) then
        !   delta_i_plus(i, j) = 0.0_rk
        ! end if
      end do
    end do
    !$omp end do

    ! Now do the j-direction
    !$omp do
    do j = jlo, jhi
      !$omp simd __DELTAQ_ALIGN__
      !dir$ vector aligned
      do i = ilo, ihi
        delta_j_plus(i, j) = q(i, j + 1) - q(i, j)
        ! abs_err = ABS_TOL + REL_TOL * max(abs(q(i, j + 1)), abs(q(i, j)))

        ! ! Error checks
        ! if(abs(delta_j_plus(i, j)) < epsilon(1.0_rk)) then
        !   delta_j_plus(i, j) = 0.0_rk
        ! else if(abs(q(i, j + 1)) < tiny(1.0_rk) .and. abs(q(i, j)) < tiny(1.0_rk)) then
        !   delta_j_plus(i, j) = 0.0_rk
        ! else if(abs(delta_j_plus(i, j)) < abs_err) then
        !   delta_j_plus(i, j) = 0.0_rk
        ! end if
      end do
    end do
    !$omp end do

    ! Since the "minus" value is just the previous cell's "+" value, loop over and copy
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

  subroutine get_r_smoothness(delta_plus, delta_minus, r, r_inv)
    !< Find the smoothness of solution based on nearest neighbor cell averages. Typically referred to
    !< as "r" in the literature. This is sent to the limiters, e.g. phi(r)
    real(rk), dimension(:, :), contiguous, intent(in) :: delta_plus  !< (i,j); Difference operators
    real(rk), dimension(:, :), contiguous, intent(in) :: delta_minus !< (i,j); Difference operators
    real(rk), dimension(:, :), contiguous, intent(inout) :: r        !< (i,j); smoothness r = delta+/delta-
    real(rk), dimension(:, :), contiguous, intent(inout) :: r_inv    !< (i,j); 1/r

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
    !$omp shared(delta_plus, delta_minus, r, r_inv)
    !$omp do
    do j = jlo, jhi
      !$omp simd __DELTAR_ALIGN__
      !dir$ vector aligned
      do i = ilo, ihi
        r(i, j) = (delta_plus(i, j) + EPS) / (delta_minus(i, j) + EPS)
        r_inv(i, j) = (delta_minus(i, j) + EPS) / (delta_plus(i, j) + EPS)
      end do
    end do
    !$omp end do

    ! !$omp do
    ! do j = jlo, jhi
    !   !$omp simd __DELTAR_ALIGN__
    !   !dir$ vector aligned
    !   do i = ilo, ihi
    !     if(abs(r(i, j)) < 1e-10_rk) r(i, j) = 0.0_rk
    !     if(abs(r_inv(i, j)) < 1e-10_rk) r_inv(i, j) = 0.0_rk
    !   end do
    ! end do
    ! !$omp end do
    !$omp end parallel
  end subroutine get_r_smoothness

  subroutine limit(r, psi, name)
    !< Apply the flux limiter based on the smoothness
    real(rk), dimension(:, :), contiguous, intent(in) :: r !< smoothness
    real(rk), dimension(:, :), contiguous, intent(inout) :: psi !< flux limiter value
    character(len=*), intent(in) :: name !< name of the limiter

    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(RK), parameter :: PSI_EPS = 1e-5_rk
    !dir$ assume_aligned r: __ALIGNBYTES__
    !dir$ assume_aligned psi: __ALIGNBYTES__

    ilo = lbound(r, dim=1)
    ihi = ubound(r, dim=1)
    jlo = lbound(r, dim=2)
    jhi = ubound(r, dim=2)
    select case(trim(name))
    case default
      write(std_err, '(a)') "Error: Unsupported limiter type in mod_edge_interpolator%limit() '"//trim(name)//"'"
      error stop "Error: Unsupported limiter type in mod_edge_interpolator%limit()"
    case('none')
      psi = 1.0_rk
    case('minmod')
      !$omp parallel default(none), &
      !$omp firstprivate(ilo, ihi, jlo, jhi) &
      !$omp private(i, j) shared(psi, r)
      !$omp do
      do j = jlo, jhi
        !dir$ vector aligned
        !$omp simd __PSI_ALIGN__
        do i = ilo, ihi
          psi(i, j) = max(0.0_rk, min(r(i, j), 1.0_rk))
        end do
      end do

      !$omp do
      do j = jlo, jhi
        !dir$ vector aligned
        !$omp simd __PSI_ALIGN__
        do i = ilo, ihi
          if(r(i, j) < PSI_EPS) then
            psi(i, j) = 0.0_rk
          else if(r(i, j) > 1.0_rk) then
            psi(i, j) = 1.0_rk
          endif
        end do
      end do
      !$omp end do
      !$omp end parallel

    case('superbee')
      !$omp parallel default(none), &
      !$omp firstprivate(ilo, ihi, jlo, jhi) &
      !$omp private(i, j) shared(psi, r)
      !$omp do
      do j = jlo, jhi
        !dir$ vector aligned
        !$omp simd __PSI_ALIGN__
        do i = ilo, ihi
          psi(i, j) = max(0.0_rk, min(2.0_rk * r(i, j), 1.0_rk), min(r(i, j), 2.0_rk))
        end do
      end do

      !$omp do
      do j = jlo, jhi
        !dir$ vector aligned
        !$omp simd __PSI_ALIGN__
        do i = ilo, ihi
          if(r(i, j) < PSI_EPS) then
            psi(i, j) = 0.0_rk
          else if(r(i, j) > 2.0_rk) then
            psi(i, j) = 2.0_rk
          endif
        end do
      end do
      !$omp end do
      !$omp end parallel
    case('van_leer')
      !$omp parallel default(none), &
      !$omp firstprivate(ilo, ihi, jlo, jhi) &
      !$omp private(i, j) shared(psi, r)
      !$omp do
      do j = jlo, jhi
        !dir$ vector aligned
        !$omp simd __PSI_ALIGN__
        do i = ilo, ihi
          psi(i, j) = (r(i, j) + abs(r(i, j))) / (1.0_rk + abs(r(i, j)))
        end do
      end do

      !$omp do
      do j = jlo, jhi
        !dir$ vector aligned
        !$omp simd __PSI_ALIGN__
        do i = ilo, ihi
          if(r(i, j) < PSI_EPS) psi(i, j) = 0.0_rk
        end do
      end do
      !$omp end do
      !$omp end parallel

    end select
  end subroutine limit

end module mod_edge_interpolator
