module mod_edge_reconstruction
  !> Summary: Provide procedures to reconstruct the cell interface values
  !> Date: 05/08/2020
  !> Author: Sam Miller
  !> Notes:
  !> References:
  !>   [1] M. Berger, M. Aftosmis, S. Muman, "Analysis of Slope Limiters on Irregular Grids",
  !>       43rd AIAA Aerospace Sciences Meeting and Exhibit (2005), https://doi.org/10.2514/6.2005-490
  !>
  !>   [2] K.H. Kim, C. Kim, "Accurate, efficient and monotonic numerical methods for multi-dimensional compressible flows Part II: Multi-dimensional limiting process",
  !>       Journal of Computational Physics 208 (2005) 570–615, https://doi.org/10.1016/j.jcp.2005.02.022

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_globals, only: debug_print, plot_limiters
  use mod_flux_limiter, only: flux_limiter_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_globals, only: MACHINE_EPS, n_ghost_layers

  implicit none

  logical, parameter :: filter_small = .true.
  real(rk), parameter :: SMALL_DIFF = 1e-9_rk

  private
  public :: reconstruct_edge_values

contains

  real(rk) function get_delta(a, b) result(delta)
    !< Find the delta in the solution, e.g. delta = a - b. This checks for numbers
    !< near 0 and when a and b are very similar in magnitude. The aim is to avoid
    !< catastrophic cancellation and very small numbers that are essentially 0 for this scenario

    real(rk), intent(in) :: a, b
    real(rk), parameter :: rel_tol = 1e-12_rk     !< relative error tolerance
    real(rk), parameter :: abs_tol = tiny(1.0_rk) !< absolute error tolerance
    real(rk) :: abs_err !< absolute error

    delta = a - b
    abs_err = abs_tol + rel_tol * max(abs(a), abs(b))

    ! write(*,'(8(es16.6))') a, b, delta, abs_err

    if(abs(delta) < epsilon(1.0_rk)) then
      delta = 0.0_rk
      return
    end if

    if(abs(a) < tiny(1.0_rk) .and. abs(b) < tiny(1.0_rk)) then
      delta = 0.0_rk
      return
    end if

    if(abs(delta) < abs_err) then
      ! write(*,'(8(es16.6))') a, b, delta, abs_err
      delta = 0.0_rk
      return
    end if
  end function

  subroutine get_smoothness(minus, current, plus, R, R_inv)
    real(rk), intent(in) :: minus   !< q(i-1)
    real(rk), intent(in) :: current !< q(i)
    real(rk), intent(in) :: plus    !< q(i+1)
    real(rk), intent(out) :: R      !< R
    real(rk), intent(out) :: R_inv  !< 1/R

    real(rk) :: delta_plus, delta_minus
    real(rk), parameter :: infinity = 1e20_rk

    delta_minus = get_delta(current, minus) ! q(i) - q(i-1)
    delta_plus = get_delta(plus, current)   ! q(i+1) - q(i)

    if(abs(delta_plus - delta_minus) < epsilon(1.0_rk) .or. & ! deltas are the same
       (abs(delta_plus) < tiny(1.0_rk) .and. abs(delta_minus) < tiny(1.0_rk)) & ! both are 0
       ) then
      R = 1.0_rk
      R_inv = 1.0_rk
    else if(abs(delta_minus) < tiny(1.0_rk)) then ! delta- is 0
      R = infinity
      R_inv = 0.0_rk
    else if(abs(delta_plus) < tiny(1.0_rk)) then ! delta+ is 0
      R = 0.0_rk
      R_inv = infinity
    else
      R = delta_plus / delta_minus
      R_inv = 1.0_rk / R
    end if

  end subroutine get_smoothness

  subroutine reconstruct_edge_values(q, lbounds, limiter, edge_values)
    !< Reconstruct the cell interface values, e.g. q_i-1/2, q_i+1/2. This assumes a cartesian
    !< structured square grid

    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
    !< (i,j); primitive variable to reconstruct at the edge

    real(rk), dimension(:, :, :), allocatable, intent(out) :: edge_values
    !<((bottom, right, top, left), i, j); reconstructed edge values

    type(slope_limiter_t), intent(in) :: limiter !< slope limiter used to reconstruct the edge interface

    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: ilo_bc, ihi_bc, jlo_bc, jhi_bc
    real(rk) :: R_i     !< R   smoothness indicator for i
    real(rk) :: R_i_inv !< 1/R smoothness indicator for i
    real(rk) :: R_j     !< R   smoothness indicator for j
    real(rk) :: R_j_inv !< 1/R smoothness indicator for j

    real(rk) :: phi_top    !< limiter for the top edge
    real(rk) :: phi_bottom !< limiter for the bottom edge
    real(rk) :: phi_left   !< limiter for the left edge
    real(rk) :: phi_right  !< limiter for the right edge

    call debug_print('Running reconstruct_edge_values()', __FILE__, __LINE__)

    ilo_bc = lbound(q, dim=1)
    ihi_bc = ubound(q, dim=1)
    jlo_bc = lbound(q, dim=2)
    jhi_bc = ubound(q, dim=2)

    ilo = ilo_bc + n_ghost_layers
    ihi = ihi_bc - n_ghost_layers
    jlo = jlo_bc + n_ghost_layers
    jhi = jhi_bc - n_ghost_layers

    allocate(edge_values(4, ilo:ihi, jlo:jhi))

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp private(R_i, R_i_inv, R_j, R_j_inv) &
    !$omp private(phi_bottom, phi_top, phi_left, phi_right) &
    !$omp shared(q, limiter, edge_values)
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi

        call get_smoothness(q(i - 1, j), q(i, j), q(i + 1, j), R_i, R_i_inv)
        call get_smoothness(q(i, j - 1), q(i, j), q(i, j + 1), R_j, R_j_inv)

        ! slope limiters
        phi_top = 0.0_rk
        phi_bottom = 0.0_rk
        phi_left = 0.0_rk
        phi_right = 0.0_rk

        ! (i, j-1/2), bottom edge
        phi_bottom = limiter%limit(R_j_inv)
        edge_values(1, i, j) = q(i, j) - 0.5_rk * phi_bottom * ((q(i, j + 1) - q(i, j - 1)) / 2.0_rk)

        ! (i+1/2, j), right edge
        phi_right = limiter%limit(R_i)
        edge_values(2, i, j) = q(i, j) + 0.5_rk * phi_right * ((q(i + 1, j) - q(i - 1, j)) / 2.0_rk)

        ! (i, j+1/2), top edge
        phi_top = limiter%limit(R_j)
        edge_values(3, i, j) = q(i, j) + 0.5_rk * phi_top * ((q(i, j + 1) - q(i, j - 1)) / 2.0_rk)

        ! (i-1/2, j), left edge
        phi_left = limiter%limit(R_i_inv)
        edge_values(4, i, j) = q(i, j) - 0.5_rk * phi_left * ((q(i + 1, j) - q(i - 1, j)) / 2.0_rk)

      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine reconstruct_edge_values

  ! subroutine reconstruct_edge_values(q, lbounds, limiter, edge_values)
  !   !< Reconstruct the cell interface values, e.g. q_i-1/2, q_i+1/2. This assumes a cartesian
  !   !< structured square grid

  !   integer(ik), dimension(2), intent(in) :: lbounds
  !   real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: q
  !   !< (i,j); primitive variable to reconstruct at the edge

  !   real(rk), dimension(:, :, :), allocatable, intent(out) :: edge_values
  !   !<((bottom, right, top, left), i, j); reconstructed edge values

  !   type(flux_limiter_t), intent(in) :: limiter !< flux limiter used to reconstruct the edge interface
  !   ! type(slope_limiter_t), intent(in) :: limiter !< flux limiter used to reconstruct the edge interface

  !   integer(ik) :: i, j
  !   integer(ik) :: ilo, ihi, jlo, jhi
  !   integer(ik) :: ilo_bc, ihi_bc, jlo_bc, jhi_bc
  !   real(rk) :: R_i !< smoothness indicator
  !   real(rk) :: R_j !< smoothness indicator
  !   real(rk) :: delta_i_minus !< q_i - q_{i-1}
  !   real(rk) :: delta_j_minus !< q_j - q_{j-1}
  !   real(rk) :: delta_i_plus  !< q_{i+1} - q_i
  !   real(rk) :: delta_j_plus  !< q_{j+1} - q_j
  !   real(rk) :: phi_top    !< limiter for the top edge
  !   real(rk) :: phi_bottom !< limiter for the bottom edge
  !   real(rk) :: phi_left   !< limiter for the left edge
  !   real(rk) :: phi_right  !< limiter for the right edge

  !   call debug_print('Running reconstruct_edge_values()', __FILE__, __LINE__)

  !   ilo_bc = lbound(q, dim=1)
  !   ihi_bc = ubound(q, dim=1)
  !   jlo_bc = lbound(q, dim=2)
  !   jhi_bc = ubound(q, dim=2)

  !   ilo = ilo_bc + n_ghost_layers
  !   ihi = ihi_bc - n_ghost_layers
  !   jlo = jlo_bc + n_ghost_layers
  !   jhi = jhi_bc - n_ghost_layers

  !   allocate(edge_values(4, ilo:ihi, jlo:jhi))

  !   !$omp parallel default(none), &
  !   !$omp firstprivate(ilo, ihi, jlo, jhi) &
  !   !$omp private(i, j, R_i, R_j) &
  !   !$omp private(delta_i_minus, delta_j_minus, delta_i_plus, delta_j_plus) &
  !   !$omp private(phi_bottom, phi_top, phi_left, phi_right) &
  !   !$omp shared(q, limiter, edge_values)
  !   !$omp do
  !   do j = jlo, jhi
  !     do i = ilo, ihi
  !       delta_i_minus = get_delta(q(i, j), q(i - 1, j))
  !       delta_j_minus = get_delta(q(i, j), q(i, j - 1))
  !       delta_i_plus = get_delta(q(i + 1, j), q(i, j))
  !       delta_j_plus = get_delta(q(i, j + 1), q(i, j))

  !       phi_top    = 0.0_rk
  !       phi_bottom  = 0.0_rk
  !       phi_left  = 0.0_rk
  !       phi_right = 0.0_rk

  !       ! (i, j-1/2), bottom
  !       if(abs(delta_j_plus) > 0.0_rk) then
  !         R_j = delta_j_minus / delta_j_plus ! this is 1/R of the normal definition
  !         phi_bottom = limiter%limit(R_j)
  !         edge_values(1, i, j) = q(i, j) - 0.5_rk * phi_bottom * delta_j_plus
  !       else
  !         edge_values(1, i, j) = q(i, j)
  !       end if

  !       ! (i+1/2, j), right
  !       if(abs(delta_i_minus) > 0.0_rk) then
  !         R_i = delta_i_plus / delta_i_minus
  !         phi_right = limiter%limit(R_i)
  !         edge_values(2, i, j) = q(i, j) + 0.5_rk * phi_right * delta_i_minus
  !       else
  !         edge_values(2, i, j) = q(i, j)
  !       end if

  !       ! (i, j+1/2), top
  !       if(abs(delta_j_minus) > 0.0_rk) then
  !         R_j = delta_j_plus / delta_j_minus
  !         phi_top = limiter%limit(R_j)
  !         edge_values(3, i, j) = q(i, j) + 0.5_rk * phi_top * delta_j_minus
  !       else
  !         edge_values(3, i, j) = q(i, j)
  !       end if

  !       ! (i-1/2, j), left
  !       if(abs(delta_i_plus) > 0.0_rk) then
  !         R_i = delta_i_minus / delta_i_plus ! this is 1/R of the normal definition
  !         phi_left = limiter%limit(R_i)
  !         edge_values(4, i, j) = q(i, j) - 0.5_rk * phi_left * delta_i_plus
  !       else
  !         edge_values(4, i, j) = q(i, j)
  !       end if

  !       ! write(*,'(4(es16.6))') delta_j_plus, delta_j_minus, delta_i_minus, delta_i_plus
  !       ! write(*,'(4(es16.6))') phi_bottom,   phi_top,       phi_right,     phi_left
  !       ! print*
  !     end do
  !   end do
  !   !$omp end do
  !   !$omp end parallel
  ! end subroutine reconstruct_edge_values

end module mod_edge_reconstruction
