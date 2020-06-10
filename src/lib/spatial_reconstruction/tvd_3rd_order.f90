module mod_tvd_3rd_order
  !< Summary: Provide class for 3rd order TVD edge interpolation
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
  use mod_flux_limiter, only: flux_limiter_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_edge_interp, only: edge_iterpolator_t
  use mod_globals, only: n_ghost_layers

  implicit none
  private
  public :: tvd_3rd_order_t

  type, extends(edge_iterpolator_t) :: tvd_3rd_order_t
    !< 3rd order edge interpolation with TVD filtering
    type(slope_limiter_t) :: limiter
  contains
    procedure, public :: initialize
    ! procedure, nopass, public :: get_beta => beta_3rd_order
    procedure, public :: interpolate_edge_values
  end type tvd_3rd_order_t

contains

  subroutine initialize(self, limiter)
    !< Constructor for tvd_3rd_order_t
    class(tvd_3rd_order_t), intent(inout) :: self
    character(len=*), intent(in) :: limiter
    self%limiter_name = trim(limiter)
    self%order = 3
    self%limiter = slope_limiter_t(trim(limiter))
  end subroutine initialize

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
    real(rk), dimension(:, :), allocatable :: r_i     !< r   smoothness indicator for i
    real(rk), dimension(:, :), allocatable :: r_i_inv !< 1/r smoothness indicator for i
    real(rk), dimension(:, :), allocatable :: r_j     !< r   smoothness indicator for j
    real(rk), dimension(:, :), allocatable :: r_j_inv !< 1/r smoothness indicator for j

    real(rk), dimension(:, :), allocatable :: beta_left  !< r   smoothness indicator for i
    real(rk), dimension(:, :), allocatable :: beta_right !< 1/r smoothness indicator for i
    real(rk), dimension(:, :), allocatable :: beta_top    !< r   smoothness indicator for j
    real(rk), dimension(:, :), allocatable :: beta_bottom  !< 1/r smoothness indicator for j

    real(rk) :: phi_top    !< limiter for the top edge
    real(rk) :: phi_bottom !< limiter for the bottom edge
    real(rk) :: phi_left   !< limiter for the left edge
    real(rk) :: phi_right  !< limiter for the right edge

    real(rk) :: delta_i_plus, delta_i_minus, delta_j_plus, delta_j_minus

    ! call debug_print('Running interpolate_edge_values()', __FILE__, __LINE__)

    ilo_bc = lbound(q, dim=1)
    ihi_bc = ubound(q, dim=1)
    jlo_bc = lbound(q, dim=2)
    jhi_bc = ubound(q, dim=2)

    ilo = ilo_bc + n_ghost_layers
    ihi = ihi_bc - n_ghost_layers
    jlo = jlo_bc + n_ghost_layers
    jhi = jhi_bc - n_ghost_layers

    allocate(edge_values(4, ilo:ihi, jlo:jhi))
    allocate(r_i(ilo:ihi, jlo:jhi))
    allocate(r_i_inv(ilo:ihi, jlo:jhi))
    allocate(r_j(ilo:ihi, jlo:jhi))
    allocate(r_j_inv(ilo:ihi, jlo:jhi))
    allocate(beta_left(ilo:ihi, jlo:jhi))
    allocate(beta_right(ilo:ihi, jlo:jhi))
    allocate(beta_top(ilo:ihi, jlo:jhi))
    allocate(beta_bottom(ilo:ihi, jlo:jhi))

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp private(phi_bottom, phi_top, phi_left, phi_right) &
    !$omp private(delta_i_plus, delta_i_minus, delta_j_plus, delta_j_minus) &
    !$omp shared(r_i, r_i_inv, r_j, r_j_inv) &
    !$omp shared(beta_bottom, beta_top, beta_left, beta_right) &
    !$omp shared(q, self, edge_values)

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        ! Get r_i, i.e. smoothness across i direction
        call self%get_smoothness(q(i - 1, j), q(i, j), q(i + 1, j), r_i(i, j), r_i_inv(i, j))

        ! Get r_j, i.e. smoothness across j direction
        call self%get_smoothness(q(i, j - 1), q(i, j), q(i, j + 1), r_j(i, j), r_j_inv(i, j))
      end do
    end do
    !$omp end do

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        ! 3rd order interpolation function
        beta_left(i, j) = (1.0_rk + 2.0_rk * r_i(i, j)) / 3.0_rk
        beta_right(i, j) = (1.0_rk + 2.0_rk * r_i_inv(i, j)) / 3.0_rk
        beta_top(i, j) = (1.0_rk + 2.0_rk * r_j(i, j)) / 3.0_rk
        beta_bottom(i, j) = (1.0_rk + 2.0_rk * r_j_inv(i, j)) / 3.0_rk

        ! write(*,'(4(es16.6))') r_i(i, j), r_i_inv(i, j), r_j(i, j), r_j_inv(i, j)
        ! write(*,'(4(es16.6))') beta_left(i, j), beta_right(i, j), beta_top(i, j), beta_bottom(i, j)
      end do
    end do
    !$omp end do
    ! error stop

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi

        ! (i, j+1/2), top edge
        phi_top = max(0.0_rk, min(2.0_rk, 2.0_rk * r_j(i, j), beta_top(i, j)))
        delta_j_plus = self%get_delta(q(i, j), q(i, j - 1)) ! q(i,j) - q(i,j-1)
        edge_values(3, i, j) = q(i, j) + 0.5_rk * phi_top * delta_j_plus

        ! (i, j-1/2), bottom edge
        phi_bottom = max(0.0_rk, min(2.0_rk, 2.0_rk * r_j_inv(i, j), beta_bottom(i, j)))
        delta_j_minus = self%get_delta(q(i, j + 1), q(i, j)) ! q(i,j+1) - q(i,j)
        edge_values(1, i, j) = q(i, j) - 0.5_rk * phi_bottom * delta_j_minus

        ! (i+1/2, j), right edge
        phi_right = max(0.0_rk, min(2.0_rk, 2.0_rk * r_i(i, j), beta_right(i, j)))
        delta_i_plus = self%get_delta(q(i, j), q(i - 1, j)) ! q(i,j) - q(i-1,j)
        edge_values(2, i, j) = q(i, j) + 0.5_rk * phi_right * delta_i_plus

        ! (i-1/2, j), left edge
        phi_left = max(0.0_rk, min(2.0_rk, 2.0_rk * r_i_inv(i, j), beta_left(i, j)))
        delta_i_minus = self%get_delta(q(i + 1, j), q(i, j)) ! q(i+1,j) - q(i,j)
        edge_values(4, i, j) = q(i, j) - 0.5_rk * phi_left * delta_i_minus

        ! write(*,'(5(es16.6))') q(i,j), edge_values(1, i, j), edge_values(2, i, j), edge_values(3, i, j), edge_values(4, i, j)
      end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(r_i)
    deallocate(r_i_inv)
    deallocate(r_j)
    deallocate(r_j_inv)

    deallocate(beta_left)
    deallocate(beta_right)
    deallocate(beta_top)
    deallocate(beta_bottom)

  end subroutine
end module mod_tvd_3rd_order
