module mod_mlp_5th_order
  !< Summary: Provide class for 5th order MLP edge interpolation
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
  use mod_flux_limiter, only: flux_limiter_t, smoothness, delta
  use mod_mlp_baseline, only: mlp_baseline_t
  use mod_globals, only: n_ghost_layers, debug_print

  implicit none
  private
  public :: mlp_5th_order_t

  type, extends(mlp_baseline_t) :: mlp_5th_order_t
    !< 5th order edge interpolation with the MLP limiter
    type(flux_limiter_t) :: limiter
  contains
    procedure, public :: initialize
    procedure, public :: interpolate_edge_values
    ! procedure, public, nopass :: get_beta => beta_5th_order
  end type mlp_5th_order_t

contains
  subroutine initialize(self, limiter)
    class(mlp_5th_order_t), intent(inout) :: self
    character(len=*), intent(in) :: limiter
    self%limiter_name = 'MLP5'
    self%order = 5
    ! self%limiter = slope_limiter_t(trim(limiter))
  end subroutine initialize

  subroutine interpolate_edge_values(self, q, lbounds, edge_values)
    class(mlp_5th_order_t), intent(in) :: self
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

    real(rk), dimension(:, :), allocatable :: alpha_L_xi   !< alpha_L,xi ; MLP alpha factor for i-direction
    real(rk), dimension(:, :), allocatable :: alpha_R_xi   !< alpha_R,xi ; MLP alpha factor for i-direction
    real(rk), dimension(:, :), allocatable :: alpha_L_eta  !< alpha_L,eta; MLP alpha factor for j-direction
    real(rk), dimension(:, :), allocatable :: alpha_R_eta  !< alpha_R,eta; MLP alpha factor for j-direction

    real(rk), dimension(:, :), allocatable :: tan_theta_i !< tan(theta_i); Representative MLP angle
    real(rk), dimension(:, :), allocatable :: tan_theta_j !< tan(theta_j); Representative MLP angle

    real(rk) :: phi_top    !< limiter for the top edge
    real(rk) :: phi_bottom !< limiter for the bottom edge
    real(rk) :: phi_left   !< limiter for the left edge
    real(rk) :: phi_right  !< limiter for the right edge

    real(rk) :: delta_i_plus, delta_i_minus, delta_j_plus, delta_j_minus

    call debug_print('Running mlp_3rd_order_t%interpolate_edge_values()', __FILE__, __LINE__)

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

    allocate(alpha_L_xi(ilo - 1:ihi + 1, jlo - 1:jhi + 1))
    allocate(alpha_R_xi(ilo - 1:ihi + 1, jlo - 1:jhi + 1))
    allocate(alpha_L_eta(ilo - 1:ihi + 1, jlo - 1:jhi + 1))
    allocate(alpha_R_eta(ilo - 1:ihi + 1, jlo - 1:jhi + 1))

    allocate(tan_theta_i(ilo - 1:ihi + 1, jlo - 1:jhi + 1))
    allocate(tan_theta_j(ilo - 1:ihi + 1, jlo - 1:jhi + 1))

    ! get tan(theta_j) and tan(theta_j)
    call self%get_rep_angles(q, lbounds=lbound(q), tan_theta_i=tan_theta_i, tan_theta_j=tan_theta_j)

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp private(phi_bottom, phi_top, phi_left, phi_right) &
    !$omp private(delta_i_plus, delta_i_minus, delta_j_plus, delta_j_minus) &
    !$omp shared(r_L_i, r_R_i, r_L_j, r_R_j) &
    !$omp shared(beta_L_i, beta_R_i, beta_L_j, beta_R_j) &
    !$omp shared(alpha_L_xi, alpha_R_xi, alpha_L_eta, alpha_R_eta) &
    !$omp shared(tan_theta_i, tan_theta_j) &
    !$omp shared(q, self, edge_values)

    !$omp do
    do j = jlo - 1, jhi + 1
      do i = ilo - 1, ihi + 1
        r_L_i(i, j) = smoothness(q(i - 1, j), q(i, j), q(i + 1, j))
        r_R_i(i, j) = 1.0_rk / r_L_i(i, j)
        r_L_j(i, j) = smoothness(q(i, j - 1), q(i, j), q(i, j + 1))
        r_R_j(i, j) = 1.0_rk / r_L_j(i, j)
      end do
    end do
    !$omp end do
    !$omp barrier

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        ! 5th order interpolation function
        beta_L_i(i, j) = ((-2.0_rk / r_L_i(i - 1, j)) + 11.0_rk + &
                          24.0_rk * r_L_i(i, j) - (3.0_rk * r_L_i(i, j) * r_L_i(i + 1, j))) / 30.0_rk
        beta_R_i(i, j) = ((-2.0_rk / r_R_i(i - 1, j)) + 11.0_rk + &
                          24.0_rk * r_R_i(i, j) - (3.0_rk * r_R_i(i, j) * r_R_i(i + 1, j))) / 30.0_rk
        beta_L_j(i, j) = ((-2.0_rk / r_L_j(i, j - 1)) + 11.0_rk + &
                          24.0_rk * r_L_j(i, j) - (3.0_rk * r_L_j(i, j) * r_L_j(i, j + 1))) / 30.0_rk
        beta_R_j(i, j) = ((-2.0_rk / r_R_j(i, j - 1)) + 11.0_rk + &
                          24.0_rk * r_R_j(i, j) - (3.0_rk * r_R_j(i, j) * r_R_j(i, j + 1))) / 30.0_rk
      end do
    end do
    !$omp end do
    !$omp barrier

    !$omp do
    do j = jlo + 1, jhi - 1
      do i = ilo + 1, ihi - 1

        ! i+1/2
        if(abs(r_R_i(i + 1, j)) > tiny(1.0_rk)) then
          alpha_L_xi(i, j) = self%g((2.0_rk * max(1.0_rk, r_L_i(i, j)) &
                                     * (1.0_rk + max(0.0_rk,(tan_theta_i(i + 1, j) / r_R_i(i + 1, j))))) / &
                                    (1.0_rk + tan_theta_i(i, j)))
        else
          alpha_L_xi(i, j) = self%g((2.0_rk * max(1.0_rk, r_L_i(i, j))) / (1.0_rk + tan_theta_i(i, j)))
        end if

        ! i-1/2
        if(abs(r_L_i(i - 1, j)) > tiny(1.0_rk)) then
          alpha_R_xi(i, j) = self%g((2.0_rk * max(1.0_rk, r_R_i(i, j)) &
                                     * (1.0_rk + max(0.0_rk,(tan_theta_i(i - 1, j) / r_L_i(i - 1, j))))) / &
                                    (1.0_rk + tan_theta_i(i, j)))
        else
          alpha_R_xi(i, j) = self%g((2.0_rk * max(1.0_rk, r_R_i(i, j))) / (1.0_rk + tan_theta_i(i, j)))

        end if

        ! j+1/2
        if(abs(r_R_j(i, j + 1)) > tiny(1.0_rk)) then
          alpha_L_eta(i, j) = self%g((2.0_rk * max(1.0_rk, r_L_j(i, j)) &
                                      * (1.0_rk + max(0.0_rk,(tan_theta_j(i, j + 1) / r_R_j(i, j + 1))))) / &
                                     (1.0_rk + tan_theta_j(i, j)))
        else
          alpha_L_eta(i, j) = self%g((2.0_rk * max(1.0_rk, r_L_j(i, j))) / (1.0_rk + tan_theta_j(i, j)))
        end if

        ! j-1/2
        if(abs(r_L_j(i, j - 1)) > tiny(1.0_rk)) then
          alpha_R_eta(i, j) = self%g((2.0_rk * max(1.0_rk, r_R_j(i, j)) &
                                      * (1.0_rk + max(0.0_rk,(tan_theta_j(i, j - 1) / r_L_j(i, j - 1))))) / &
                                     (1.0_rk + tan_theta_j(i, j)))
        else
          alpha_R_eta(i, j) = self%g((2.0_rk * max(1.0_rk, r_R_j(i, j))) / (1.0_rk + tan_theta_j(i, j)))
        end if
      end do
    end do
    !$omp end do
    !$omp barrier

    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi
        ! (i, j-1/2), bottom edge -> corresponds to the "R" side of the interface, thus the "R" terms
        phi_bottom = max(0.0_rk, min(alpha_R_eta(i, j) * r_R_j(i, j), alpha_R_eta(i, j), beta_R_j(i, j)))
        delta_j_plus = delta(q(i, j + 1), q(i, j)) ! q(i,j+1) - q(i,j)
        edge_values(1, i, j) = q(i, j) - 0.5_rk * phi_bottom * delta_j_plus

        ! (i+1/2, j), right edge -> corresponds to the "L" side of the interface, thus the "L" terms
        phi_right = max(0.0_rk, min(alpha_L_xi(i, j) * r_L_i(i, j), alpha_L_xi(i, j), beta_L_i(i, j)))
        delta_i_minus = delta(q(i, j), q(i - 1, j)) ! q(i,j) - q(i-1,j)
        edge_values(2, i, j) = q(i, j) + 0.5_rk * phi_right * delta_i_minus

        ! (i, j+1/2), top edge -> corresponds to the "L" side of the interface, thus the "L" terms
        phi_top = max(0.0_rk, min(alpha_L_eta(i, j) * r_L_j(i, j), alpha_L_eta(i, j), beta_L_j(i, j)))
        delta_j_minus = delta(q(i, j), q(i, j - 1)) ! q(i,j) - q(i,j-1)
        edge_values(3, i, j) = q(i, j) + 0.5_rk * phi_top * delta_j_minus

        ! (i-1/2, j), left edge -> corresponds to the "R" side of the interface, thus the "R" terms
        phi_left = max(0.0_rk, min(alpha_R_xi(i, j) * r_R_i(i, j), alpha_R_xi(i, j), beta_R_i(i, j)))
        delta_i_plus = delta(q(i + 1, j), q(i, j)) ! q(i+1,j) - q(i,j)
        edge_values(4, i, j) = q(i, j) - 0.5_rk * phi_left * delta_i_plus
      end do
    end do
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

    deallocate(alpha_L_xi)
    deallocate(alpha_R_xi)
    deallocate(alpha_L_eta)
    deallocate(alpha_R_eta)

    deallocate(tan_theta_i)
    deallocate(tan_theta_j)

  end subroutine interpolate_edge_values

end module mod_mlp_5th_order
