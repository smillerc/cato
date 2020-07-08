module mod_flux_array
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_vector, only: vector_t
  use mod_eos, only: eos
  use mod_floating_point_utils, only: neumaier_sum

  implicit none

  private
  public :: get_fluxes, flux_array_t

  type :: flux_array_t
    real(rk), dimension(:, :, :), allocatable :: F
    real(rk), dimension(:, :, :), allocatable :: G
  contains
  end type flux_array_t

contains

  function get_fluxes(rho, u, v, p, lbounds) result(flux)
    !< Implementation of the flux tensor (H) construction. This requires the primitive variables
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: rho !< (i,j)
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: u   !< (i,j)
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: v   !< (i,j)
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: p   !< (i,j)

    type(flux_array_t) :: flux

    ! real(rk), dimension(1:, 1:, lbounds(1):, lbounds(2):), intent(inout) :: H !< ((Fi, Gj), (1:4), i, j) flux array

    real(rk), dimension(:, :), allocatable :: E !< total energy
    integer(ik) :: i, j
    integer(ik) :: ilo
    integer(ik) :: ihi
    integer(ik) :: jlo
    integer(ik) :: jhi

    real(rk) :: rho_u, rho_v, rho_u_v, rho_E, rho_E_plus_p
    real(rk) :: rho_u_sq, rho_u_sq_plus_p
    real(rk) :: rho_v_sq, rho_v_sq_plus_p

    ilo = lbound(rho, dim=1)
    ihi = ubound(rho, dim=1)
    jlo = lbound(rho, dim=2)
    jhi = ubound(rho, dim=2)

    allocate(flux%F(4, ilo:ihi, jlo:jhi))
    allocate(flux%G(4, ilo:ihi, jlo:jhi))

    ! get the total energy
    allocate(E, mold=rho)
    call eos%total_energy(rho, u, v, p, E)

    ! The flux tensor is H = Fi + Gj

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j, rho_u, rho_u_sq_plus_p, rho_u_v, rho_E_plus_p, rho_u_sq, rho_v_sq, rho_v_sq_plus_p) &
    !$omp shared(flux, rho, u, v, p, E)
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi

        rho_u = rho(i, j) * u(i, j)
        rho_u_v = rho(i, j) * u(i, j) * v(i, j)
        rho_E_plus_p = rho(i, j) * E(i, j) + p(i, j)

        ! rho u^2 + p
        rho_u_sq = rho(i, j) * u(i, j)**2
        if(abs(rho_u_sq - p(i, j)) < epsilon(1.0_rk)) then
          rho_u_sq_plus_p = 0.0_rk
        else
          rho_u_sq_plus_p = rho_u_sq + p(i, j)
        end if

        ! rho v^2 + p
        rho_v_sq = rho(i, j) * v(i, j)**2
        if(abs(rho_v_sq - p(i, j)) < epsilon(1.0_rk)) then
          rho_v_sq_plus_p = 0.0_rk
        else
          rho_v_sq_plus_p = rho_v_sq + p(i, j)
        end if

        flux%F(1, i, j) = rho_u
        flux%F(2, i, j) = rho_u_sq_plus_p
        flux%F(3, i, j) = rho_u_v
        flux%F(4, i, j) = u(i, j) * rho_E_plus_p

        flux%G(1, i, j) = rho(i, j) * v(i, j)
        flux%G(2, i, j) = rho_u_v
        flux%G(3, i, j) = rho_v_sq_plus_p
        flux%G(4, i, j) = v(i, j) * rho_E_plus_p

      end do
    end do
    !$omp end do
    !$omp end parallel
  end function
end module mod_flux_array
