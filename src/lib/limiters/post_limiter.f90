module mod_post_limiter
  !<
  !< [1] Keiichi Kitamura, Atsushi Hashimoto, Simple a posteriori slope limiter (Post Limiter) for high resolution and efficient flow computations
  !<     https://doi.org/10.1016/j.jcp.2017.04.002

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_gnoffo_aux_limiter, only: gnoffo_aux_limiter
  use mod_slope_limiter, only: minmod_limit
  implicit none

contains

  function limit(U, U_neighbor) result(q)
    real(rk), dimension(4), intent(in) :: U !< primitive variables of the current cell
    real(rk), dimension(4), intent(in) :: U_neighbor !< primitive variables of the neighbor cell

    real(rk), dimension(4) :: q_unlimited !< non-limited primitive variables
    real(rk), dimension(4) :: q_limited   !< limited primitive variables
    real(rk), dimension(4) :: q           !< final version of the primitive variables

    logical :: positivity !< density and pressure positive
    logical :: dmp        !< discrete maximum principle

    positivity = .false.
    if(rho_ij > 0.0_rk .and. p_ij > 0.0_rk) positivity = .true.

    dmp = .false.
    if(min(rho_i, rho_j) < rho_ij .and. &
       min(p_i, p_j) < p_ij .and. &
       rho_ij < max(rho_i, rho_j) .and. &
       p_ij < max(p_i, p_j)) then
      dmp = .true.
    end if

    phi_final = phi_g + (1.0_rk - phi_g) * phi_limited
  end function

end module mod_post_limiter
