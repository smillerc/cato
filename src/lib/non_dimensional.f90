module mod_nondimensionalization
  !> Summary: Provide scaling factors to non-dimensionalize the Euler equations
  !> Date: 04/25/2020
  !> Author: Sam Miller
  !> Notes:
  !> References:
  !      [1]

  use iso_fortran_env, only: rk => real64

  implicit none

  real(rk), protected :: t_0 = 1.0_rk   !< scale factor for time
  real(rk), protected :: rho_0 = 1.0_rk !< scale factor for density
  real(rk), protected :: l_0 = 1.0_rk   !< scale factor for length
  real(rk), protected :: v_0 = 1.0_rk   !< scale factor for velocity
  real(rk), protected :: p_0 = 1.0_rk   !< scale factor for pressure
  real(rk), protected :: e_0 = 1.0_rk   !< scale factor for energy (same as pressure)

contains
  subroutine set_scale_factors(time_scale, length_scale, density_scale)
    !< Set the non-dimensional scale factors based on the provided scales. This
    !< allows the user to set these protected factors outside of the module.

    real(rk), intent(in) :: time_scale
    real(rk), intent(in) :: length_scale
    real(rk), intent(in) :: density_scale

    if(time_scale < 0.0_rk) then
      error stop "Error in mod_nondimensionalization::set_scale_factors(), time_scale < 0"
    end if

    if(length_scale < 0.0_rk) then
      error stop "Error in mod_nondimensionalization::set_scale_factors(), length_scale < 0"
    end if

    if(density_scale < 0.0_rk) then
      error stop "Error in mod_nondimensionalization::set_scale_factors(), density_scale < 0"
    end if

    t_0 = time_scale
    rho_0 = density_scale
    l_0 = length_scale
    v_0 = l_0 / t_0
    p_0 = rho_0 * v_0**2
    e_0 = p_0

    print *
    write(*, '(a)') "Scale factors (CATO is non-dimensional)"
    write(*, '(a)') "========================================"
    write(*, '(a, es10.3)') "Time scale factor (t_0):      ", t_0
    write(*, '(a, es10.3)') "Density scale factor (rho_0): ", rho_0
    write(*, '(a, es10.3)') "Length scale factor (l_0):    ", l_0
    write(*, '(a, es10.3)') "Pressure scale factor (p_0):  ", p_0
    write(*, '(a, es10.3)') "Velocity scale factor (v_0):  ", v_0
    write(*, '(a, es10.3)') "Energy scale factor (e_0):    ", e_0
    write(*, '(a)') "========================================"
    print *

  end subroutine set_scale_factors
end module mod_nondimensionalization
