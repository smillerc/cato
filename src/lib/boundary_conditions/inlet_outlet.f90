module mod_inlet_outlet
  !> Summary: Define the supersonic and subsonic inlet/outlet functions
  !> Date: 04/23/2020
  !> Author: Sam Miller
  !> Notes:
  !> References:
  !      [1] Computational Fluid Dynamics: Principles and Applications, Jiri Blazek

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_eos, only: eos, universal_gas_constant

  implicit none
  private
  public :: subsonic_inlet, subsonic_outlet, supersonic_inlet, supersonic_outlet

contains
  pure function subsonic_inlet(domain_prim_vars, boundary_norm, inlet_total_temp, &
                               inlet_total_press, inlet_flow_angle) result(boundary_prim_vars)
    !< Set the subsonic inlet conditions by using the outgoing Riemann invariant and sound speed

    real(rk), dimension(4) :: boundary_prim_vars !< (rho, u, v, p); primitive vars at the boundary
    real(rk), dimension(4), intent(in) :: domain_prim_vars !< (rho, u, v, p); primitive vars w/in the domain (real values)
    real(rk), dimension(2), intent(in) :: boundary_norm  !< normal vector to the boundary edge
    real(rk), intent(in) :: inlet_total_temp  !< specified total temperature (T0) at the inlet
    real(rk), intent(in) :: inlet_total_press !< specified total pressure (p0) at the inlet
    real(rk), intent(in) :: inlet_flow_angle  !< specified flow angle at the inlet; this is with respect to the x-axis

    real(rk) :: riemann_inv !< Riemann invariant
    real(rk) :: cs_d !< sound speed of the domain
    real(rk) :: cs_b !< sound speed at the boundary
    real(rk) :: cs_0_sq !< stagnation sound speed squared
    real(rk) :: gamma !< ideal gas gamma
    real(rk) :: cp !< specific head capactity at constant pressure of the gas
    real(rk) :: cos_theta
    real(rk) :: tot_vel !< total velocity in the domain, e.g. norm2(vel vector)
    real(rk) :: v_dot_n !< domain velocity vector dotted with the normal of the boundary face
    real(rk) :: T_b   !< temperature at the boundary
    real(rk) :: p_b   !< pressure at the boundary
    real(rk) :: rho_b !< density at the boundary
    real(rk) :: u_b   !< x-velocity at the boundary
    real(rk) :: v_b   !< y-velocity at the boundary

    gamma = eos%get_gamma()
    cp = eos%get_cp()
    cos_theta = 0.0_rk

    if(inlet_total_temp < 0.0_rk) error stop "Error in subsonic_inlet, inlet temperature < 0"
    if(inlet_total_press < 0.0_rk) error stop "Error in subsonic_inlet, inlet pressure < 0"

    associate(rho=>domain_prim_vars(1), u=>domain_prim_vars(2), &
              v=>domain_prim_vars(3), p=>domain_prim_vars(4), &
              n_x=>boundary_norm(1), n_y=>boundary_norm(2), &
              T_0=>inlet_total_temp, p_0=>inlet_total_press, &
              theta=>inlet_flow_angle, &
              R_minus=>riemann_inv, R=>universal_gas_constant)

      v_dot_n = u * nx + v * n_y

      tot_vel = sqrt(u**2 + v**2)
      if(tot_vel > 0.0_rk) cos_theta = -v_dot_n / tot_vel

      ! domain sound speed
      cs_d = eos%sound_speed(pressure=p, density=rho)

      ! stagnation sound speed
      cs_0_sq = cs_d**2 + ((gamma - 1.0_rk) / 2.0_rk) * tot_vel**2

      ! outgoing Riemann invariant R-
      R_minus = v_dot_n - (2.0_rk * cs_d / (gamma - 1.0_rk))

      ! boundary sound speed
      cs_b = (-R_minus * (gamma - 1.0_rk) / ((gamma - 1.0_rk) * cos_theta**2 + 2.0_rk)) * &
             (1.0_rk + cos_theta * &
              sqrt(((((gamma - 1.0_rk) * cos_theta**2 + 2.0_rk) * cs_0_sq) / ((gamma - 1.0_rk) * R_minus**2)) - &
                   ((gamma - 1.0_rk) / 2.0_rk)))

      ! boundary state
      T_b = T_0 * (cs_b**2 / cs_0_sq)
      p_b = p_0 * (T_b / T_0)**(gamma / (gamma - 1.0_rk))
      rho_b = p_b / (R * T_b)
      total_v_b = sqrt(2.0_rk * cp * (T_0 - T_b))

      u_b = total_v_b * cos(theta)
      v_b = total_v_b * sin(theta)
    end associate

    boundary_prim_vars = [rho_b, u_b, v_b, p_b]

  end function subsonic_inlet

  pure function subsonic_outlet(domain_prim_vars, exit_pressure) result(boundary_prim_vars)
    !< Set the boundary to be a subsonic outlet. The only input parameter is the exit pressure
    !< Refer to Section 8.4 in Ref [1]

    real(rk), intent(in) :: exit_pressure !< specified pressure at the exit boundary
    real(rk), dimension(4), intent(in) :: domain_prim_vars !< (rho, u, v, p); primitive vars w/in the domain (real values)

    real(rk) :: T_b     !< temperature at the boundary
    real(rk) :: p_b     !< pressure at the boundary
    real(rk) :: rho_b   !< density at the boundary
    real(rk) :: u_b     !< x-velocity at the boundary
    real(rk) :: v_b     !< y-velocity at the boundary
    real(rk) :: cs_d    !< domain sound speed

    p_b = exit_pressure

    associate(rho_d=>domain_prim_vars(1), u_d=>domain_prim_vars(2), &
              v_d=>domain_prim_vars(3), p_d=>domain_prim_vars(4))

      cs_d = eos%sound_speed(pressure=p_d, density=rho_d)

      rho_b = rho_d + (p_b - p_d) / cs_d**2
      u_b = u_d + n_x * (p_d - p_b) / (rho_d * cs_d)
      v_b = v_d + n_y * (p_d - p_b) / (rho_d * cs_d)
    end associate

    boundary_prim_vars = [rho_b, u_b, v_b, p_b]

  end function subsonic_outlet

  pure function supersonic_inlet(inlet_prim_vars) result(boundary_prim_vars)
    real(rk), dimension(4) :: boundary_prim_vars !< (rho, u, v, p); primitive vars at the boundary
    real(rk), dimension(4), intent(in) :: inlet_prim_vars !< (rho, u, v, p); primitive vars specified at the inlet

    boundary_prim_vars = inlet_prim_vars
  end function supersonic_inlet

  pure function supersonic_outlet(domain_prim_vars) result(boundary_prim_vars)
    real(rk), dimension(4) :: boundary_prim_vars !< (rho, u, v, p); primitive vars at the boundary
    real(rk), dimension(4), intent(in) :: domain_prim_vars !< (rho, u, v, p); primitive vars w/in the domain (real values)

    boundary_prim_vars = domain_prim_vars
  end function supersonic_outlet
end module mod_inlet_outlet
