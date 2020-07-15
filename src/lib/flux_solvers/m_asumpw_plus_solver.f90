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

module mod_mausmpw_plus_solver
  !< Summary: Provide a solver based on the M-AUSMPW+ family of schemes
  !< Date: 07/15/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !<     [1] K.H. Kim, C. Kim, "Accurate, efficient and monotonic numerical methods for multi-dimensional compressible flows Part I: Spatial discretization"
  !<         Journal of Computational Physics 208 (2005) 570–615, https://doi.org/10.1016/j.jcp.2005.02.021
  !<     [2] K.H. Kim, C. Kim, "Accurate, efficient and monotonic numerical methods for multi-dimensional compressible flows Part II: Multi-dimensional limiting process"
  !<         Journal of Computational Physics 208 (2005) 570–615, https://doi.org/10.1016/j.jcp.2005.02.022

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use mod_globals, only: debug_print
  use mod_floating_point_utils, only: neumaier_sum_4
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_edge_interpolator_factory, only: edge_interpolator_factory
  use mod_edge_interpolator, only: edge_iterpolator_t
  use mod_flux_solver, only: flux_solver_t
  use mod_eos, only: eos
  use mod_grid, only: grid_t
  use mod_input, only: input_t

  implicit none

  private
  public :: m_ausmpw_plus_solver_t

  type, extends(flux_solver_t) :: m_ausmpw_plus_solver_t
    !< Implementation of the M-AUSMPW+ scheme
    private
  contains
    ! Public methods
    procedure, public :: initialize => initialize_m_ausmpw_plus
    procedure, public :: solve => solve_m_ausmpw_plus
    procedure, public, pass(lhs) :: copy => copy_m_ausmpw_plus

    ! Private methods
    procedure, private :: flux_edges
    ! procedure, private :: interface_state
    final :: finalize

    ! Operators
    generic :: assignment(=) => copy
  end type m_ausmpw_plus_solver_t
contains
  subroutine initialize_m_ausmpw_plus(self, input)
    !< Constructor for the M-AUSMPW+ solver
    class(m_ausmpw_plus_solver_t), intent(inout) :: self
    class(input_t), intent(in) :: input
  end subroutine initialize_m_ausmpw_plus

  subroutine copy_m_ausmpw_plus(lhs, rhs)
    !< Implement LHS = RHS
    class(m_ausmpw_plus_solver_t), intent(inout) :: lhs
    type(m_ausmpw_plus_solver_t), intent(in) :: rhs

    call debug_print('Running copy_m_ausmpw_plus%copy()', __FILE__, __LINE__)

    lhs%iteration = rhs%iteration
    lhs%time = rhs%time
    lhs%dt = rhs%dt

  end subroutine copy_m_ausmpw_plus

  subroutine solve_m_ausmpw_plus(self, dt, grid, lbounds, rho, u, v, p, d_rho_dt, d_rho_u_dt, d_rho_v_dt, d_rho_E_dt)
    !< Solve and flux the edges
    class(m_ausmpw_plus_solver_t), intent(inout) :: self
    class(grid_t), intent(in) :: grid
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), intent(in) :: dt !< timestep (not really used in this solver, but needed for others)
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(inout) :: rho !< density
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(inout) :: u   !< x-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(inout) :: v   !< y-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(inout) :: p   !< pressure

    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(out) ::   d_rho_dt    !< d/dt of the density field
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(out) :: d_rho_u_dt    !< d/dt of the rhou field
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(out) :: d_rho_v_dt    !< d/dt of the rhov field
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(out) :: d_rho_E_dt    !< d/dt of the rhoE field

  end subroutine solve_m_ausmpw_plus

  subroutine flux_edges(self, grid, lbounds, d_rho_dt, d_rho_u_dt, d_rho_v_dt, d_rho_E_dt)
    !< Flux the edges to get the residuals, e.g. 1/vol * d/dt U
    class(m_ausmpw_plus_solver_t), intent(in) :: self
    class(grid_t), intent(in) :: grid
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) ::   d_rho_dt    !< d/dt of the density field
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: d_rho_u_dt    !< d/dt of the rhou field
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: d_rho_v_dt    !< d/dt of the rhov field
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: d_rho_E_dt    !< d/dt of the rhoE field
  end subroutine flux_edges

  subroutine finalize(self)
    !< Cleanup the M-AUSMPW+ solver
    type(m_ausmpw_plus_solver_t), intent(inout) :: self
  end subroutine finalize

  pure real(rk) function pressure_split_plus(M) result(P_plus)
    !< The pressure splitting function (Eq. 4 in Ref[1]). This is the P+ version. This is kept
    !< simple for inlining
    real(rk), intent(in) :: M !< interface Mach number

    if(abs(M) > 1.0_rk) then
      P_plus = 0.5_rk * (1.0_rk + sign(1.0_rk, M))
    else ! |M| <= 1
      P_plus = 0.25_rk * (2.0_rk - M) * (M + 1.0_rk)**2
    endif
  end function pressure_split_plus

  pure real(rk) function pressure_split_minus(M) result(P_minus)
    !< The pressure splitting function (Eq. 4 in Ref[1]). This is the P- version. This is kept
    !< simple for inlining
    real(rk), intent(in) :: M !< interface Mach number

    if(abs(M) > 1.0_rk) then
      P_minus = 0.5_rk * (1.0_rk - sign(1.0_rk, M))
    else ! |M| <= 1
      P_minus = 0.25_rk * (2.0_rk + M) * (M - 1.0_rk)**2
    endif
  end function pressure_split_minus

  pure real(rk) function mach_split_plus(M) result(M_plus)
    !< The Mach splitting function (Eq. 3 in Ref[1]). This is the M+ version. This is kept
    !< simple for inlining
    real(rk), intent(in) :: M !< interface Mach number

    if(abs(M) > 1.0_rk) then
      M_plus = 0.5_rk * (M + abs(M))
    else ! |M| <= 1
      M_plus = 0.25_rk * (M + 1.0_rk)**2
    endif
  end function mach_split_plus

  pure real(rk) function mach_split_minus(M) result(M_minus)
    !< The Mach splitting function (Eq. 3 in Ref[1]). This is the M- version. This is kept
    !< simple for inlining
    real(rk), intent(in) :: M !< interface Mach number

    if(abs(M) > 1.0_rk) then
      M_minus = 0.5_rk * (M - abs(M))
    else ! |M| <= 1
      M_minus = -0.25_rk * (M - 1.0_rk)**2
    endif
  end function mach_split_minus

  pure real(rk) function smoothness(plus, current, minus) result(r)
    real(rk), intent(in) :: plus, current, minus
    real(rk) :: delta_plus, delta_minus
    delta_plus = plus - current
    if(abs(delta_plus) < EPS) delta_plus = EPS

    delta_minus = current - minus
    if(abs(delta_minus) < EPS) delta_minus = EPS

    r = (delta_minus + EPS) / (delta_plus + EPS)
  end function smoothness

  elemental real(rk) function delta(a, b)
    !< Find the delta in the solution, e.g. delta = a - b. This checks for numbers
    !< near 0 and when a and b are very similar in magnitude. The aim is to avoid
    !< catastrophic cancellation and very small numbers that are essentially 0 for this scenario

    real(rk), intent(in) :: a, b
    real(rk), parameter :: rel_tol = 1e-12_rk     !< relative error tolerance
    real(rk), parameter :: abs_tol = tiny(1.0_rk) !< absolute error tolerance
    real(rk) :: abs_err !< absolute error

    delta = a - b
    abs_err = abs_tol + rel_tol * max(abs(a), abs(b))

    if(abs(delta) < epsilon(1.0_rk)) then
      delta = 0.0_rk
      return
    end if

    if(abs(a) < tiny(1.0_rk) .and. abs(b) < tiny(1.0_rk)) then
      delta = 0.0_rk
      return
    end if

    if(abs(delta) < abs_err) then
      delta = 0.0_rk
      return
    end if
  end function delta

  pure real(rk) function superbee(R) result(psi_lim)
    !< Superbee flux limiter. While this is also defined in the flux limiter module, it
    !< is a key part of M-AUSMPW+, so its here again and inlined for speed.
    real(rk), intent(in) :: R !< smoothness indicator
    real(rk), parameter :: PHI_EPS = 1e-5_rk

    if(R < 0.0_rk .or. abs(R) < PHI_EPS) then
      psi_lim = 0.0_rk
    else if(R > 2.0_rk) then
      psi_lim = 2.0_rk
    else
      psi_lim = max(0.0_rk, min(2.0_rk * R, 1.0_rk), min(R, 2.0_rk))
    end if
  end function superbee

  pure real(rk) function f(p, p_s, w_2)
    !< Implementation of the f function (Eq. 33 in Ref [1]). This is used to help determine
    !< where shock discontinuities are
    real(rk), intent(in) :: p   !< L/R pressure
    real(rk), intent(in) :: p_s !< p_s = p_L * P_L(+) + p_R * P_R(-)
    real(rk), intent(in) :: w_2 !< discontinuity sensor

    if(abs(p_s) < tiny(1.0_rk)) then ! p_s == 0
      f = 0.0_rk
    else
      f = ((p / p_s) - 1.0_rk) * (1.0_rk - w_2)
    end if
  end function f

  pure real(rk) function w_1(p_L, p_R)
    !< Discontinuity sensor w_1 (Eq. 33b in Ref [1]). This detects whether a shock
    !< exists in the normal direction to the cell-interface or not
    real(rk), intent(in) :: p_L !< left interface pressure
    real(rk), intent(in) :: p_R !< right interface pressure

    w_1 = 1.0_rk - min((p_L / p_R),(p_R / p_L))**3
  end function w_1

  pure real(rk) function w_2_xi(p, i, j, lbounds)
    !< Discontinuity sensor w_2 in the xi (i) direction. Equation 31a in Ref [1]
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: p
    integer(ik), intent(in) :: i, j

    w_2_xi = (1.0_rk - min(1.0_rk, &
                           (p(i + 1, j) - p(i, j)) / (0.25_rk * (p(i + 1, j + 1) + p(i, j + 1) - p(i + 1, j - 1) - p(i, j - 1)))) &
              )**2 * &
             (1.0_rk - min(p(i, j) / p(i + 1, j), p(i + 1, j) / p(i, j)))**2
  end function w_2_xi

  pure real(rk) function w_2_eta(p, i, j, lbounds)
    !< Discontinuity sensor w_2 in the eta (j) direction. Equation 31b in Ref [1]
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: p
    integer(ik), intent(in) :: i, j

    w_2_eta = (1.0_rk - min(1.0_rk, &
                            (p(i, j + 1) - p(i, j)) / (0.25_rk * (p(i + 1, j + 1) + p(i + 1, j) - p(i - 1, j + 1) - p(i - 1, j)))) &
               )**2 * &
              (1.0_rk - min(p(i, j) / p(i, j + 1), p(i, j + 1) / p(i, j)))**2
  end function w_2_eta

  pure real(rk) function phi_L_half(phi_L, phi_R, a)
    real(rk), intent(in) :: phi_L !<
    real(rk), intent(in) :: phi_R !<
    real(rk), intent(in) :: a     !< a = 1 - min(1, max(|M_L|, |M_R|))^2

    if(abs(a) < tiny(1.0_rk)) then
      phi_L_half = phi_L
    else
      phi_R_minus_L = phi_L - phi_R
      if(abs(phi_R_minus_L) < epsilon(1.0_rk)) phi_R_minus_L = 0.0_rk

      phi_superbee_minus_L = phi_L_superbee - phi_L
      if(abs(phi_superbee_minus_L) < epsilon(1.0_rk)) phi_superbee_minus_L = 0.0_rk

      phi_L_half = phi_L + (max(0.0_rk, phi_R_minus_L * phi_superbee_minus_L) / (phi_R_minus_L*))
    end if
  end function phi_L_half

  pure real(rk) function phi_R_half()

    phi_L_minus_R = phi_R - phi_L
    if(abs(phi_R_minus_L) < epsilon(1.0_rk)) phi_R_minus_L = 0.0_rk

    phi_superbee_minus_R = phi_R_superbee - phi_R
    if(abs(phi_superbee_minus_R) < epsilon(1.0_rk)) phi_superbee_minus_R = 0.0_rk
  end function phi_R_half

end module mod_mausmpw_plus_solver
