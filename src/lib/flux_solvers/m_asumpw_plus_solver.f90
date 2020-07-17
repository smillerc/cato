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

    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: ilo_bc, ihi_bc, jlo_bc, jhi_bc

    ! Use the baseline interpolation
    ! call edge_interpolator%interpolate_edge_values(q=rho, lbounds=lbounds, edge_values=rho_interface_values)
    ! call edge_interpolator%interpolate_edge_values(q=u, lbounds=lbounds, edge_values=u_interface_values)
    ! call edge_interpolator%interpolate_edge_values(q=v, lbounds=lbounds, edge_values=v_interface_values)
    ! call edge_interpolator%interpolate_edge_values(q=p, lbounds=lbounds, edge_values=p_interface_values)

    do j = jlo, jhi
      do i = ilo, ihi

        if(0.5_rk * (U_L + U_R) < 0.0_rk) then
          c_half = c_s**2 / max(abs(U_L), c_s)
        else
          c_half = c_s**2 / max(abs(U_R), c_s)
        end if

        M_L = U_L / c_half
        M_R = U_R / c_half

        c_s = sqrt(2.0_rk * ((gamma - 1.0_rk) / (gamma + 1.0_rk)) * H_normal)
        H_normal = min(H_L - 0.5_rk * V_L**2, H_R - 0.5_rk * V_R**2)

        M_L_plus = mach_split_plus(M_L)
        M_R_minus = mach_split_minus(M_R)

        P_L_plus = pressure_split_plus(M_L)
        P_R_minus = pressure_split_minus(M_R)

        if(m_half < 0.0_rk) then
          M_bar_L_plus = M_L_plus + M_R_minus * ((1.0_rk - w) * (1.0_rk + f_R) - f_L)
          M_bar_R_minus = M_R_minus * w * (1.0_rk + f_R)
        else
          M_bar_L_plus = M_L_plus * w * (1.0_rk + f_L)
          M_bar_R_minus = M_R_minus + M_L_plus * ((1.0_rk - w) * (1.0_rk + f_L) - f_R)
        end if

        m_half = M_plus_L + M_minus_R
      end do
    end do

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

  pure real(rk) function phi_L_half(phi_L, phi_R, phi_L_superbee, a)
    !$omp declare simd(phi_L_half) linear(ref(phi_L, phi_R, phi_L_superbee, a)) simdlen(__ALIGNBYTES__)
    real(rk), intent(in) :: phi_L          !< limited primitive interface LHS variable
    real(rk), intent(in) :: phi_R          !< limited primitive interface RHS variable
    real(rk), intent(in) :: phi_L_superbee !< limited primitive interface LHS variable w/ the superbee limiter
    real(rk), intent(in) :: a              !< supersonic function (a=0 in subsonic); a = 1 - min(1, max(|M_L|, |M_R|))^2

    ! Locals
    real(rk) :: phi_R_minus_L        !< phi_R - phi_L
    real(rk) :: phi_superbee_minus_L !< phi_L_superbee - phi_L

    if(abs(a) < tiny(1.0_rk)) then
      phi_L_half = phi_L
    else
      phi_R_minus_L = phi_R - phi_L
      phi_superbee_minus_L = phi_L_superbee - phi_L

      if(abs(phi_R_minus_L) < epsilon(1.0_rk)) phi_R_minus_L = 0.0_rk
      if(abs(phi_superbee_minus_L) < epsilon(1.0_rk)) phi_superbee_minus_L = 0.0_rk

      ! Eq. 27a in Ref [1]
      phi_L_half = phi_L + (max(0.0_rk, phi_R_minus_L * phi_superbee_minus_L) / &
                            (phi_R_minus_L * abs(phi_superbee_minus_L))) * &
                   min(a * 0.5_rk * abs(phi_R_minus_L), abs(phi_superbee_minus_L))
    end if
  end function phi_L_half

  pure real(rk) function phi_R_half(phi_L, phi_R, phi_R_superbee, a)
    !$omp declare simd(phi_R_half) linear(ref(phi_L, phi_R, phi_R_superbee, a)) simdlen(__ALIGNBYTES__)
    real(rk), intent(in) :: phi_L          !< limited primitive interface LHS variable
    real(rk), intent(in) :: phi_R          !< limited primitive interface RHS variable
    real(rk), intent(in) :: phi_R_superbee !< limited primitive interface RHS variable w/ the superbee limiter
    real(rk), intent(in) :: a              !< supersonic function (a=0 in subsonic); a = 1 - min(1, max(|M_L|, |M_R|))^2

    ! Locals
    real(rk) :: phi_L_minus_R        !< phi_L - phi_R
    real(rk) :: phi_superbee_minus_R !< phi_R_superbee - phi_R

    if(abs(a) < tiny(1.0_rk)) then
      phi_R_half = phi_R
    else
      phi_L_minus_R = phi_L - phi_R
      phi_superbee_minus_R = phi_R_superbee - phi_R

      if(abs(phi_L_minus_R) < epsilon(1.0_rk)) phi_L_minus_R = 0.0_rk
      if(abs(phi_superbee_minus_R) < epsilon(1.0_rk)) phi_superbee_minus_R = 0.0_rk

      ! Eq. 27b in Ref [1]
      phi_R_half = phi_R + (max(0.0_rk, phi_L_minus_R * phi_superbee_minus_R) / &
                            (phi_L_minus_R * abs(phi_superbee_minus_R))) * &
                   min(a * 0.5_rk * abs(phi_L_minus_R), abs(phi_superbee_minus_R))
    end if
  end function phi_R_half

end module mod_mausmpw_plus_solver
