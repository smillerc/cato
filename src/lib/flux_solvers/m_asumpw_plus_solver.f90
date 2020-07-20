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
#define __INTERP_ALIGN__ aligned(rho, u, v, p, rho_sb, u_sb, v_sb, p_sb, edge_flux:__ALIGNBYTES__)
#else
#define __INTERP_ALIGN__
#endif

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
  !<     [3] K.H. Kim, C. Kim, O.H. Rho, "Methods for the Accurate Computations of Hypersonic Flows I. AUSMPW+ Scheme"
  !<         Journal of Computational Physics 174, (2001) 38–80, https://doi.org/10.1006/jcph.2001.6873

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
    real(rk), dimension(:, :, :), allocatable :: i_edge_flux !< ((1:4), i, j) edge flux of the i-direction edges
    real(rk), dimension(:, :, :), allocatable :: j_edge_flux !< ((1:4), i, j) edge flux of the j-direction edges
    real(rk) :: gamma = 0.0_rk

  contains
    ! Public methods
    procedure, public :: initialize => initialize_m_ausmpw_plus
    procedure, public :: solve => solve_m_ausmpw_plus
    procedure, public, pass(lhs) :: copy => copy_m_ausmpw_plus

    ! Private methods
    procedure, private :: get_edge_flux
    procedure, private :: flux_edges
    final :: finalize

    ! Operators
    generic :: assignment(=) => copy
  end type m_ausmpw_plus_solver_t
contains
  subroutine initialize_m_ausmpw_plus(self, input)
    !< Constructor for the M-AUSMPW+ solver
    class(m_ausmpw_plus_solver_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    self%gamma = eos%get_gamma()
  end subroutine initialize_m_ausmpw_plus

  subroutine copy_m_ausmpw_plus(lhs, rhs)
    !< Implement LHS = RHS
    class(m_ausmpw_plus_solver_t), intent(inout) :: lhs
    type(m_ausmpw_plus_solver_t), intent(in) :: rhs

    call debug_print('Running copy_m_ausmpw_plus%copy()', __FILE__, __LINE__)

    lhs%iteration = rhs%iteration
    lhs%time = rhs%time
    lhs%dt = rhs%dt
    lhs%gamma = rhs%gamma

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
    real(rk) :: gamma

    real(rk), dimension(:, :, :), allocatable :: rho_interface_values !< interface (L/R state) values for density
    real(rk), dimension(:, :, :), allocatable :: u_interface_values   !< interface (L/R state) values for x-velocity
    real(rk), dimension(:, :, :), allocatable :: v_interface_values   !< interface (L/R state) values for y-velocity
    real(rk), dimension(:, :, :), allocatable :: p_interface_values   !< interface (L/R state) values for pressure

    real(rk), dimension(:, :, :), allocatable :: rho_superbee_values !< interface (L/R state) values for density
    real(rk), dimension(:, :, :), allocatable :: u_superbee_values   !< interface (L/R state) values for x-velocity
    real(rk), dimension(:, :, :), allocatable :: v_superbee_values   !< interface (L/R state) values for y-velocity
    real(rk), dimension(:, :, :), allocatable :: p_superbee_values   !< interface (L/R state) values for pressure

    real(rk), dimension(:, :), allocatable :: w_2_xi
    real(rk), dimension(:, :), allocatable :: w_2_eta

    ! FIXME:
    allocate(w_2_xi(ilo:ihi, jlo:jhi))
    allocate(w_2_eta(ilo:ihi, jlo:jhi))

    do j = jlo, jhi
      do i = ilo, ihi
        w_2_xi(i, j) = (1.0_rk - min(1.0_rk,(p(i + 1, j) - p(i, j)) / &
                                     (0.25_rk * (p(i + 1, j + 1) + p(i, j + 1) - p(i + 1, j - 1) - p(i, j - 1)))))**2 * &
                       (1.0_rk - min(p(i, j) / p(i + 1, j), p(i + 1, j) / p(i, j)))**2

        w_2_eta(i, j) = (1.0_rk - min(1.0_rk,(p(i, j + 1) - p(i, j)) / &
                                      (0.25_rk * (p(i + 1, j + 1) + p(i + 1, j) - p(i - 1, j + 1) - p(i - 1, j)))))**2 * &
                        (1.0_rk - min(p(i, j) / p(i, j + 1), p(i, j + 1) / p(i, j)))**2
      end do
    end do

    allocate(self%i_edge_flux(4, ilo - 1:ihi, jlo:jhi))
    call self%get_edge_flux(direction='i', grid=grid, &
                            rho=rho_interface_values, u=u_interface_values, &
                            v=v_interface_values, p=p_interface_values, &
                            rho_sb=rho_superbee_values, u_sb=u_superbee_values, &
                            v_sb=v_superbee_values, p_sb=p_superbee_values, &
                            w2=w_2_xi, &
                            bounds=[ilo, ihi, jlo - 1, jhi], edge_flux=self%i_edge_flux)

    allocate(self%j_edge_flux(4, ilo:ihi, jlo - 1:jhi))
    call self%get_edge_flux(direction='j', grid=grid, &
                            rho=rho_interface_values, u=u_interface_values, &
                            v=v_interface_values, p=p_interface_values, &
                            rho_sb=rho_superbee_values, u_sb=u_superbee_values, &
                            v_sb=v_superbee_values, p_sb=p_superbee_values, &
                            w2=w_2_eta, &
                            bounds=[ilo - 1, ihi, jlo, jhi], edge_flux=self%j_edge_flux)

  end subroutine solve_m_ausmpw_plus

  subroutine get_edge_flux(self, direction, grid, rho, u, v, p, rho_sb, u_sb, v_sb, p_sb, w2, bounds, edge_flux)
    !< Solve and flux the edges
    class(m_ausmpw_plus_solver_t), intent(inout) :: self
    class(grid_t), intent(in) :: grid
    integer(ik), dimension(4), intent(in) :: bounds !< (ilo, ihi, jlo, jhi)
    character(len=1), intent(in) :: direction !> 'i', or 'j'

    real(rk), dimension(:, :, :), contiguous, intent(in) :: rho    !< (1:4, i,j); interpolated w/limiter of choice (L/R state) values for density
    real(rk), dimension(:, :, :), contiguous, intent(in) :: u      !< (1:4, i,j); interpolated w/limiter of choice (L/R state) values for x-velocity
    real(rk), dimension(:, :, :), contiguous, intent(in) :: v      !< (1:4, i,j); interpolated w/limiter of choice (L/R state) values for y-velocity
    real(rk), dimension(:, :, :), contiguous, intent(in) :: p      !< (1:4, i,j); interpolated w/limiter of choice (L/R state) values for pressure
    real(rk), dimension(:, :, :), contiguous, intent(in) :: rho_sb !< (1:4, i,j); interpolated w/superbee (L/R state) values for density
    real(rk), dimension(:, :, :), contiguous, intent(in) :: u_sb   !< (1:4, i,j); interpolated w/superbee (L/R state) values for x-velocity
    real(rk), dimension(:, :, :), contiguous, intent(in) :: v_sb   !< (1:4, i,j); interpolated w/superbee (L/R state) values for y-velocity
    real(rk), dimension(:, :, :), contiguous, intent(in) :: p_sb   !< (1:4, i,j); interpolated w/superbee (L/R state) values for pressure
    real(rk), dimension(:, :), contiguous, intent(in) :: w2 !< shock sensing function (determines if shock exists in transversal direction to the interface)
    real(rk), dimension(:, :, :), allocatable, intent(out) :: edge_flux

    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: ilo_bc, ihi_bc, jlo_bc, jhi_bc
    real(rk) :: gamma      !< polytropic gas index
    real(rk) :: rho_L      !< density left state w/ limiter of choice
    real(rk) :: rho_R      !< density right state w/ limiter of choice
    real(rk) :: u_LHS      !< x-velocity left state w/ limiter of choice
    real(rk) :: u_RHS      !< x-velocity right state w/ limiter of choice
    real(rk) :: v_LHS      !< y-velocity left state w/ limiter of choice
    real(rk) :: v_RHS      !< y-velocity right state w/ limiter of choice
    real(rk) :: p_L        !< pressure left state w/ limiter of choice
    real(rk) :: p_R        !< pressure right state w/ limiter of choice
    real(rk) :: rho_L_sb   !< density left state w/ limiter of choice
    real(rk) :: rho_R_sb   !< density right state w/ limiter of choice
    real(rk) :: u_LHS_sb   !< x-velocity left state w/ limiter of choice
    real(rk) :: u_RHS_sb   !< x-velocity right state w/ limiter of choice
    real(rk) :: v_LHS_sb   !< y-velocity left state w/ limiter of choice
    real(rk) :: v_RHS_sb   !< y-velocity right state w/ limiter of choice
    real(rk) :: p_L_sb     !< pressure left state w/ limiter of choice
    real(rk) :: p_R_sb     !< pressure right state w/ limiter of choice
    real(rk) :: p_s        !< p_s = p_L * P_L_plus + p_R * P_R_minus
    real(rk) :: U_L        !< left state velocity component normal to the interface
    real(rk) :: U_R        !< right state velocity component normal to the interface
    real(rk) :: V_L        !< left state velocity component parallel to the interface
    real(rk) :: V_R        !< right state velocity component parallel to the interface
    real(rk) :: rho_L_half !< final interface density left state
    real(rk) :: rho_R_half !< final interface density right state
    real(rk) :: u_L_half   !< final interface x-velocity left state
    real(rk) :: u_R_half   !< final interface x-velocity right state
    real(rk) :: v_L_half   !< final interface y-velocity left state
    real(rk) :: v_R_half   !< final interface y-velocity right state
    real(rk) :: p_L_half   !< final interface pressure left state
    real(rk) :: p_R_half   !< final interface pressure right state
    real(rk) :: H_L_half   !< final interface enthalpy left state
    real(rk) :: H_R_half   !< final interface enthalpy right state
    real(rk) :: H_normal   !< total enthalpy in the normal direction to the interface
    real(rk) :: M_L        !< initial Mach number left state
    real(rk) :: M_R        !< initial Mach number right state
    real(rk) :: H_L        !< initial total enthalpy left state
    real(rk) :: H_R        !< initial total enthalpy right state
    real(rk) :: f_L        !< left state shock sensing function (Eq 33 in Ref [1])
    real(rk) :: f_R        !< right state shock sensing function (Eq 33 in Ref [1])
    real(rk) :: w1         !< pressure sensing function: this detects if there is a shock in the normal direction to the interface
    real(rk) :: w          !< shock sensing function: max(w1, w2); w1: normal direction, w2: transverse direction
    real(rk) :: c_s        !< transversal interface sound speed
    real(rk) :: c_half     !< final interface sound speed
    real(rk) :: a               !< supersonic sensor
    real(rk) :: M_bar_L_plus    !< left final split Mach
    real(rk) :: M_bar_R_minus   !< right final split Mach
    real(rk) :: P_L_plus        !< left split pressure function
    real(rk) :: P_R_minus       !< right split pressure function
    real(rk) :: M_L_plus        !< left split Mach function
    real(rk) :: M_R_minus       !< right split Mach function
    real(rk) :: n_x             !< normal vectors of each face
    real(rk) :: n_y             !< normal vectors of each face
    real(rk), dimension(4) :: psi_L_half ! Convective vectors, Psi
    real(rk), dimension(4) :: psi_R_half ! Convective vectors, Psi
    real(rk), dimension(4) :: p_L_vec! Pressure vectors
    real(rk), dimension(4) :: p_R_vec! Pressure vectors

    integer(ik), parameter :: BOTTOM_IDX = 1 !< edge index for the bottom edge of the current cell
    integer(ik), parameter ::  RIGHT_IDX = 2 !< edge index for the right edge of the current cell
    integer(ik), parameter ::    TOP_IDX = 3 !< edge index for the top edge of the current cell
    integer(ik), parameter ::   LEFT_IDX = 4 !< edge index for the left edge of the current cell

    integer(ik) :: left, right

    !dir$ assume_aligned rho: __ALIGNBYTES__
    !dir$ assume_aligned u: __ALIGNBYTES__
    !dir$ assume_aligned v: __ALIGNBYTES__
    !dir$ assume_aligned p: __ALIGNBYTES__
    !dir$ assume_aligned rho_sb: __ALIGNBYTES__
    !dir$ assume_aligned u_sb: __ALIGNBYTES__
    !dir$ assume_aligned v_sb: __ALIGNBYTES__
    !dir$ assume_aligned p_sb: __ALIGNBYTES__
    !dir$ assume_aligned w2: __ALIGNBYTES__
    !dir$ assume_aligned edge_flux: __ALIGNBYTES__

    select case(direction)
    case('i')
      left = LEFT_IDX
      right = RIGHT_IDX
    case('j')
      left = TOP_IDX
      right = BOTTOM_IDX
    case default
      error stop "Invalid direction in m_ausmpw_plus_solver_t%get_edge_flux(); must be 'i' or 'j'"
    end select

    gamma = self%gamma

    allocate(edge_flux(4, ilo:ihi, jlo:jhi))

    !$omp parallel default(none), &
    !$omp firstprivate(gamma, left, right, ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp private(rho_L, rho_R, u_LHS, u_RHS, v_LHS, v_RHS, p_L, p_R) &
    !$omp private(rho_L_sb, rho_R_sb, u_LHS_sb, u_RHS_sb, v_LHS_sb, v_RHS_sb, p_L_sb, p_R_sb) &
    !$omp private(n_x, n_y, p_s, U_L, U_R, V_L, V_R) &
    !$omp private(rho_L_half, rho_R_half ,u_L_half, u_R_half, v_L_half, v_R_half, p_L_half, p_R_half) &
    !$omp private(H_L_half, H_R_half, H_normal, M_L, M_R, H_L, H_R, f_L, f_R, w1, w,c_s ,c_half, a) &
    !$omp private(M_bar_L_plus, M_bar_R_minus, P_L_plus, P_R_minus, M_L_plus, M_R_minus) &
    !$omp private(psi_L_half, psi_R_half, p_L_vec, p_R_vec) &
    !$omp shared(grid, w2, rho, u, v, p, rho_sb, u_sb, v_sb, p_sb, edge_flux)
    !$omp do
    do j = jlo, jhi
      !$omp simd __INTERP_ALIGN__
      !dir$ vector aligned
      do i = ilo, ihi
        ! use the normal vector of the right edge of the current cell
        n_x = grid%cell_edge_norm_vectors(1, right, i, j)
        n_y = grid%cell_edge_norm_vectors(2, right, i, j)

        ! right edge (i + 1/2)
        rho_L = rho(right, i, j)    ! right edge of the current cell
        rho_R = rho(left, i + 1, j) ! left edge of the cell to the right
        u_LHS = u(right, i, j)      ! right edge of the current cell
        u_RHS = u(left, i + 1, j)   ! left edge of the cell to the right
        v_LHS = v(right, i, j)      ! right edge of the current cell
        v_RHS = v(left, i + 1, j)   ! left edge of the cell to the right
        p_L = p(right, i, j)        ! right edge of the current cell
        p_R = p(left, i + 1, j)     ! left edge of the cell to the right

        rho_L_sb = rho_sb(right, i, j)    ! right edge of the current cell
        rho_R_sb = rho_sb(left, i + 1, j) ! left edge of the cell to the right
        u_LHS_sb = u_sb(right, i, j)        ! right edge of the current cell
        u_RHS_sb = u_sb(left, i + 1, j)     ! left edge of the cell to the right
        v_LHS_sb = v_sb(right, i, j)        ! right edge of the current cell
        v_RHS_sb = v_sb(left, i + 1, j)     ! left edge of the cell to the right
        p_L_sb = p_sb(right, i, j)        ! right edge of the current cell
        p_R_sb = p_sb(left, i + 1, j)     ! left edge of the cell to the right

        ! Velocity normal to the edge, see Fig 2 in Ref[3]
        ! _RHS/_LHS is to avoid naming conflicts with _L and _R (slightly different meaning)
        U_L = u_LHS * n_x + v_LHS * n_y
        U_R = u_RHS * n_x + v_RHS * n_y

        ! Velocity component parallel to the edge
        V_L = u_LHS * (-n_y) + v_LHS * n_x
        V_R = u_RHS * (-n_y) + v_RHS * n_x

        H_L = (gamma / (gamma - 1.0_rk)) * (p_L / rho_L) + 0.5_rk * (U_L**2 + V_L**2)
        H_R = (gamma / (gamma - 1.0_rk)) * (p_R / rho_R) + 0.5_rk * (U_R**2 + V_R**2)

        ! Total enthalpy normal to the edge
        H_normal = min(H_L - 0.5_rk * V_L**2, H_R - 0.5_rk * V_R**2)

        ! Speed of sound normal to the edge
        c_s = sqrt(2.0_rk * ((gamma - 1.0_rk) / (gamma + 1.0_rk)) * H_normal)

        ! Interface sound speed
        if(0.5_rk * (U_L + U_R) < 0.0_rk) then
          c_half = c_s**2 / max(abs(U_L), c_s)
        else
          c_half = c_s**2 / max(abs(U_R), c_s)
        end if

        ! Left/Right Mach number
        M_L = U_L / c_half
        M_R = U_R / c_half
        a = 1.0_rk - min(1.0_rk, max(abs(M_L), abs(M_R)))**2

        ! Mach splitting functions
        M_L_plus = mach_split_plus(M_L)
        M_R_minus = mach_split_minus(M_R)

        ! Pressure splitting functions
        P_L_plus = pressure_split_plus(M_L)
        P_R_minus = pressure_split_minus(M_R)

        p_s = p_L * P_L_plus + p_R * P_R_minus
        w1 = w_1(p_L=p_L, p_R=p_R)
        w = max(w1, w2(i, j))

        f_L = f(p=p_L, p_s=p_s, w_2=w2(i, j))
        f_R = f(p=p_R, p_s=p_s, w_2=w2(i, j))

        if(M_L_plus + M_R_minus < 0.0_rk) then
          M_bar_L_plus = M_L_plus + M_R_minus * ((1.0_rk - w) * (1.0_rk + f_R) - f_L)
          M_bar_R_minus = M_R_minus * w * (1.0_rk + f_R)
        else
          M_bar_L_plus = M_L_plus * w * (1.0_rk + f_L)
          M_bar_R_minus = M_R_minus + M_L_plus * ((1.0_rk - w) * (1.0_rk + f_L) - f_R)
        end if

        rho_L_half = phi_L_half(rho_L, rho_R, rho_L_sb, a)
        rho_R_half = phi_R_half(rho_L, rho_R, rho_R_sb, a)

        u_L_half = phi_L_half(u_LHS, u_RHS, u_LHS_sb, a)
        u_R_half = phi_R_half(u_LHS, u_RHS, u_RHS_sb, a)

        v_L_half = phi_L_half(v_LHS, v_RHS, v_LHS_sb, a)
        v_R_half = phi_R_half(v_LHS, v_RHS, v_RHS_sb, a)

        p_L_half = phi_L_half(p_L, p_R, p_L_sb, a)
        p_R_half = phi_R_half(p_L, p_R, p_R_sb, a)

        ! Enthalpy
        H_L_half = (gamma / (gamma - 1.0_rk)) * (p_L_half / rho_L_half) + 0.5_rk * (u_L_half**2 + v_L_half**2)
        H_R_half = (gamma / (gamma - 1.0_rk)) * (p_R_half / rho_R_half) + 0.5_rk * (u_R_half**2 + v_R_half**2)

        ! Convective vectors, Psi
        psi_L_half = rho_L_half*[1.0_rk, u_L_half, v_L_half, H_L_half]
        psi_R_half = rho_R_half*[1.0_rk, u_R_half, v_R_half, H_R_half]

        ! Pressure vectors
        p_L_vec = [0.0_rk, n_x * p_L, n_y * p_L, 0.0_rk]
        p_R_vec = [0.0_rk, n_x * p_R, n_y * p_R, 0.0_rk]

        ! F 1/2, See Eq. 1 in Ref [1]
        edge_flux(:, i, j) = (M_bar_L_plus * c_half * psi_L_half) + &
                             (M_bar_R_minus * c_half * psi_R_half) + &
                             (P_L_plus * p_L_vec) + (P_R_minus * p_R_vec)
      end do
    end do
    !$omp end do

    !$omp end parallel

  end subroutine get_edge_flux

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

  pure subroutine project_vector(a, b, b_parallel, b_perpendicular)
    !< Project vector "b" onto vector "a" and get the parallel and perpenticular components
    real(rk), dimension(2), intent(in) :: a                !<(x,y); the vector getting projected
    real(rk), dimension(2), intent(in) :: b                !<(x,y); the vector to project onto
    real(rk), dimension(2), intent(out) :: b_parallel      !<(x,y); vector parallel to a
    real(rk), dimension(2), intent(out) :: b_perpendicular !<(x,y); vector perpendicular to a

    b_parallel = (dot_product(b, a) / dot_product(a, a)) * a
    b_perpendicular = b - b_parallel
  end subroutine project_vector

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

    real(rk) :: min_term

    min_term = min((p_L / p_R),(p_R / p_L))
    w_1 = 1.0_rk - min_term * min_term * min_term
  end function w_1

  ! pure real(rk) function w_2_xi(p, i, j, lbounds)
  !   !< Discontinuity sensor w_2 in the xi (i) direction. Equation 31a in Ref [1]
  !   integer(ik), dimension(2), intent(in) :: lbounds
  !   real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: p
  !   integer(ik), intent(in) :: i, j

  !   w_2_xi = (1.0_rk - min(1.0_rk, &
  !                          (p(i + 1, j) - p(i, j)) / (0.25_rk * (p(i + 1, j + 1) + p(i, j + 1) - p(i + 1, j - 1) - p(i, j - 1)))) &
  !             )**2 * &
  !            (1.0_rk - min(p(i, j) / p(i + 1, j), p(i + 1, j) / p(i, j)))**2
  ! end function w_2_xi

  ! pure real(rk) function w_2_eta(p, i, j, lbounds)
  !   !< Discontinuity sensor w_2 in the eta (j) direction. Equation 31b in Ref [1]
  !   integer(ik), dimension(2), intent(in) :: lbounds
  !   real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: p
  !   integer(ik), intent(in) :: i, j

  !   w_2_eta = (1.0_rk - min(1.0_rk, &
  !                           (p(i, j + 1) - p(i, j)) / (0.25_rk * (p(i + 1, j + 1) + p(i + 1, j) - p(i - 1, j + 1) - p(i - 1, j)))) &
  !              )**2 * &
  !             (1.0_rk - min(p(i, j) / p(i, j + 1), p(i, j + 1) / p(i, j)))**2
  ! end function w_2_eta

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
