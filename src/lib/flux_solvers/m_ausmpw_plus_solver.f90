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

! Fypp variables. This allows us to generate an edge flux subroutine for each direction
! and still allow the compiler to optimize

#ifdef __SIMD_ALIGN_OMP__
#define __INTERP_ALIGN__ aligned(rho, u, v, p, rho_sb, u_sb, v_sb, p_sb, edge_flux:__ALIGNBYTES__)
#define __W2_ALIGN__ aligned(p, w_2_xi, w_2_eta:__ALIGNBYTES__)
#else
#define __INTERP_ALIGN__
#define __W2_ALIGN__
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
  use mod_flux_solver, only: edge_split_flux_solver_t
  use mod_eos, only: eos
  use mod_grid, only: grid_t
  use mod_input, only: input_t

  implicit none

  private
  public :: m_ausmpw_plus_solver_t

  type, extends(edge_split_flux_solver_t) :: m_ausmpw_plus_solver_t
    !< Implementation of the M-AUSMPW+ scheme
    private
    real(rk) :: gamma = 0.0_rk

  contains
    ! Public methods
    procedure, public :: initialize => initialize_m_ausmpw_plus
    procedure, public :: solve => solve_m_ausmpw_plus
    procedure, public, pass(lhs) :: copy => copy_m_ausmpw_plus
    ! procedure, public :: flux_split_edges
    ! Private methods
    procedure, private :: get_iflux
    procedure, private :: get_jflux
    final :: finalize

    ! Operators
    generic :: assignment(=) => copy
  end type m_ausmpw_plus_solver_t
contains
  subroutine initialize_m_ausmpw_plus(self, input)
    !< Constructor for the M-AUSMPW+ solver
    class(m_ausmpw_plus_solver_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    self%input = input
    self%gamma = eos%get_gamma()
    self%name = 'M-AUSMPW+_'//input%limiter
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

    class(edge_iterpolator_t), pointer :: edge_interpolator => null()
    class(edge_iterpolator_t), pointer :: sb_edge_interpolator => null()
    class(boundary_condition_t), allocatable:: bc_plus_x
    class(boundary_condition_t), allocatable:: bc_plus_y
    class(boundary_condition_t), allocatable:: bc_minus_x
    class(boundary_condition_t), allocatable:: bc_minus_y

    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: ilo_bc, ihi_bc, jlo_bc, jhi_bc

    real(rk), dimension(:, :, :), allocatable :: rho_interface_values !< interface (L/R state) values for density
    real(rk), dimension(:, :, :), allocatable :: u_interface_values   !< interface (L/R state) values for x-velocity
    real(rk), dimension(:, :, :), allocatable :: v_interface_values   !< interface (L/R state) values for y-velocity
    real(rk), dimension(:, :, :), allocatable :: p_interface_values   !< interface (L/R state) values for pressure

    real(rk), dimension(:, :, :), allocatable :: rho_superbee_values !< interface (L/R state) values for density
    real(rk), dimension(:, :, :), allocatable :: u_superbee_values   !< interface (L/R state) values for x-velocity
    real(rk), dimension(:, :, :), allocatable :: v_superbee_values   !< interface (L/R state) values for y-velocity
    real(rk), dimension(:, :, :), allocatable :: p_superbee_values   !< interface (L/R state) values for pressure

    real(rk) :: denom_xi           !< scalar used to aid in the compuation of w_2_xi
    real(rk) :: numer_xi           !< scalar used to aid in the compuation of w_2_xi
    real(rk) :: w_2_xi_first_term  !< scalar used to aid in the compuation of w_2_xi
    real(rk) :: denom_eta          !< scalar used to aid in the compuation of w_2_eta
    real(rk) :: numer_eta          !< scalar used to aid in the compuation of w_2_eta
    real(rk) :: w_2_eta_first_term !< scalar used to aid in the compuation of w_2_eta
    real(rk), dimension(:, :), allocatable :: w_2_xi
    real(rk), dimension(:, :), allocatable :: w_2_eta

    call debug_print('Running m_ausmpw_plus_solver_t%solve_m_ausmpw_plus()', __FILE__, __LINE__)

    if(dt < tiny(1.0_rk)) then
      write(std_err, '(a, es16.6)') "Error in m_ausmpw_plus_solver_t%solve_m_ausmpw_plus(), "// &
        "the timestep dt is < tiny(1.0_rk): dt = ", dt
      error stop "Error in m_ausmpw_plus_solver_t%solve_m_ausmpw_plus(), the timestep dt is < tiny(1.0_rk)"
    end if
    self%time = self%time + dt
    self%dt = dt
    self%iteration = self%iteration + 1

    call self%init_boundary_conditions(grid, &
                                       bc_plus_x=bc_plus_x, bc_minus_x=bc_minus_x, &
                                       bc_plus_y=bc_plus_y, bc_minus_y=bc_minus_y)
    call self%apply_primitive_bc(rho=rho, u=u, v=v, p=p, lbounds=lbounds, &
                                 bc_plus_x=bc_plus_x, bc_minus_x=bc_minus_x, &
                                 bc_plus_y=bc_plus_y, bc_minus_y=bc_minus_y)

    edge_interpolator => edge_interpolator_factory(self%input)
    call edge_interpolator%interpolate_edge_values(q=rho, lbounds=lbounds, edge_values=rho_interface_values)
    call edge_interpolator%interpolate_edge_values(q=u, lbounds=lbounds, edge_values=u_interface_values)
    call edge_interpolator%interpolate_edge_values(q=v, lbounds=lbounds, edge_values=v_interface_values)
    call edge_interpolator%interpolate_edge_values(q=p, lbounds=lbounds, edge_values=p_interface_values)

    sb_edge_interpolator => edge_interpolator_factory(self%input, limiter_name='superbee')
    call sb_edge_interpolator%interpolate_edge_values(q=rho, lbounds=lbounds, edge_values=rho_superbee_values)
    call sb_edge_interpolator%interpolate_edge_values(q=u, lbounds=lbounds, edge_values=u_superbee_values)
    call sb_edge_interpolator%interpolate_edge_values(q=v, lbounds=lbounds, edge_values=v_superbee_values)
    call sb_edge_interpolator%interpolate_edge_values(q=p, lbounds=lbounds, edge_values=p_superbee_values)

    call self%apply_reconstructed_bc(recon_rho=rho_interface_values, recon_u=u_interface_values, &
                                     recon_v=v_interface_values, recon_p=p_interface_values, &
                                     lbounds=lbounds, &
                                     bc_plus_x=bc_plus_x, bc_minus_x=bc_minus_x, &
                                     bc_plus_y=bc_plus_y, bc_minus_y=bc_minus_y)

    call self%apply_reconstructed_bc(recon_rho=rho_superbee_values, recon_u=u_superbee_values, &
                                     recon_v=v_superbee_values, recon_p=p_superbee_values, &
                                     lbounds=lbounds, &
                                     bc_plus_x=bc_plus_x, bc_minus_x=bc_minus_x, &
                                     bc_plus_y=bc_plus_y, bc_minus_y=bc_minus_y)

    !                   edge(3)
    !                  jflux(i,j+1)  'R'
    !               o--------------------o
    !               |                'L' |
    !               |                    |
    !   iflux(i, j) |     cell (i,j)     | iflux(i+1, j)
    !    edge(4)    |                    |  edge(2)
    !               |                'L' | 'R'
    !               o--------------------o
    !                  jflux(i,j)
    !                   edge(1)
    !  For jflux: "R" is the bottom of the cell above, "L" is the top of the current cell
    !  For iflux: "R" is the left of the cell to the right, "L" is the right of the current cell

    ilo = grid%ilo_cell
    ihi = grid%ihi_cell
    jlo = grid%jlo_cell
    jhi = grid%jhi_cell

    if(allocated(self%iflux)) deallocate(self%iflux)
    if(allocated(self%jflux)) deallocate(self%jflux)

    ! allocate(self%iflux(4, ilo:ihi + 1, jlo:jhi))
    ! allocate(self%jflux(4, ilo:ihi, jlo:jhi + 1))
    allocate(self%iflux(4, ilo - 1:ihi, jlo:jhi))
    allocate(self%jflux(4, ilo:ihi, jlo - 1:jhi))

    ilo_bc = grid%ilo_bc_cell
    ihi_bc = grid%ihi_bc_cell
    jlo_bc = grid%jlo_bc_cell
    jhi_bc = grid%jhi_bc_cell

    ! Find the shock sensor values w_2_xi and w_2_eta
    allocate(w_2_xi(ilo_bc:ihi_bc, jlo_bc:jhi_bc))  !dir$ attributes align:__ALIGNBYTES__ :: w_2_xi
    allocate(w_2_eta(ilo_bc:ihi_bc, jlo_bc:jhi_bc)) !dir$ attributes align:__ALIGNBYTES__ :: w_2_eta
    w_2_xi = 0.0_rk
    w_2_eta = 0.0_rk

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j, denom_xi, numer_xi, w_2_xi_first_term, denom_eta, numer_eta, w_2_eta_first_term) &
    !$omp shared(p, w_2_xi, w_2_eta)
    !$omp do
    do j = jlo, jhi
      !$omp simd __W2_ALIGN__
      !dir$ vector aligned
      do i = ilo, ihi
        ! Xi direction (i)
        denom_xi = p(i + 1, j + 1) + p(i, j + 1) - p(i + 1, j - 1) - p(i, j - 1)
        numer_xi = p(i + 1, j) - p(i, j)
        if(abs(numer_xi) < epsilon(1.0_rk)) then
          w_2_xi_first_term = 0.0_rk
        else if(abs(denom_xi) < epsilon(1.0_rk)) then
          w_2_xi_first_term = 1.0_rk
        else
          w_2_xi_first_term = min(1.0_rk, numer_xi / (0.25_rk * denom_xi))
        end if
        w_2_xi(i, j) = (1.0_rk - w_2_xi_first_term)**2 * (1.0_rk - min(p(i, j) / p(i + 1, j), p(i + 1, j) / p(i, j)))**2

        ! Eta direction (j)
        denom_eta = (p(i + 1, j + 1) + p(i + 1, j) - p(i - 1, j + 1) - p(i - 1, j))
        numer_eta = p(i, j + 1) - p(i, j)
        if(abs(numer_eta) < epsilon(1.0_rk)) then
          w_2_eta_first_term = 0.0_rk
        else if(abs(denom_eta) < epsilon(1.0_rk)) then
          w_2_eta_first_term = 1.0_rk
        else
          w_2_eta_first_term = min(1.0_rk, numer_eta / (0.25_rk * denom_eta))
        end if
        w_2_eta(i, j) = (1.0_rk - w_2_eta_first_term)**2 * (1.0_rk - min(p(i, j) / p(i, j + 1), p(i, j + 1) / p(i, j)))**2
      end do
    end do
    !$omp end do
    !$omp end parallel

    call self%get_iflux(grid=grid, &
                        lbounds=lbound(rho_interface_values), &
                        eflux_lbounds=lbound(self%iflux), &
                        rho_ave=rho, u_ave=u, v_ave=v, p_ave=p, &
                        rho=rho_interface_values, u=u_interface_values, &
                        v=v_interface_values, p=p_interface_values, &
                        rho_sb=rho_superbee_values, u_sb=u_superbee_values, &
                        v_sb=v_superbee_values, p_sb=p_superbee_values, &
                        w2=w_2_xi, edge_flux=self%iflux)

    call self%get_jflux(grid=grid, &
                        lbounds=lbound(rho_interface_values), &
                        eflux_lbounds=lbound(self%jflux), &
                        rho_ave=rho, u_ave=u, v_ave=v, p_ave=p, &
                        rho=rho_interface_values, u=u_interface_values, &
                        v=v_interface_values, p=p_interface_values, &
                        rho_sb=rho_superbee_values, u_sb=u_superbee_values, &
                        v_sb=v_superbee_values, p_sb=p_superbee_values, &
                        w2=w_2_eta, edge_flux=self%jflux)

    ! Now flux the edges to get the next solution
    call self%flux_split_edges(grid, lbounds, d_rho_dt, d_rho_u_dt, d_rho_v_dt, d_rho_E_dt)

    deallocate(edge_interpolator)
    deallocate(sb_edge_interpolator)
    deallocate(w_2_eta)
    deallocate(w_2_xi)

    deallocate(rho_interface_values)
    deallocate(u_interface_values)
    deallocate(v_interface_values)
    deallocate(p_interface_values)

    deallocate(rho_superbee_values)
    deallocate(u_superbee_values)
    deallocate(v_superbee_values)
    deallocate(p_superbee_values)

    deallocate(bc_plus_x)
    deallocate(bc_plus_y)
    deallocate(bc_minus_x)
    deallocate(bc_minus_y)

    deallocate(self%iflux)
    deallocate(self%jflux)
  end subroutine solve_m_ausmpw_plus

  subroutine get_iflux(self, grid, lbounds, rho, u, v, p, rho_sb, u_sb, v_sb, p_sb, rho_ave, u_ave, v_ave, p_ave, w2,&
      & eflux_lbounds, edge_flux)
    !< Construct the fluxes for each edge in the i direction. This is templated via the Fypp pre-processor
    class(m_ausmpw_plus_solver_t), intent(inout) :: self
    class(grid_t), intent(in) :: grid
    integer(ik), dimension(3), intent(in) :: lbounds !< bounds of the primitive variable arrays
    real(rk), dimension(lbounds(2):, lbounds(3):), contiguous, intent(in) :: rho_ave    !< (i,j); cell averaged value ofdensity; needed for critical Mach number calcs
    real(rk), dimension(lbounds(2):, lbounds(3):), contiguous, intent(in) :: u_ave    !< (i,j); cell averaged value of x-velocity; needed for critical Mach number calcs
    real(rk), dimension(lbounds(2):, lbounds(3):), contiguous, intent(in) :: v_ave    !< (i,j); cell averaged value of y-velocity; needed for critical Mach number calcs
    real(rk), dimension(lbounds(2):, lbounds(3):), contiguous, intent(in) :: p_ave    !< (i,j); cell averaged value of pressure; needed for critical Mach number calcs
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: rho    !< (1:4, i,j); interpolated w/limiter of choice (L/R state) values for density
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: u      !< (1:4, i,j); interpolated w/limiter of choice (L/R state) values for x-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: v      !< (1:4, i,j); interpolated w/limiter of choice (L/R state) values for y-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: p      !< (1:4, i,j); interpolated w/limiter of choice (L/R state) values for pressure
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: rho_sb !< (1:4, i,j); interpolated w/superbee (L/R state) values for density
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: u_sb   !< (1:4, i,j); interpolated w/superbee (L/R state) values for x-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: v_sb   !< (1:4, i,j); interpolated w/superbee (L/R state) values for y-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: p_sb   !< (1:4, i,j); interpolated w/superbee (L/R state) values for pressure
    real(rk), dimension(lbounds(2):, lbounds(3):), contiguous, intent(in) :: w2 !< shock sensing function (determines if shock exists in transversal direction to the interface)
    integer(ik), dimension(3), intent(in) :: eflux_lbounds !< bounds of the primitive variable arrays
    real(rk), dimension(eflux_lbounds(1):, eflux_lbounds(2):, eflux_lbounds(3):), contiguous, intent(inout) :: edge_flux

    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi
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
    real(rk) :: m_half     !< interface Mach number
    real(rk) :: a               !< supersonic sensor
    real(rk) :: M_bar_L_plus    !< left final split Mach
    real(rk) :: M_bar_R_minus   !< right final split Mach
    real(rk) :: P_L_plus        !< left split pressure function
    real(rk) :: P_R_minus       !< right split pressure function
    real(rk) :: M_L_plus        !< left split Mach function
    real(rk) :: M_R_minus       !< right split Mach function
    real(rk) :: n_x             !< normal vectors of each face
    real(rk) :: n_y             !< normal vectors of each face
    real(rk) :: M_star, M_star_2, M_star_1
    real(rk) :: vel_ave, vel_ave_2
    real(rk) :: mass_flux_L, mass_flux_R

    integer(ik), parameter :: BOTTOM_IDX = 1 !< edge index for the bottom edge of the current cell
    integer(ik), parameter ::  RIGHT_IDX = 2 !< edge index for the right edge of the current cell
    integer(ik), parameter ::    TOP_IDX = 3 !< edge index for the top edge of the current cell
    integer(ik), parameter ::   LEFT_IDX = 4 !< edge index for the left edge of the current cell

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

    gamma = self%gamma

    ilo = lbound(edge_flux, dim=2)
    ihi = ubound(edge_flux, dim=2)
    jlo = lbound(edge_flux, dim=3)
    jhi = ubound(edge_flux, dim=3)

    !                   edge(3)
    !                  jflux(i,j+1)  'R'
    !               o--------------------o
    !               |                'L' |
    !               |                    |
    !   iflux(i, j) |     cell (i,j)     | iflux(i+1, j)
    !    edge(4)    |                    |  edge(2)
    !               |                'L' | 'R'
    !               o--------------------o
    !                  jflux(i,j)
    !                   edge(1)
    !  For jflux: "R" is the bottom of the cell above, "L" is the top of the current cell
    !  For iflux: "R" is the left of the cell to the right, "L" is the right of the current cell

    !$omp parallel default(none), &
    !$omp firstprivate(gamma, ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp private(rho_L, rho_R, u_LHS, u_RHS, v_LHS, v_RHS, p_L, p_R, m_half) &
    !$omp private(rho_L_sb, rho_R_sb, u_LHS_sb, u_RHS_sb, v_LHS_sb, v_RHS_sb, p_L_sb, p_R_sb) &
    !$omp private(n_x, n_y, p_s, U_L, U_R, V_L, V_R) &
    !$omp private(rho_L_half, rho_R_half ,u_L_half, u_R_half, v_L_half, v_R_half, p_L_half, p_R_half) &
    !$omp private(H_L_half, H_R_half, H_normal, M_L, M_R, H_L, H_R, f_L, f_R, w1, w,c_s ,c_half, a) &
    !$omp private(M_bar_L_plus, M_bar_R_minus, P_L_plus, P_R_minus, M_L_plus, M_R_minus) &
    !$omp private(mass_flux_L, mass_flux_R) &
    !$omp private(M_star, M_star_2, M_star_1, vel_ave, vel_ave_2) &
    !$omp shared(rho_ave, u_ave, v_ave, p_ave) &
    !$omp shared(grid, w2, rho, u, v, p, rho_sb, u_sb, v_sb, p_sb, edge_flux)
    !$omp do
    do j = jlo, jhi
      !$omp simd __INTERP_ALIGN__
      !dir$ vector aligned
      do i = ilo, ihi

        n_x = grid%cell_edge_norm_vectors(1, RIGHT_IDX, i, j)
        n_y = grid%cell_edge_norm_vectors(2, RIGHT_IDX, i, j)

        rho_L = rho(RIGHT_IDX, i, j)
        u_LHS = u(RIGHT_IDX, i, j)
        v_LHS = v(RIGHT_IDX, i, j)
        p_L = p(RIGHT_IDX, i, j)
        rho_R = rho(LEFT_IDX, i + 1, j)
        u_RHS = u(LEFT_IDX, i + 1, j)
        v_RHS = v(LEFT_IDX, i + 1, j)
        p_R = p(LEFT_IDX, i + 1, j)

        ! Superbee values
        rho_L_sb = rho_sb(RIGHT_IDX, i, j)
        u_LHS_sb = u_sb(RIGHT_IDX, i, j)
        v_LHS_sb = v_sb(RIGHT_IDX, i, j)
        p_L_sb = p_sb(RIGHT_IDX, i, j)
        rho_R_sb = rho_sb(LEFT_IDX, i + 1, j)
        u_RHS_sb = u_sb(LEFT_IDX, i + 1, j)
        v_RHS_sb = v_sb(LEFT_IDX, i + 1, j)
        p_R_sb = p_sb(LEFT_IDX, i + 1, j)

        U_L = u_LHS
        U_R = u_RHS
        V_L = v_LHS
        V_R = v_RHS

        ! Velocity normal to the edge, see Fig 2 in Ref[3]
        ! _RHS/_LHS is to avoid naming conflicts with _L and _R (slightly different meaning)
        ! U_L = u_LHS * n_x + v_LHS * n_y
        ! U_R = u_RHS * n_x + v_RHS * n_y

        ! ! Velocity component parallel to the edge
        ! V_L = u_LHS * (-n_y) + v_LHS * n_x
        ! V_R = u_RHS * (-n_y) + v_RHS * n_x

        H_L = (gamma / (gamma - 1.0_rk)) * (p_L / rho_L) + 0.5_rk * (u_LHS**2 + v_LHS**2)
        H_R = (gamma / (gamma - 1.0_rk)) * (p_R / rho_R) + 0.5_rk * (u_RHS**2 + v_RHS**2)

        ! Total enthalpy normal to the edge
        ! H_normal = min(H_L - 0.5_rk * V_L**2, H_R - 0.5_rk * V_R**2)
        H_normal = min(H_L - 0.5_rk * V_L**2, H_R - 0.5_rk * V_R**2)

        ! Speed of sound normal to the edge, also like the critical sound speed
        ! across a normal shock
        c_s = sqrt(2.0_rk * ((gamma - 1.0_rk) / (gamma + 1.0_rk)) * H_normal)

        ! Intermediate charactersistic Mach numbers M*
        M_star = sqrt(u_ave(i, j)**2 + v_ave(i, j)**2) / c_s ! M*(i)

        M_star_2 = sqrt(u_ave(i + 1, j)**2 + v_ave(i + 1, j)**2) / c_s ! M*(i+1)

        ! Interface sound speed
        if(0.5_rk * (U_L + U_R) < 0.0_rk) then  ! part (ii) in the paper after Eq 3
          c_half = c_s**2 / max(abs(U_R), c_s)
        else
          c_half = c_s**2 / max(abs(U_L), c_s) ! part (i)
        end if

        ! Left/Right Mach number
        M_L = U_L / c_half
        M_R = U_R / c_half
        a = 1.0_rk - min(1.0_rk, max(abs(M_L), abs(M_R)))**2

        ! Mach splitting functions
        M_L_plus = mach_split_plus(M_L)
        M_R_minus = mach_split_minus(M_R)

        vel_ave = sqrt(u_ave(i, j)**2 + v_ave(i, j)**2)
        vel_ave_2 = sqrt(u_ave(i + 1, j)**2 + v_ave(i + 1, j)**2)

        ! Pressure splitting functions
        ! Eq. 38a in Ref [1]
        if(M_star > 1.0_rk .and. M_star_2 < 1.0_rk .and. 0.0_rk < M_star * M_star_2 .and. M_star * M_star_2 < 1.0_rk) then
          P_R_minus = 1.0_rk - ((rho_ave(i, j) * vel_ave * (vel_ave - vel_ave_2) + p_ave(i, j)) / p_ave(i + 1, j))
        else
          P_R_minus = pressure_split_minus(M_R)
        end if

        ! Eq, 38b in Ref [1]
        if(M_star < -1.0_rk .and. M_star_2 < -1.0_rk .and. 0.0_rk < M_star * M_star_2 .and. M_star * M_star_2 < 1.0_rk) then
          P_L_plus = 1.0_rk - ((rho_ave(i + 1, j) * vel_ave_2 * (vel_ave_2 - vel_ave) + p_ave(i + 1, j)) / p_ave(i, j))
        else
          P_L_plus = pressure_split_plus(M_L)
        end if

        p_s = p_L * P_L_plus + p_R * P_R_minus
        w1 = w_1(p_L=p_L, p_R=p_R)
        w = max(w1, w2(i, j))

        f_L = f(p=p_L, p_s=p_s, w_2=w2(i, j))
        f_R = f(p=p_R, p_s=p_s, w_2=w2(i, j))

        ! Eq. 2b in Ref [1]
        if(M_L_plus + M_R_minus < 0.0_rk) then
          M_bar_L_plus = M_L_plus * w * (1.0_rk + f_L)
          M_bar_R_minus = M_R_minus + M_L_plus * ((1.0_rk - w) * (1.0_rk + f_L) - f_R)
        else ! Eq. 2a in Ref [1]
          M_bar_L_plus = M_L_plus + M_R_minus * ((1.0_rk - w) * (1.0_rk + f_R) - f_L)
          M_bar_R_minus = M_R_minus * w * (1.0_rk + f_R)
        end if

        rho_L_half = phi_L_half(phi_L=rho_L, phi_R=rho_R, phi_L_superbee=rho_L_sb, a=a)
        rho_R_half = phi_R_half(phi_L=rho_L, phi_R=rho_R, phi_R_superbee=rho_R_sb, a=a)

        u_L_half = phi_L_half(phi_L=u_LHS, phi_R=u_RHS, phi_L_superbee=u_LHS_sb, a=a)
        u_R_half = phi_R_half(phi_L=u_LHS, phi_R=u_RHS, phi_R_superbee=u_RHS_sb, a=a)

        v_L_half = phi_L_half(phi_L=v_LHS, phi_R=v_RHS, phi_L_superbee=v_LHS_sb, a=a)
        v_R_half = phi_R_half(phi_L=v_LHS, phi_R=v_RHS, phi_R_superbee=v_RHS_sb, a=a)

        p_L_half = phi_L_half(phi_L=p_L, phi_R=p_R, phi_L_superbee=p_L_sb, a=a)
        p_R_half = phi_R_half(phi_L=p_L, phi_R=p_R, phi_R_superbee=p_R_sb, a=a)

        ! Enthalpy
        H_L_half = (gamma / (gamma - 1.0_rk)) * (p_L_half / rho_L_half) + 0.5_rk * (u_L_half**2 + v_L_half**2)
        H_R_half = (gamma / (gamma - 1.0_rk)) * (p_R_half / rho_R_half) + 0.5_rk * (u_R_half**2 + v_R_half**2)

        mass_flux_L = M_bar_L_plus * c_half * rho_L_half
        mass_flux_R = M_bar_R_minus * c_half * rho_R_half
        edge_flux(1, i, j) = (mass_flux_L) + (mass_flux_R)
       edge_flux(2, i, j) = (mass_flux_L * u_L_half) + (mass_flux_R * u_R_half) + ((P_L_plus * n_x * p_L) + (P_R_minus * n_x * p_R))
       edge_flux(3, i, j) = (mass_flux_L * v_L_half) + (mass_flux_R * v_R_half) + ((P_L_plus * n_y * p_L) + (P_R_minus * n_y * p_R))
        edge_flux(4, i, j) = (mass_flux_L * H_L_half) + (mass_flux_R * H_R_half)

        ! if (i == 51 .or. i == 52) then

        !   ! write(*, '(a, i3, 2(es16.6))') 'BOTTOM_IDX: ', BOTTOM_IDX, grid%cell_edge_norm_vectors(1:2, BOTTOM_IDX, i, j)
        !   ! write(*, '(a, i3, 2(es16.6))') 'RIGHT_IDX : ', RIGHT_IDX , grid%cell_edge_norm_vectors(1:2, RIGHT_IDX, i, j)
        !   ! write(*, '(a, i3, 2(es16.6))') 'TOP_IDX   : ', TOP_IDX   , grid%cell_edge_norm_vectors(1:2, TOP_IDX, i, j)
        !   ! write(*, '(a, i3, 2(es16.6))') 'LEFT_IDX  : ', LEFT_IDX  , grid%cell_edge_norm_vectors(1:2, LEFT_IDX, i, j)

        !   ! call project_vector(a=[n_x, n_y], b=[u_LHS, v_LHS], b_parallel=LHS_par, b_perpendicular=LHS_per)
        !   ! error stop

        !   write(*,'(a, 2(i4), 10(es16.6))') "n_x, n_y        : ", i, j, n_x, n_y
        !   write(*,'(a, 2(i4), 10(es16.6))') "rho L/R 1/2     : ", i, j, rho_L_half, rho_R_half
        !   write(*,'(a, 2(i4), 10(es16.6))') "u   L/R         : ", i, j, u_LHS, u_RHS
        !   write(*,'(a, 2(i4), 10(es16.6))') "v   L/R         : ", i, j, v_LHS, v_RHS
        !   ! print*, "parallel     ", LHS_par, 'norm: ', norm2(LHS_par)
        !   ! print*, "perpendicular", LHS_per, 'norm: ', norm2(LHS_per)
        !   write(*,'(a, 2(i4), 10(es16.6))') "u   L/R 1/2     : ", i, j, u_L_half, u_R_half
        !   write(*,'(a, 2(i4), 10(es16.6))') "v   L/R 1/2     : ", i, j, v_L_half, v_R_half
        !   write(*,'(a, 2(i4), 10(es16.6))') "U_L, U_R        : ", i, j, U_L, U_R
        !   write(*,'(a, 2(i4), 10(es16.6))') "V_L, V_R        : ", i, j, V_L, V_R
        !   write(*,'(a, 2(i4), 10(es16.6))') "p   L/R 1/2     : ", i, j, p_L_half, p_R_half
        !   write(*,'(a, 2(i4), 10(es16.6))') "H   L/R 1/2     : ", i, j, H_L_half, H_R_half
        !   write(*,'(a, 2(i4), 10(es16.6))') "P_L+, P_R-      : ", i, j, P_L_plus, P_R_minus
        !   write(*,'(a, 2(i4), 10(es16.6))') "M_L, M_R        : ", i, j, M_L, M_R
        !   write(*,'(a, 2(i4), 10(es16.6))') "M_L+, M_R-      : ", i, j, M_L_plus, M_R_minus
        !   write(*,'(a, 2(i4), 10(es16.6))') "M_b_L_+, M_b_R_-: ", i, j, M_bar_L_plus, M_bar_R_minus
        !   write(*,'(a, 2(i4), 10(es16.6))') "f L/R, w, w1, w2: ", i, j, f_L, f_R, w, w1, w2(i,j)
        !   write(*,'(a, 2(i4), 10(es16.6))') "Edge flux       : ", i, j, edge_flux(:, i, j)
        !   print*
        !   ! error stop
        ! end if
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine get_iflux

  subroutine get_jflux(self, grid, lbounds, rho, u, v, p, rho_sb, u_sb, v_sb, p_sb, rho_ave, u_ave, v_ave, p_ave, w2,&
      & eflux_lbounds, edge_flux)
    !< Construct the fluxes for each edge in the j direction. This is templated via the Fypp pre-processor
    class(m_ausmpw_plus_solver_t), intent(inout) :: self
    class(grid_t), intent(in) :: grid
    integer(ik), dimension(3), intent(in) :: lbounds !< bounds of the primitive variable arrays
    real(rk), dimension(lbounds(2):, lbounds(3):), contiguous, intent(in) :: rho_ave    !< (i,j); cell averaged value ofdensity; needed for critical Mach number calcs
    real(rk), dimension(lbounds(2):, lbounds(3):), contiguous, intent(in) :: u_ave    !< (i,j); cell averaged value of x-velocity; needed for critical Mach number calcs
    real(rk), dimension(lbounds(2):, lbounds(3):), contiguous, intent(in) :: v_ave    !< (i,j); cell averaged value of y-velocity; needed for critical Mach number calcs
    real(rk), dimension(lbounds(2):, lbounds(3):), contiguous, intent(in) :: p_ave    !< (i,j); cell averaged value of pressure; needed for critical Mach number calcs
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: rho    !< (1:4, i,j); interpolated w/limiter of choice (L/R state) values for density
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: u      !< (1:4, i,j); interpolated w/limiter of choice (L/R state) values for x-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: v      !< (1:4, i,j); interpolated w/limiter of choice (L/R state) values for y-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: p      !< (1:4, i,j); interpolated w/limiter of choice (L/R state) values for pressure
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: rho_sb !< (1:4, i,j); interpolated w/superbee (L/R state) values for density
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: u_sb   !< (1:4, i,j); interpolated w/superbee (L/R state) values for x-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: v_sb   !< (1:4, i,j); interpolated w/superbee (L/R state) values for y-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: p_sb   !< (1:4, i,j); interpolated w/superbee (L/R state) values for pressure
    real(rk), dimension(lbounds(2):, lbounds(3):), contiguous, intent(in) :: w2 !< shock sensing function (determines if shock exists in transversal direction to the interface)
    integer(ik), dimension(3), intent(in) :: eflux_lbounds !< bounds of the primitive variable arrays
    real(rk), dimension(eflux_lbounds(1):, eflux_lbounds(2):, eflux_lbounds(3):), contiguous, intent(inout) :: edge_flux

    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi
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
    real(rk) :: m_half     !< interface Mach number
    real(rk) :: a               !< supersonic sensor
    real(rk) :: M_bar_L_plus    !< left final split Mach
    real(rk) :: M_bar_R_minus   !< right final split Mach
    real(rk) :: P_L_plus        !< left split pressure function
    real(rk) :: P_R_minus       !< right split pressure function
    real(rk) :: M_L_plus        !< left split Mach function
    real(rk) :: M_R_minus       !< right split Mach function
    real(rk) :: n_x             !< normal vectors of each face
    real(rk) :: n_y             !< normal vectors of each face
    real(rk) :: M_star, M_star_2, M_star_1
    real(rk) :: vel_ave, vel_ave_2
    real(rk) :: mass_flux_L, mass_flux_R

    integer(ik), parameter :: BOTTOM_IDX = 1 !< edge index for the bottom edge of the current cell
    integer(ik), parameter ::  RIGHT_IDX = 2 !< edge index for the right edge of the current cell
    integer(ik), parameter ::    TOP_IDX = 3 !< edge index for the top edge of the current cell
    integer(ik), parameter ::   LEFT_IDX = 4 !< edge index for the left edge of the current cell

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

    gamma = self%gamma

    ilo = lbound(edge_flux, dim=2)
    ihi = ubound(edge_flux, dim=2)
    jlo = lbound(edge_flux, dim=3)
    jhi = ubound(edge_flux, dim=3)

    !                   edge(3)
    !                  jflux(i,j+1)  'R'
    !               o--------------------o
    !               |                'L' |
    !               |                    |
    !   iflux(i, j) |     cell (i,j)     | iflux(i+1, j)
    !    edge(4)    |                    |  edge(2)
    !               |                'L' | 'R'
    !               o--------------------o
    !                  jflux(i,j)
    !                   edge(1)
    !  For jflux: "R" is the bottom of the cell above, "L" is the top of the current cell
    !  For iflux: "R" is the left of the cell to the right, "L" is the right of the current cell

    !$omp parallel default(none), &
    !$omp firstprivate(gamma, ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp private(rho_L, rho_R, u_LHS, u_RHS, v_LHS, v_RHS, p_L, p_R, m_half) &
    !$omp private(rho_L_sb, rho_R_sb, u_LHS_sb, u_RHS_sb, v_LHS_sb, v_RHS_sb, p_L_sb, p_R_sb) &
    !$omp private(n_x, n_y, p_s, U_L, U_R, V_L, V_R) &
    !$omp private(rho_L_half, rho_R_half ,u_L_half, u_R_half, v_L_half, v_R_half, p_L_half, p_R_half) &
    !$omp private(H_L_half, H_R_half, H_normal, M_L, M_R, H_L, H_R, f_L, f_R, w1, w,c_s ,c_half, a) &
    !$omp private(M_bar_L_plus, M_bar_R_minus, P_L_plus, P_R_minus, M_L_plus, M_R_minus) &
    !$omp private(mass_flux_L, mass_flux_R) &
    !$omp private(M_star, M_star_2, M_star_1, vel_ave, vel_ave_2) &
    !$omp shared(rho_ave, u_ave, v_ave, p_ave) &
    !$omp shared(grid, w2, rho, u, v, p, rho_sb, u_sb, v_sb, p_sb, edge_flux)
    !$omp do
    do j = jlo, jhi
      !$omp simd __INTERP_ALIGN__
      !dir$ vector aligned
      do i = ilo, ihi

        n_x = grid%cell_edge_norm_vectors(1, TOP_IDX, i, j)
        n_y = grid%cell_edge_norm_vectors(2, TOP_IDX, i, j)

        rho_L = rho(TOP_IDX, i, j)
        u_LHS = u(TOP_IDX, i, j)
        v_LHS = v(TOP_IDX, i, j)
        p_L = p(TOP_IDX, i, j)
        rho_R = rho(BOTTOM_IDX, i, j + 1)
        u_RHS = u(BOTTOM_IDX, i, j + 1)
        v_RHS = v(BOTTOM_IDX, i, j + 1)
        p_R = p(BOTTOM_IDX, i, j + 1)

        ! Superbee values
        rho_L_sb = rho_sb(TOP_IDX, i, j)
        u_LHS_sb = u_sb(TOP_IDX, i, j)
        v_LHS_sb = v_sb(TOP_IDX, i, j)
        p_L_sb = p_sb(TOP_IDX, i, j)
        rho_R_sb = rho_sb(BOTTOM_IDX, i, j + 1)
        u_RHS_sb = u_sb(BOTTOM_IDX, i, j + 1)
        v_RHS_sb = v_sb(BOTTOM_IDX, i, j + 1)
        p_R_sb = p_sb(BOTTOM_IDX, i, j + 1)

        U_L = v_LHS
        U_R = v_RHS
        V_L = u_LHS
        V_R = u_RHS

        ! Velocity normal to the edge, see Fig 2 in Ref[3]
        ! _RHS/_LHS is to avoid naming conflicts with _L and _R (slightly different meaning)
        ! U_L = u_LHS * n_x + v_LHS * n_y
        ! U_R = u_RHS * n_x + v_RHS * n_y

        ! ! Velocity component parallel to the edge
        ! V_L = u_LHS * (-n_y) + v_LHS * n_x
        ! V_R = u_RHS * (-n_y) + v_RHS * n_x

        H_L = (gamma / (gamma - 1.0_rk)) * (p_L / rho_L) + 0.5_rk * (u_LHS**2 + v_LHS**2)
        H_R = (gamma / (gamma - 1.0_rk)) * (p_R / rho_R) + 0.5_rk * (u_RHS**2 + v_RHS**2)

        ! Total enthalpy normal to the edge
        ! H_normal = min(H_L - 0.5_rk * V_L**2, H_R - 0.5_rk * V_R**2)
        H_normal = min(H_L - 0.5_rk * V_L**2, H_R - 0.5_rk * V_R**2)

        ! Speed of sound normal to the edge, also like the critical sound speed
        ! across a normal shock
        c_s = sqrt(2.0_rk * ((gamma - 1.0_rk) / (gamma + 1.0_rk)) * H_normal)

        ! Intermediate charactersistic Mach numbers M*
        M_star = sqrt(u_ave(i, j)**2 + v_ave(i, j)**2) / c_s ! M*(i)

        M_star_2 = sqrt(u_ave(i, j + 1)**2 + v_ave(i, j + 1)**2) / c_s ! M*(j+1)

        ! Interface sound speed
        if(0.5_rk * (U_L + U_R) < 0.0_rk) then  ! part (ii) in the paper after Eq 3
          c_half = c_s**2 / max(abs(U_R), c_s)
        else
          c_half = c_s**2 / max(abs(U_L), c_s) ! part (i)
        end if

        ! Left/Right Mach number
        M_L = U_L / c_half
        M_R = U_R / c_half
        a = 1.0_rk - min(1.0_rk, max(abs(M_L), abs(M_R)))**2

        ! Mach splitting functions
        M_L_plus = mach_split_plus(M_L)
        M_R_minus = mach_split_minus(M_R)

        vel_ave = sqrt(u_ave(i, j)**2 + v_ave(i, j)**2)
        vel_ave_2 = sqrt(u_ave(i, j + 1)**2 + v_ave(i, j + 1)**2)

        ! Pressure splitting functions
        ! Eq. 38a in Ref [1]
        if(M_star > 1.0_rk .and. M_star_2 < 1.0_rk .and. 0.0_rk < M_star * M_star_2 .and. M_star * M_star_2 < 1.0_rk) then
          P_R_minus = 1.0_rk - ((rho_ave(i, j) * vel_ave * (vel_ave - vel_ave_2) + p_ave(i, j)) / p_ave(i, j + 1))
        else
          P_R_minus = pressure_split_minus(M_R)
        end if

        ! Eq, 38b in Ref [1]
        if(M_star < -1.0_rk .and. M_star_2 < -1.0_rk .and. 0.0_rk < M_star * M_star_2 .and. M_star * M_star_2 < 1.0_rk) then
          P_L_plus = 1.0_rk - ((rho_ave(i, j + 1) * vel_ave_2 * (vel_ave_2 - vel_ave) + p_ave(i, j + 1)) / p_ave(i, j))
        else
          P_L_plus = pressure_split_plus(M_L)
        end if

        p_s = p_L * P_L_plus + p_R * P_R_minus
        w1 = w_1(p_L=p_L, p_R=p_R)
        w = max(w1, w2(i, j))

        f_L = f(p=p_L, p_s=p_s, w_2=w2(i, j))
        f_R = f(p=p_R, p_s=p_s, w_2=w2(i, j))

        ! Eq. 2b in Ref [1]
        if(M_L_plus + M_R_minus < 0.0_rk) then
          M_bar_L_plus = M_L_plus * w * (1.0_rk + f_L)
          M_bar_R_minus = M_R_minus + M_L_plus * ((1.0_rk - w) * (1.0_rk + f_L) - f_R)
        else ! Eq. 2a in Ref [1]
          M_bar_L_plus = M_L_plus + M_R_minus * ((1.0_rk - w) * (1.0_rk + f_R) - f_L)
          M_bar_R_minus = M_R_minus * w * (1.0_rk + f_R)
        end if

        rho_L_half = phi_L_half(phi_L=rho_L, phi_R=rho_R, phi_L_superbee=rho_L_sb, a=a)
        rho_R_half = phi_R_half(phi_L=rho_L, phi_R=rho_R, phi_R_superbee=rho_R_sb, a=a)

        u_L_half = phi_L_half(phi_L=u_LHS, phi_R=u_RHS, phi_L_superbee=u_LHS_sb, a=a)
        u_R_half = phi_R_half(phi_L=u_LHS, phi_R=u_RHS, phi_R_superbee=u_RHS_sb, a=a)

        v_L_half = phi_L_half(phi_L=v_LHS, phi_R=v_RHS, phi_L_superbee=v_LHS_sb, a=a)
        v_R_half = phi_R_half(phi_L=v_LHS, phi_R=v_RHS, phi_R_superbee=v_RHS_sb, a=a)

        p_L_half = phi_L_half(phi_L=p_L, phi_R=p_R, phi_L_superbee=p_L_sb, a=a)
        p_R_half = phi_R_half(phi_L=p_L, phi_R=p_R, phi_R_superbee=p_R_sb, a=a)

        ! Enthalpy
        H_L_half = (gamma / (gamma - 1.0_rk)) * (p_L_half / rho_L_half) + 0.5_rk * (u_L_half**2 + v_L_half**2)
        H_R_half = (gamma / (gamma - 1.0_rk)) * (p_R_half / rho_R_half) + 0.5_rk * (u_R_half**2 + v_R_half**2)

        mass_flux_L = M_bar_L_plus * c_half * rho_L_half
        mass_flux_R = M_bar_R_minus * c_half * rho_R_half
        edge_flux(1, i, j) = (mass_flux_L) + (mass_flux_R)
       edge_flux(2, i, j) = (mass_flux_L * u_L_half) + (mass_flux_R * u_R_half) + ((P_L_plus * n_x * p_L) + (P_R_minus * n_x * p_R))
       edge_flux(3, i, j) = (mass_flux_L * v_L_half) + (mass_flux_R * v_R_half) + ((P_L_plus * n_y * p_L) + (P_R_minus * n_y * p_R))
        edge_flux(4, i, j) = (mass_flux_L * H_L_half) + (mass_flux_R * H_R_half)

        ! if (i == 51 .or. i == 52) then

        !   ! write(*, '(a, i3, 2(es16.6))') 'BOTTOM_IDX: ', BOTTOM_IDX, grid%cell_edge_norm_vectors(1:2, BOTTOM_IDX, i, j)
        !   ! write(*, '(a, i3, 2(es16.6))') 'RIGHT_IDX : ', RIGHT_IDX , grid%cell_edge_norm_vectors(1:2, RIGHT_IDX, i, j)
        !   ! write(*, '(a, i3, 2(es16.6))') 'TOP_IDX   : ', TOP_IDX   , grid%cell_edge_norm_vectors(1:2, TOP_IDX, i, j)
        !   ! write(*, '(a, i3, 2(es16.6))') 'LEFT_IDX  : ', LEFT_IDX  , grid%cell_edge_norm_vectors(1:2, LEFT_IDX, i, j)

        !   ! call project_vector(a=[n_x, n_y], b=[u_LHS, v_LHS], b_parallel=LHS_par, b_perpendicular=LHS_per)
        !   ! error stop

        !   write(*,'(a, 2(i4), 10(es16.6))') "n_x, n_y        : ", i, j, n_x, n_y
        !   write(*,'(a, 2(i4), 10(es16.6))') "rho L/R 1/2     : ", i, j, rho_L_half, rho_R_half
        !   write(*,'(a, 2(i4), 10(es16.6))') "u   L/R         : ", i, j, u_LHS, u_RHS
        !   write(*,'(a, 2(i4), 10(es16.6))') "v   L/R         : ", i, j, v_LHS, v_RHS
        !   ! print*, "parallel     ", LHS_par, 'norm: ', norm2(LHS_par)
        !   ! print*, "perpendicular", LHS_per, 'norm: ', norm2(LHS_per)
        !   write(*,'(a, 2(i4), 10(es16.6))') "u   L/R 1/2     : ", i, j, u_L_half, u_R_half
        !   write(*,'(a, 2(i4), 10(es16.6))') "v   L/R 1/2     : ", i, j, v_L_half, v_R_half
        !   write(*,'(a, 2(i4), 10(es16.6))') "U_L, U_R        : ", i, j, U_L, U_R
        !   write(*,'(a, 2(i4), 10(es16.6))') "V_L, V_R        : ", i, j, V_L, V_R
        !   write(*,'(a, 2(i4), 10(es16.6))') "p   L/R 1/2     : ", i, j, p_L_half, p_R_half
        !   write(*,'(a, 2(i4), 10(es16.6))') "H   L/R 1/2     : ", i, j, H_L_half, H_R_half
        !   write(*,'(a, 2(i4), 10(es16.6))') "P_L+, P_R-      : ", i, j, P_L_plus, P_R_minus
        !   write(*,'(a, 2(i4), 10(es16.6))') "M_L, M_R        : ", i, j, M_L, M_R
        !   write(*,'(a, 2(i4), 10(es16.6))') "M_L+, M_R-      : ", i, j, M_L_plus, M_R_minus
        !   write(*,'(a, 2(i4), 10(es16.6))') "M_b_L_+, M_b_R_-: ", i, j, M_bar_L_plus, M_bar_R_minus
        !   write(*,'(a, 2(i4), 10(es16.6))') "f L/R, w, w1, w2: ", i, j, f_L, f_R, w, w1, w2(i,j)
        !   write(*,'(a, 2(i4), 10(es16.6))') "Edge flux       : ", i, j, edge_flux(:, i, j)
        !   print*
        !   ! error stop
        ! end if
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine get_jflux

  pure subroutine project_vector(a, b, b_parallel, b_perpendicular)
    !< Project vector "b" onto vector "a" and get the parallel and perpenticular components
    real(rk), dimension(2), intent(in) :: a                !<(x,y); the vector getting projected
    real(rk), dimension(2), intent(in) :: b                !<(x,y); the vector to project onto
    real(rk), dimension(2), intent(out) :: b_parallel      !<(x,y); vector parallel to a
    real(rk), dimension(2), intent(out) :: b_perpendicular !<(x,y); vector perpendicular to a

    b_parallel = (dot_product(b, a) / dot_product(a, a)) * a
    b_perpendicular = b - b_parallel
  end subroutine project_vector

  subroutine finalize(self)
    !< Cleanup the M-AUSMPW+ solver
    type(m_ausmpw_plus_solver_t), intent(inout) :: self
    call debug_print('Running m_ausmpw_plus_solver_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%iflux)) deallocate(self%iflux) ! these should already be deallocated
    if(allocated(self%jflux)) deallocate(self%jflux) ! these should already be deallocated
  end subroutine finalize

  pure real(rk) function pressure_split_plus(M) result(P_plus)
    !< The pressure splitting function (Eq. 4 in Ref[1]). This is the P+ version
    real(rk), intent(in) :: M !< interface Mach number

    if(abs(M) > 1.0_rk) then
      P_plus = 0.5_rk * (1.0_rk + sign(1.0_rk, M))
    else ! |M| <= 1
      P_plus = 0.25_rk * (M + 1.0_rk)**2 * (2.0_rk - M)
    endif
  end function pressure_split_plus

  pure real(rk) function pressure_split_minus(M) result(P_minus)
    !< The pressure splitting function (Eq. 4 in Ref[1]). This is the P- version
    real(rk), intent(in) :: M !< interface Mach number

    if(abs(M) > 1.0_rk) then
      P_minus = 0.5_rk * (1.0_rk - sign(1.0_rk, M))
    else ! |M| <= 1
      P_minus = 0.25_rk * (M - 1.0_rk)**2 * (2.0_rk + M)
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
    !< The Mach splitting function (Eq. 3 in Ref [1]). This is the M- version. This is kept
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

      if(abs(phi_superbee_minus_L) < epsilon(1.0_rk) .or. abs(phi_R_minus_L) < epsilon(1.0_rk)) then
        phi_L_half = phi_L
      else
        ! Eq. 27a in Ref [1]
        phi_L_half = phi_L + (max(0.0_rk, phi_R_minus_L * phi_superbee_minus_L) / &
                              (phi_R_minus_L * abs(phi_superbee_minus_L))) * &
                     min(a * 0.5_rk * abs(phi_R_minus_L), abs(phi_superbee_minus_L))
      end if
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

      if(abs(phi_superbee_minus_R) < epsilon(1.0_rk) .or. abs(phi_L_minus_R) < epsilon(1.0_rk)) then
        phi_R_half = phi_R
      else
        ! Eq. 27b in Ref [1]
        phi_R_half = phi_R + (max(0.0_rk, phi_L_minus_R * phi_superbee_minus_R) / &
                              (phi_L_minus_R * abs(phi_superbee_minus_R))) * &
                     min(a * 0.5_rk * abs(phi_L_minus_R), abs(phi_superbee_minus_R))
      end if
    end if
  end function phi_R_half

end module mod_mausmpw_plus_solver
