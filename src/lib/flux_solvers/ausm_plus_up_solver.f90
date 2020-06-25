module mod_ausm_plus_up_solver
  !< Summary: Provide a base ausm_plus_up solver class structure
  !< Date: 06/22/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !<     [1] M. S. Liou "A sequel to AUSM, Part II AUSM+-up for all speeds",
  !<         Journal of Computational Physics 214 (2006) 137–170, https://doi.org/10.1016/j.jcp.2005.09.020

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_edge_interpolator_factory, only: edge_interpolator_factory
  use mod_edge_interpolator, only: edge_iterpolator_t
  use mod_flux_solver, only: flux_solver_t
  use mod_eos, only: eos
  use mod_grid, only: grid_t
  use mod_input, only: input_t

  implicit none

  private
  public :: ausm_plus_up_solver_t

  type, extends(flux_solver_t) :: ausm_plus_up_solver_t
    !< Implementation of the AUSM+-up scheme
    private
    real(rk) :: mach_split_beta = 1.0_rk / 8.0_rk        !< beta parameter in Eq. 20 in Ref [1]
    ! real(rk) :: pressure_split_alpha = 3.0_rk / 16.0_rk  !< alpha parameter in Eq. 24 in Ref [1]
    real(rk) :: pressure_diffusion_coeff = 0.25_rk       !< K_p; pressure diffusion coefficient in Eq 21 in Ref [1]
    real(rk) :: sigma = 1.0_rk                           !< sigma; another pressure diffusion coefficient in Eq 21 in Ref [1]
    real(rk) :: pressure_flux_coeff = 0.75_rk            !< K_u; pressure flux coefficient in Eq 26 in Ref [1]
    real(rk) :: M_inf_sq = 1.0_rk                        !< Freestream Mach number squared

    real(rk), dimension(:, :, :), allocatable :: i_edge_flux !< ((1:4), i, j) edge flux of the i-direction edges
    real(rk), dimension(:, :, :), allocatable :: j_edge_flux !< ((1:4), i, j) edge flux of the j-direction edges
  contains
    ! Public methods
    procedure, public :: initialize => initialize_ausm_plus_up
    procedure, public :: solve => solve_ausm_plus_up
    procedure, public, pass(lhs) :: copy => copy_ausm_plus_up

    ! Private methods
    procedure, private :: flux_edges
    procedure, private :: apply_primitive_bc
    procedure, private, nopass :: split_mach_deg_1
    procedure, private, nopass :: split_mach_deg_2
    procedure, private :: split_mach_deg_4
    procedure, private :: split_pressure_deg_5
    procedure, private :: pressure_diffusion_term
    final :: finalize

    ! Operators
    generic :: assignment(=) => copy
  end type ausm_plus_up_solver_t

contains

  subroutine initialize_ausm_plus_up(self, grid, input)
    class(ausm_plus_up_solver_t), intent(inout) :: self
    class(grid_t), intent(in), target :: grid
    class(input_t), intent(in) :: input

    call debug_print('Running ausm_plus_up_solver_t%initialize_ausm_plus_up()', __FILE__, __LINE__)

    ! if (beta < (-1.0_rk / 16.0_rk) .or. beta > 0.5_rk) then
    !   error stop "Invalid value of beta in ausm_plus_up_solver_t%initialize_ausm_plus_up(); beta must be in the interval -1/16 <= beta <= 1/2"
    ! end if

    ! if (K_p < 0.0_rk or K_p > 1.0_rk) then
    !   error stop "Invalid value of K_p in ausm_plus_up_solver_t%pressure_diffusion_term();"//
    !             " K_p must be in the interval 0 <= K_p <= 1"
    !   end if

    ! if sigma <= 1 error stop
  end subroutine initialize_ausm_plus_up

  subroutine copy_ausm_plus_up(lhs, rhs)
    !< Implement LHS = RHS
    class(ausm_plus_up_solver_t), intent(inout) :: lhs
    type(ausm_plus_up_solver_t), intent(in) :: rhs

    call debug_print('Running ausm_plus_up_solver_t%copy()', __FILE__, __LINE__)

    ! allocate(lhs%bc_plus_x, source=rhs%bc_plus_x)
    ! allocate(lhs%bc_plus_y, source=rhs%bc_plus_y)
    ! allocate(lhs%bc_minus_x, source=rhs%bc_minus_x)
    ! allocate(lhs%bc_minus_y, source=rhs%bc_minus_y)

    lhs%iteration = rhs%iteration
    lhs%time = rhs%time
    lhs%dt = rhs%dt

  end subroutine copy_ausm_plus_up

  subroutine solve_ausm_plus_up(self, dt, grid, lbounds, rho, u, v, p, d_rho_dt, d_rho_u_dt, d_rho_v_dt, d_rho_E_dt)
    !< Solve and flux the edges
    class(ausm_plus_up_solver_t), intent(inout) :: self
    class(grid_t), intent(in) :: grid
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), intent(in) :: dt !< timestep (not really used in this solver, but needed for others)
    real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: rho !< density
    real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: u   !< x-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: v   !< y-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: p   !< pressure

    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) ::   d_rho_dt    !< d/dt of the density field
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: d_rho_u_dt    !< d/dt of the rhou field
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: d_rho_v_dt    !< d/dt of the rhov field
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: d_rho_E_dt    !< d/dt of the rhoE field

    ! Locals
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk), dimension(:, :, :), allocatable :: rho_interface_values !< interface (L/R state) values for density
    real(rk), dimension(:, :, :), allocatable :: u_interface_values   !< interface (L/R state) values for x-velocity
    real(rk), dimension(:, :, :), allocatable :: v_interface_values   !< interface (L/R state) values for y-velocity
    real(rk), dimension(:, :, :), allocatable :: p_interface_values   !< interface (L/R state) values for pressure

    class(edge_iterpolator_t), pointer :: edge_interpolator => null()

    real(rk) :: rho_L !< density of the left side
    real(rk) :: rho_R !< density of the right side
    real(rk) :: u_L   !< x-velocity of the left side
    real(rk) :: u_R   !< x-velocity of the right side
    real(rk) :: v_L   !< y-velocity of the left side
    real(rk) :: v_R   !< y-velocity of the right side
    real(rk) :: vel_L !< total-velocity of the left side
    real(rk) :: vel_R !< total-velocity of the right side
    real(rk) :: p_L   !< pressure of the left side
    real(rk) :: p_R   !< pressure of the right side

    real(rk) :: H_L            !< total enthalpy of the left side
    real(rk) :: H_R            !< total enthalpy of the right side
    real(rk) :: a_circumflex_L !< â_L; used to find the left interface sound speed see Eq 28 in Ref [1]
    real(rk) :: a_circumflex_R !< â_R; used to find the rgight interface sound speed see Eq 28 in Ref [1]
    real(rk) :: a_crit_L       !< a*; critical sound speed of the left side
    real(rk) :: a_crit_R       !< a*; critical sound speed of the right state

    real(rk) :: M_p      !< Pressure diffusion term
    real(rk) :: M_plus   !< Split Mach term M+, see Eq. 73
    real(rk) :: M_minus  !< Split Mach term M-, see Eq. 73
    real(rk) :: M_bar_sq !< Mean Mach number ^2
    real(rk) :: M_0      !< Reference Mach
    real(rk) :: M_0_sq   !< Reference Mach squared
    real(rk) :: M_L      !< Mach number of the left side of the interface
    real(rk) :: M_R      !< Mach number of the right side of the interface
    real(rk) :: M_half   !< Final interface Mach number
    real(rk) :: rho_half !< Final interface density
    real(rk) :: a_half   !< Final interface sound speed
    real(rk) :: f_a      !< Scaling factor: f_a(M_0) = M0(2-M0) -> [0,1]
    real(rk) :: alpha    !< Pressure split alpha factor, see Eq. 76

    real(rk) :: p_u      !< Velocity difference (diffusion) term
    real(rk) :: p_half   !< Final interface pressure
    real(rk) :: P_plus   !< Split pressure term P+, see Eq. 75
    real(rk) :: P_minus  !< Split pressure term P-, see Eq. 75

    integer(ik), parameter :: bottom = 1 !< edge index for the bottom edge of the current cell
    integer(ik), parameter :: right = 2 !< edge index for the right edge of the current cell
    integer(ik), parameter :: top = 3 !< edge index for the top edge of the current cell
    integer(ik), parameter :: left = 4 !< edge index for the left edge of the current cell

    real(rk) :: n_x  !< normal vectors of each face
    real(rk) :: n_y  !< normal vectors of each face
    real(rk) :: gamma        !< EOS polytropic index
    real(rk) :: gamma_factor !< Precomputed gamma factor = 2*(gamma-1)/(gamma+1)

    ! Split the convective flux (F_c) into the convective part and the pressure part
    !         [rho  ] + [0    ]
    !         [rho u] + [n_x p]
    ! F_c = V [rho v] + [n_y p]
    !         [rho H] + [0    ]
    !
    ! Which becomes the following: (where a is sound speed, M is Mach number)
    !                         [rho a  ]     +    [0    ]
    !                         [rho a u]     +    [n_x p]
    ! F_c_(i+1/2) = M_(i+1/2) [rho a v]     +    [n_y p]
    !                         [rho a H]     +    [0    ]
    !                                 (L/R)            (i+1/2)
    !
    ! with (L/R) = [L  if M_(i+1/2) >= 0]
    !              [R  otherwise        ]

    call debug_print('Running ausm_plus_up_solver_t%solve_ausm_plus_up()', __FILE__, __LINE__)

    if(dt < tiny(1.0_rk)) error stop "Error in ausm_plus_up_solver_t%solve_ausm_plus_up(), the timestep dt is < tiny(1.0_rk)"
    self%time = self%time + dt
    self%dt = dt
    self%iteration = self%iteration + 1

    ! Precompute a few re-used scalars
    gamma = eos%get_gamma()
    gamma_factor = ((2.0_rk * (gamma - 1.0_rk)) / (gamma + 1.0_rk))

    ! call self%apply_primitive_bc() !FIXME:

    ! Get the edge values (left/right states for each edge)
    edge_interpolator => edge_interpolator_factory(self%input)
    call edge_interpolator%interpolate_edge_values(q=rho, lbounds=lbounds, edge_values=rho_interface_values)
    call edge_interpolator%interpolate_edge_values(q=u, lbounds=lbounds, edge_values=u_interface_values)
    call edge_interpolator%interpolate_edge_values(q=v, lbounds=lbounds, edge_values=v_interface_values)
    call edge_interpolator%interpolate_edge_values(q=p, lbounds=lbounds, edge_values=p_interface_values)

    ! The cell interfaces are shared by all the cells, so the indexing here is key.
    ! We want to avoid any i-1 or j-1 indexing for vectorization dependence, so we then
    ! make sure to use i or i+1
    !
    !             (i-1/2)  (i+1/2)         interface locations
    !        |       |       |       |
    !        |       |     L | R     |
    !        |       |       |       |
    !        | (i-1) |  (i)  | (i+1) |     cell indexing
    !        |       |       |       |
    !        |       |       |       |
    !        |       |       |       |
    !              (i-1)    (i)    (i+1)   edge indexing
    !           0        1       2         e.g. on left boundary
    !         ilo_bc    ilo      i         i indices

    ilo = grid%ilo_cell
    ihi = grid%ihi_cell
    jlo = grid%jlo_cell
    jhi = grid%jhi_cell

    if(.not. allocated(self%i_edge_flux)) allocate(self%i_edge_flux(4, ilo:ihi + 1, jlo:jhi))
    if(.not. allocated(self%j_edge_flux)) allocate(self%j_edge_flux(4, ilo:ihi, jlo:jhi + 1))

    ! FIXME: make sure to get the indexing right!!!
    ! Find fluxes in the i-direction
    do j = jlo, jhi
      do i = ilo - 1, ihi ! start at ilo-1 b/c of the first i edge

        ! use the normal vector of the right edge of the current cell
        n_x = grid%cell_edge_norm_vectors(1, right, i, j)
        n_y = grid%cell_edge_norm_vectors(2, right, i, j)

        ! right edge (i + 1/2)
        rho_L = rho_interface_values(right, i, j)    ! right edge of the current cell
        rho_R = rho_interface_values(left, i + 1, j) ! left edge of the cell to the right
        u_L = rho_interface_values(right, i, j)      ! right edge of the current cell
        u_R = rho_interface_values(left, i + 1, j)   ! left edge of the cell to the right
        v_L = rho_interface_values(right, i, j)      ! right edge of the current cell
        v_R = rho_interface_values(left, i + 1, j)   ! left edge of the cell to the right
        p_L = rho_interface_values(right, i, j)      ! right edge of the current cell
        p_R = rho_interface_values(left, i + 1, j)   ! left edge of the cell to the right

        ! Dot-product to get total velocity based on the edge normals
        vel_L = u_L * n_x + v_L * n_y
        vel_R = u_R * n_x + v_R * n_y

        ! Find the total enthalpy in order to get the critical sound speed of the interface
        call eos%total_enthalpy(rho=rho_L, u=u_L, v=v_L, p=p_L, H=H_L)
        call eos%total_enthalpy(rho=rho_R, u=u_R, v=v_R, p=p_R, H=H_R)

        ! Critical sound speed, Eq. 29
        a_crit_L = gamma_factor * H_L
        a_crit_R = gamma_factor * H_R

        ! Eq. 28
        a_circumflex_L = a_crit_L**2 / max(a_crit_L, abs(vel_L))
        a_circumflex_R = a_crit_R**2 / max(a_crit_R, abs(vel_R))

        ! Interface sound speed, Eq. 28
        a_half = min(a_circumflex_L, a_circumflex_R)

        ! Left/Right Mach, Eq. 69
        M_L = vel_L / a_half
        M_R = vel_R / a_half

        ! Mean Mach Eq. 70
        M_bar_sq = (vel_L**2 + vel_R**2) / (2.0_rk * a_half**2)

        ! Reference Mach number, see Eq. 71
        M_0_sq = min(1.0_rk, max(M_bar_sq, self%M_inf_sq))
        M_0 = sqrt(M_0_sq)

        ! Scaling function, Eq. 72
        f_a = M_0 * (2.0_rk - M_0)

        ! alpha parameter, Eq. 76
        alpha = (3.0_rk / 16.0_rk) * (-4.0_rk + 5.0_rk * f_a**2)

        ! Interface Mach number, Eq. 73
        M_plus = self%split_mach_deg_4(M=M_L, plus_or_minus='+')
        M_minus = self%split_mach_deg_4(M=M_R, plus_or_minus='-')

        ! Interface density
        rho_half = 0.5_rk * (rho_R + rho_L)

        !< The pressure diffusion term, Mp, as defined in Eq. 21 in Ref [1]
        associate(K_p=>self%pressure_diffusion_coeff, sigma=>self%sigma)
          M_p = -(K_p / f_a) * max(1.0_rk - sigma * M_bar_sq, 0.0_rk) * ((p_R - p_L) / (rho_half * a_half**2))
        end associate

        ! Interface Mach number
        M_half = M_plus + M_minus + M_p

        ! Pressure splitting functions in Eq. 75
        P_plus = self%split_pressure_deg_5(M=M_L, plus_or_minus='+', alpha=alpha)
        P_minus = self%split_pressure_deg_5(M=M_R, plus_or_minus='-', alpha=alpha)

        ! Velocity difference (diffusion) term p_u, NOT, Eq. 26, but rather the last term in Eq 75
        associate(K_u=>self%pressure_flux_coeff)
          p_u = -K_u * P_plus * P_minus * (rho_L + rho_R) * (f_a * a_half) * (u_R - u_L)
        end associate

        ! Interface pressure, Eq. 75
        p_half = P_plus * p_L + P_minus * p_R + p_u

        ! Edge flux values, via upwinding
        if(M_half > 0.0_rk) then
          self%i_edge_flux(1, i, j) = M_half * a_half * rho_L
          self%i_edge_flux(2, i, j) = M_half * a_half * rho_L * u_L + n_x * p_half
          self%i_edge_flux(3, i, j) = M_half * a_half * rho_L * v_L + n_y * p_half
          self%i_edge_flux(4, i, j) = M_half * a_half * rho_L * H_L
        else
          self%i_edge_flux(1, i, j) = M_half * a_half * rho_R
          self%i_edge_flux(2, i, j) = M_half * a_half * rho_R * u_R + n_x * p_half
          self%i_edge_flux(3, i, j) = M_half * a_half * rho_R * v_R + n_y * p_half
          self%i_edge_flux(4, i, j) = M_half * a_half * rho_R * H_R
        end if
      end do
    end do

    ! Find fluxes in the i-direction
    do j = jlo - 1, jhi
      do i = ilo, ihi
        n_x = grid%cell_edge_norm_vectors(1, right, i, j) ! FIXME:
        n_y = grid%cell_edge_norm_vectors(2, right, i, j) ! FIXME:

        ! right edge (i + 1/2)
        !FIXME:
        rho_L = rho_interface_values(right, i, j)     ! right edge of the current cell
        rho_R = rho_interface_values(left, i + 1, j)  ! left edge of the cell to the right

        u_L = rho_interface_values(right, i, j)     ! right edge of the current cell
        u_R = rho_interface_values(left, i + 1, j)  ! left edge of the cell to the right

        v_L = rho_interface_values(right, i, j)     ! right edge of the current cell
        v_R = rho_interface_values(left, i + 1, j)  ! left edge of the cell to the right

        ! Dot-product to get total velocity based on the edge normals
        vel_L = u_L * n_x + v_L * n_y
        vel_R = u_R * n_x + v_R * n_y

        p_L = rho_interface_values(right, i, j)     ! right edge of the current cell
        p_R = rho_interface_values(left, i + 1, j)  ! left edge of the cell to the right

        ! Find the total enthalpy in order to get the critical sound speed of the interface
        call eos%total_enthalpy(rho=rho_L, u=u_L, v=v_L, p=p_L, H=H_L)
        call eos%total_enthalpy(rho=rho_R, u=u_R, v=v_R, p=p_R, H=H_R)

        ! Critical sound speed, Eq. 29
        a_crit_L = gamma_factor * H_L
        a_crit_R = gamma_factor * H_R

        ! Eq. 28
        a_circumflex_L = a_crit_L**2 / max(a_crit_L, abs(vel_L))
        a_circumflex_R = a_crit_R**2 / max(a_crit_R, abs(vel_R))

        ! Interface sound speed, Eq. 28
        a_half = min(a_circumflex_L, a_circumflex_R)

        ! Left/Right Mach, Eq. 69
        M_L = vel_L / a_half
        M_R = vel_R / a_half

        ! Mean Mach Eq. 70
        M_bar_sq = (vel_L**2 + vel_R**2) / (2.0_rk * a_half**2)

        ! Reference Mach number, see Eq. 71
        M_0_sq = min(1.0_rk, max(M_bar_sq, self%M_inf_sq))
        M_0 = sqrt(M_0_sq)

        ! Scaling function, Eq. 72
        f_a = M_0 * (2.0_rk - M_0)

        ! alpha parameter, Eq. 76
        alpha = (3.0_rk / 16.0_rk) * (-4.0_rk + 5.0_rk * f_a**2)

        ! Interface Mach number, Eq. 73
        M_plus = self%split_mach_deg_4(M=M_L, plus_or_minus='+')
        M_minus = self%split_mach_deg_4(M=M_R, plus_or_minus='-')

        ! Interface density
        rho_half = 0.5_rk * (rho_R + rho_L)

        !< The pressure diffusion term, Mp, as defined in Eq. 21 in Ref [1]
        associate(K_p=>self%pressure_diffusion_coeff, sigma=>self%sigma)
          M_p = -(K_p / f_a) * max(1.0_rk - sigma * M_bar_sq, 0.0_rk) * ((p_R - p_L) / (rho_half * a_half**2))
        end associate

        ! Interface Mach number
        M_half = M_plus + M_minus + M_p

        ! Pressure splitting functions in Eq. 75
        P_plus = self%split_pressure_deg_5(M=M_L, plus_or_minus='+', alpha=alpha)
        P_minus = self%split_pressure_deg_5(M=M_R, plus_or_minus='-', alpha=alpha)

        ! Velocity difference (diffusion) term p_u, NOT, Eq. 26, but rather the last term in Eq 75
        associate(K_u=>self%pressure_flux_coeff)
          p_u = -K_u * P_plus * P_minus * (rho_L + rho_R) * (f_a * a_half) * (u_R - u_L)
        end associate

        ! Interface pressure, Eq. 75
        p_half = P_plus * p_L + P_minus * p_R + p_u

        ! Edge flux values, via upwinding
        if(M_half > 0.0_rk) then
          self%j_edge_flux(1, i, j) = M_half * a_half * rho_L
          self%j_edge_flux(2, i, j) = M_half * a_half * rho_L * u_L + n_x * p_half
          self%j_edge_flux(3, i, j) = M_half * a_half * rho_L * v_L + n_y * p_half
          self%j_edge_flux(4, i, j) = M_half * a_half * rho_L * H_L
        else
          self%j_edge_flux(1, i, j) = M_half * a_half * rho_R
          self%j_edge_flux(2, i, j) = M_half * a_half * rho_R * u_R + n_x * p_half
          self%j_edge_flux(3, i, j) = M_half * a_half * rho_R * v_R + n_y * p_half
          self%j_edge_flux(4, i, j) = M_half * a_half * rho_R * H_R
        end if
      end do
    end do

    call self%flux_edges(grid, lbounds, d_rho_dt, d_rho_u_dt, d_rho_v_dt, d_rho_E_dt)
    deallocate(rho_interface_values)
    deallocate(u_interface_values)
    deallocate(v_interface_values)
    deallocate(p_interface_values)
    deallocate(edge_interpolator)
    if(allocated(self%i_edge_flux)) deallocate(self%i_edge_flux)
    if(allocated(self%j_edge_flux)) deallocate(self%j_edge_flux)
  end subroutine solve_ausm_plus_up

  subroutine flux_edges(self, grid, lbounds, d_rho_dt, d_rho_u_dt, d_rho_v_dt, d_rho_E_dt)
    !< Flux the edges to get the residuals, e.g. 1/vol * d/dt U
    class(ausm_plus_up_solver_t), intent(in) :: self
    class(grid_t), intent(in) :: grid
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) ::   d_rho_dt    !< d/dt of the density field
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: d_rho_u_dt    !< d/dt of the rhou field
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: d_rho_v_dt    !< d/dt of the rhov field
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: d_rho_E_dt    !< d/dt of the rhoE field

  end subroutine flux_edges

  pure function split_mach_deg_1(M, plus_or_minus) result(M_split)
    !< Implementation of the M +- of degree 1. See Eq. 18 in Ref [1]
    real(rk), intent(in) :: M         !< Mach number
    character(len=1), intent(in) :: plus_or_minus !< which split? '+' or '-'
    real(rk) :: M_split !< Split Mach number

    if(plus_or_minus == '+') then
      M_split = 0.5_rk * (M + abs(M))
    else
      M_split = 0.5_rk * (M - abs(M))
    end if
  end function split_mach_deg_1

  pure function split_mach_deg_2(M, plus_or_minus) result(M_split)
    !< Implementation of the M +- of degree 2. See Eq. 19 in Ref [1]
    real(rk), intent(in) :: M         !< Mach number
    character(len=1), intent(in) :: plus_or_minus !< which split? '+' or '-'
    real(rk) :: M_split !< Split Mach number

    if(plus_or_minus == '+') then
      M_split = 0.25_rk * (M + 1.0_rk)**2
    else
      M_split = 0.25_rk * (M - 1.0_rk)**2
    end if

  end function split_mach_deg_2

  pure function split_mach_deg_4(self, M, plus_or_minus) result(M_split)
    !< Implementation of the M +- of degree 4. See Eq. 20 in Ref [1]

    class(ausm_plus_up_solver_t), intent(in) :: self
    real(rk), intent(in) :: M         !< Mach number
    character(len=1), intent(in) :: plus_or_minus !< which split? '+' or '-'
    real(rk) :: M_split !< Split Mach number

    ! Locals
    real(rk) :: M_split_deg_2 !< Split Mach number using a 2nd deg polynomial

    if(plus_or_minus == '+') then
      if(abs(M) >= 1.0_rk) then
        M_split = self%split_mach_deg_1(M, '+')
      else
        M_split_deg_2 = self%split_mach_deg_2(M, '+')
        associate(beta=>self%mach_split_beta)
          M_split = M_split_deg_2 * (1.0_rk + 16.0_rk * beta * M_split_deg_2)
        end associate
      end if
    else ! '-'
      if(abs(M) >= 1.0_rk) then
        M_split = self%split_mach_deg_1(M, '-')
      else
        M_split_deg_2 = self%split_mach_deg_2(M, '-')
        associate(beta=>self%mach_split_beta)
          M_split = M_split_deg_2 * (1.0_rk - 16.0_rk * beta * M_split_deg_2)
        end associate
      end if
    end if
  end function split_mach_deg_4

  pure function split_pressure_deg_5(self, M, plus_or_minus, alpha) result(P_split)
    !< Implementation of the pressure P +- split factor of degree 5. See Eq. 24 in Ref [1]

    class(ausm_plus_up_solver_t), intent(in) :: self
    real(rk), intent(in) :: M                     !< Mach number
    real(rk), intent(in) :: alpha
    character(len=1), intent(in) :: plus_or_minus !< Which split? '+' or '-'
    real(rk) :: P_split                           !< Split pressure factor

    ! Locals
    real(rk) :: M_split_deg_2 !< Split Mach number using a 2nd deg polynomial

    if(plus_or_minus == '+') then
      if(abs(M) >= 1.0_rk) then
        P_split = self%split_mach_deg_1(M, '+') / M
      else
        M_split_deg_2 = self%split_mach_deg_2(M, '+')
        P_split = M_split_deg_2 * ((2.0_rk - M) + 16.0_rk * alpha * M * M_split_deg_2)
      end if
    else ! '-'
      if(abs(M) >= 1.0_rk) then
        P_split = self%split_mach_deg_1(M, '-') / M
      else
        M_split_deg_2 = self%split_mach_deg_2(M, '-')
        P_split = M_split_deg_2 * ((-2.0_rk - M) - 16.0_rk * alpha * M * M_split_deg_2)
      end if
    end if
  end function split_pressure_deg_5

  pure real(rk) function pressure_diffusion_term(self, rho_R, rho_L, p_R, p_L, a_half, M_bar_sq) result(M_p)
    !< The pressure diffusion term, Mp, as defined in Eq. 21 in Ref [1]

    class(ausm_plus_up_solver_t), intent(in) :: self
    real(rk), intent(in) :: rho_R    !< density of the RHS of the interface
    real(rk), intent(in) :: rho_L    !< density of the LHS of the interface
    real(rk), intent(in) :: p_R      !< pressure of the RHS of the interface
    real(rk), intent(in) :: p_L      !< pressure of the LHS of the interface
    real(rk), intent(in) :: a_half   !< sound speed of the interface
    real(rk), intent(in) :: M_bar_sq !< mean local Mach number, e.g. \bar(M)^2

    real(rk) :: rho_half !< mean density

    rho_half = 0.5_rk * (rho_R + rho_L)

    associate(K_p=>self%pressure_diffusion_coeff, sigma=>self%sigma)
      M_p = -K_p * max(1.0_rk - sigma * M_bar_sq, 0.0_rk) * ((p_R - p_L) / (rho_half * a_half**2))
    end associate
  end function pressure_diffusion_term

  subroutine apply_primitive_bc(self, lbounds, rho, u, v, p, &
                                bc_plus_x, bc_minus_x, bc_plus_y, bc_minus_y)
    class(ausm_plus_up_solver_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: rho
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: u
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: v
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: p
    class(boundary_condition_t), intent(inout):: bc_plus_x
    class(boundary_condition_t), intent(inout):: bc_plus_y
    class(boundary_condition_t), intent(inout):: bc_minus_x
    class(boundary_condition_t), intent(inout):: bc_minus_y

    integer(ik) :: priority
    integer(ik) :: max_priority_bc !< highest goes first

    call debug_print('Running fvleg_solver_t%apply_primitive_var_bc()', __FILE__, __LINE__)

    max_priority_bc = max(bc_plus_x%priority, bc_plus_y%priority, &
                          bc_minus_x%priority, bc_minus_y%priority)

    do priority = max_priority_bc, 0, -1

      if(bc_plus_x%priority == priority) then
        call bc_plus_x%apply_primitive_var_bc(rho=rho, u=u, v=v, p=p, lbounds=lbounds)
      end if

      if(bc_plus_y%priority == priority) then
        call bc_plus_y%apply_primitive_var_bc(rho=rho, u=u, v=v, p=p, lbounds=lbounds)
      end if

      if(bc_minus_x%priority == priority) then
        call bc_minus_x%apply_primitive_var_bc(rho=rho, u=u, v=v, p=p, lbounds=lbounds)
      end if

      if(bc_minus_y%priority == priority) then
        call bc_minus_y%apply_primitive_var_bc(rho=rho, u=u, v=v, p=p, lbounds=lbounds)
      end if

    end do

  end subroutine apply_primitive_bc

  subroutine finalize(self)
    !< Class finalizer
    type(ausm_plus_up_solver_t), intent(inout) :: self

    call debug_print('Running ausm_plus_up_solver_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%i_edge_flux)) deallocate(self%i_edge_flux) ! these should already be deallocated
    if(allocated(self%j_edge_flux)) deallocate(self%j_edge_flux) ! these should already be deallocated

    ! if(allocated(self%reconstructor)) deallocate(self%reconstructor)
    ! if(allocated(self%bc_plus_x)) deallocate(self%bc_plus_x)
    ! if(allocated(self%bc_plus_y)) deallocate(self%bc_plus_y)
    ! if(allocated(self%bc_minus_x)) deallocate(self%bc_minus_x)
    ! if(allocated(self%bc_minus_y)) deallocate(self%bc_minus_y)
  end subroutine finalize
end module mod_ausm_plus_up_solver
