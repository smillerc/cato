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

module mod_ausm_plus_solver
  !< Summary: Provide a solver based on the AUSM+ family of schemes
  !< Date: 06/22/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !<     [1] M. S. Liou "A sequel to AUSM, Part II AUSM+-up for all speeds",
  !<         Journal of Computational Physics 214 (2006) 137–170, https://doi.org/10.1016/j.jcp.2005.09.020

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
  public :: ausm_plus_solver_t

  type, extends(edge_split_flux_solver_t) :: ausm_plus_solver_t
    !< Implementation of the AUSM+-up scheme
    private
    real(rk) :: mach_split_beta = 1.0_rk / 8.0_rk        !< beta parameter in Eq. 20 in Ref [1]
    real(rk) :: pressure_split_alpha = 3.0_rk / 16.0_rk  !< beta parameter in Eq. 20 in Ref [1]
    real(rk) :: pressure_diffusion_coeff = 0.25_rk       !< K_p; pressure diffusion coefficient in Eq 21 in Ref [1]
    real(rk) :: pressure_flux_coeff = 0.75_rk            !< K_u; pressure flux coefficient in Eq 26 in Ref [1]
    real(rk) :: sigma = 1.0_rk                           !< sigma; another pressure diffusion coefficient in Eq 21 in Ref [1]
    real(rk) :: M_inf_sq = 1.0_rk                        !< Freestream Mach number squared

    ! Variants of the AUSM+ family
    logical :: ausm_plus_u = .false.            !< Enable the AUSM+-u scheme
    logical :: ausm_plus_up_basic = .false.     !< Enable the AUSM+-up basic scheme
    logical :: ausm_plus_up_all_speed = .false. !< Enable the AUSM+-up all-speed scheme

  contains
    ! Public methods
    procedure, public :: initialize => initialize_ausm_plus
    procedure, public :: solve => solve_ausm_plus
    procedure, public, pass(lhs) :: copy => copy_ausm_plus

    ! Private methods
    procedure, private :: interface_state
    final :: finalize

    ! Operators
    generic :: assignment(=) => copy
  end type ausm_plus_solver_t

contains

  subroutine initialize_ausm_plus(self, input)
    class(ausm_plus_solver_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    call debug_print('Running ausm_plus_solver_t%initialize_ausm_plus()', __FILE__, __LINE__)

    self%input = input
    self%M_inf_sq = input%reference_mach**2

    select case(trim(input%flux_solver))
    case('AUSM+-u')
      self%name = 'AUSM+-u_'//input%limiter
      self%ausm_plus_u = .true.
    case('AUSM+-up', 'AUSM+-up_basic')
      self%name = 'AUSM+-up'//input%limiter
      self%ausm_plus_up_basic = .true.
    case('AUSM+-up_all_speed')
      self%name = 'AUSM+-up all-speed'//input%limiter
      self%ausm_plus_up_all_speed = .true.
    case default
      write(std_err, '(a)') "Unknown variant of the AUSM+ family, must be one of the following"// &
        "['AUSM+-u','AUSM+-up', 'AUSM+-up_all_speed'], input was '"//trim(input%flux_solver)//"'"

      error stop "Unknown variant of the AUSM+ family, must be one of the following"// &
        "['AUSM+-u', 'AUSM+-up','AUSM+-up_all_speed']"
    end select

    if(input%ausm_beta < (-1.0_rk / 16.0_rk) .or. input%ausm_beta > 0.5_rk) then
      error stop "Invalid value of input%ausm_beta in ausm_plus_solver_t%initialize_ausm_plus(); "// &
        "beta must be in the interval -1/16 <= beta <= 1/2"
    else
      self%mach_split_beta = input%ausm_beta
    end if

    if(input%ausm_pressure_diffusion_coeff < 0.0_rk .or. input%ausm_pressure_diffusion_coeff > 1.0_rk) then
      error stop "Invalid value of the input%ausm_pressure_diffusion_coeff (K_p) "// &
        "in ausm_plus_solver_t%initialize_ausm_plus(); "// &
        "K_p must be in the interval 0 <= K_p <= 1"
    else
      self%pressure_diffusion_coeff = input%ausm_pressure_diffusion_coeff
    end if

    if(input%ausm_sonic_point_sigma < 0.25_rk .or. input%ausm_sonic_point_sigma > 1.0_rk) then
      error stop "Invalid value of the sonic point resolution parameter input%ausm_sonic_point_sigma in "// &
        "ausm_plus_solver_t%initialize_ausm_plus(); sigma must be in "// &
        "the interval  0.25 <= sigma < 1"
    else
      self%sigma = input%ausm_sonic_point_sigma
    end if
  end subroutine initialize_ausm_plus

  subroutine copy_ausm_plus(lhs, rhs)
    !< Implement LHS = RHS
    class(ausm_plus_solver_t), intent(inout) :: lhs
    type(ausm_plus_solver_t), intent(in) :: rhs

    call debug_print('Running ausm_plus_solver_t%copy()', __FILE__, __LINE__)

    ! allocate(lhs%bc_plus_x, source=rhs%bc_plus_x)
    ! allocate(lhs%bc_plus_y, source=rhs%bc_plus_y)
    ! allocate(lhs%bc_minus_x, source=rhs%bc_minus_x)
    ! allocate(lhs%bc_minus_y, source=rhs%bc_minus_y)

    lhs%iteration = rhs%iteration
    lhs%time = rhs%time
    lhs%dt = rhs%dt

  end subroutine copy_ausm_plus

  subroutine solve_ausm_plus(self, dt, grid, lbounds, rho, u, v, p, d_rho_dt, d_rho_u_dt, d_rho_v_dt, d_rho_E_dt)
    !< Solve and flux the edges
    class(ausm_plus_solver_t), intent(inout) :: self
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

    ! Locals
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk), dimension(:, :, :), allocatable :: rho_interface_values !< interface (L/R state) values for density
    real(rk), dimension(:, :, :), allocatable :: u_interface_values   !< interface (L/R state) values for x-velocity
    real(rk), dimension(:, :, :), allocatable :: v_interface_values   !< interface (L/R state) values for y-velocity
    real(rk), dimension(:, :, :), allocatable :: p_interface_values   !< interface (L/R state) values for pressure

    class(edge_iterpolator_t), pointer :: edge_interpolator => null()
    class(boundary_condition_t), allocatable:: bc_plus_x
    class(boundary_condition_t), allocatable:: bc_plus_y
    class(boundary_condition_t), allocatable:: bc_minus_x
    class(boundary_condition_t), allocatable:: bc_minus_y

    real(rk) :: rho_L !< density of the left side
    real(rk) :: rho_R !< density of the right side
    real(rk) :: u_L   !< x-velocity of the left side
    real(rk) :: u_R   !< x-velocity of the right side
    real(rk) :: v_L   !< y-velocity of the left side
    real(rk) :: v_R   !< y-velocity of the right side
    real(rk) :: p_L   !< pressure of the left side
    real(rk) :: p_R   !< pressure of the right side
    real(rk), dimension(2) :: H !< (L,R); Total enthalpy

    real(rk) :: M_half    !< Final interface Mach number
    real(rk) :: p_half    !< Final interface pressure
    real(rk) :: a_half    !< Final interface sound speed
    real(rk) :: mass_flux !< Final mass flux across the interface

    integer(ik), parameter :: BOTTOM_IDX = 1 !< edge index for the bottom edge of the current cell
    integer(ik), parameter ::  RIGHT_IDX = 2 !< edge index for the right edge of the current cell
    integer(ik), parameter ::    TOP_IDX = 3 !< edge index for the top edge of the current cell
    integer(ik), parameter ::   LEFT_IDX = 4 !< edge index for the left edge of the current cell

    real(rk) :: n_x  !< normal vectors of each face
    real(rk) :: n_y  !< normal vectors of each face

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

    call debug_print('Running ausm_plus_solver_t%solve_ausm_plus()', __FILE__, __LINE__)

    if(dt < tiny(1.0_rk)) then
      write(std_err, '(a, es16.6)') "Error in ausm_plus_solver_t%solve_ausm_plus(), the timestep dt is < tiny(1.0_rk): dt = ", dt
      error stop "Error in ausm_plus_solver_t%solve_ausm_plus(), the timestep dt is < tiny(1.0_rk)"
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

    ! Get the edge values (left/right states for each edge)
    edge_interpolator => edge_interpolator_factory(self%input)
    call edge_interpolator%interpolate_edge_values(q=rho, lbounds=lbounds, edge_values=rho_interface_values)
    call edge_interpolator%interpolate_edge_values(q=u, lbounds=lbounds, edge_values=u_interface_values)
    call edge_interpolator%interpolate_edge_values(q=v, lbounds=lbounds, edge_values=v_interface_values)
    call edge_interpolator%interpolate_edge_values(q=p, lbounds=lbounds, edge_values=p_interface_values)

    call self%apply_reconstructed_bc(recon_rho=rho_interface_values, recon_u=u_interface_values, &
                                     recon_v=v_interface_values, recon_p=p_interface_values, &
                                     lbounds=lbounds, &
                                     bc_plus_x=bc_plus_x, bc_minus_x=bc_minus_x, &
                                     bc_plus_y=bc_plus_y, bc_minus_y=bc_minus_y)

    ! The interfaces where the flux occurs is shared by 2 cells. The indexing for the
    ! edges is staggered so that when reading/writing data to those arrays, we avoid using i-1
    ! as much as possible to avoid vectorization dependence flags from the compiler
    !
    !             (i-1/2)  (i+1/2)         interface locations
    !        |       |       |       |
    !        |       |     L | R     |     L/R locations
    !        |       |       |       |
    !        | (i-1) |  (i)  | (i+1) |     cell indexing
    !              (i-1)    (i)    (i+1)   edge indexing
    !        |   0   |   1   |   2   |     cell index example
    !        |       0       1       2     edge indexing example
    !        |       |       |       |
    !

    ilo = grid%ilo_cell
    ihi = grid%ihi_cell
    jlo = grid%jlo_cell
    jhi = grid%jhi_cell

    if(allocated(self%iflux)) deallocate(self%iflux)
    if(allocated(self%jflux)) deallocate(self%jflux)

    allocate(self%iflux(4, ilo - 1:ihi, jlo:jhi))
    allocate(self%jflux(4, ilo:ihi, jlo - 1:jhi))

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi), &
    !$omp shared(self, grid), &
    !$omp shared(rho_interface_values, u_interface_values, v_interface_values, p_interface_values), &
    !$omp private(i, j, n_x, n_y, mass_flux), &
    !$omp private(rho_L, u_L, v_L, p_L), &
    !$omp private(rho_R, u_R, v_R, p_R), &
    !$omp private(a_half, p_half, M_half, H)

    ! Find fluxes in the i-direction
    !$omp do
    do j = jlo, jhi       ! the iflux is the same size in the j direction as the cell quantities
      do i = ilo - 1, ihi ! start at ilo-1 b/c of the first i edge, e.g. i = 0

        ! use the normal vector of the right edge of the current cell
        n_x = grid%cell_edge_norm_vectors(1, RIGHT_IDX, i, j)
        n_y = grid%cell_edge_norm_vectors(2, RIGHT_IDX, i, j)

        ! right edge (i + 1/2)
        rho_L = rho_interface_values(RIGHT_IDX, i, j)    ! right edge of the current cell
        rho_R = rho_interface_values(LEFT_IDX, i + 1, j) ! left edge of the cell to the right
        u_L = u_interface_values(RIGHT_IDX, i, j)
        u_R = u_interface_values(LEFT_IDX, i + 1, j)
        v_L = v_interface_values(RIGHT_IDX, i, j)
        v_R = v_interface_values(LEFT_IDX, i + 1, j)
        p_L = p_interface_values(RIGHT_IDX, i, j)
        p_R = p_interface_values(LEFT_IDX, i + 1, j)

        call self%interface_state(rho=[rho_L, rho_R], &
                                  u=[u_L, u_R], &
                                  v=[v_L, v_R], &
                                  p=[p_L, p_R], &
                                  n=[n_x, n_y], &
                                  M_half=M_half, &
                                  p_half=p_half, &
                                  a_half=a_half, H=H)

        ! Edge flux values, via upwinding
        if(M_half > 0.0_rk) then
          mass_flux = M_half * a_half * rho_L
          self%iflux(1, i, j) = mass_flux
          self%iflux(2, i, j) = mass_flux * u_L + (n_x * p_half)
          self%iflux(3, i, j) = mass_flux * v_L + (n_y * p_half)
          self%iflux(4, i, j) = mass_flux * H(1)
        else
          mass_flux = M_half * a_half * rho_R
          self%iflux(1, i, j) = mass_flux
          self%iflux(2, i, j) = mass_flux * u_R + (n_x * p_half)
          self%iflux(3, i, j) = mass_flux * v_R + (n_y * p_half)
          self%iflux(4, i, j) = mass_flux * H(2)
        end if
      end do
    end do
    !$omp end do

    !$omp do
    ! Find fluxes in the j-direction
    do j = jlo - 1, jhi
      do i = ilo, ihi
        n_x = grid%cell_edge_norm_vectors(1, TOP_IDX, i, j)
        n_y = grid%cell_edge_norm_vectors(2, TOP_IDX, i, j)

        ! top edge (j + 1/2)
        rho_L = rho_interface_values(TOP_IDX, i, j)        ! top edge of the current cell
        rho_R = rho_interface_values(BOTTOM_IDX, i, j + 1) ! bottom edge of the cell above
        u_L = u_interface_values(TOP_IDX, i, j)
        u_R = u_interface_values(BOTTOM_IDX, i, j + 1)
        v_L = v_interface_values(TOP_IDX, i, j)
        v_R = v_interface_values(BOTTOM_IDX, i, j + 1)
        p_L = p_interface_values(TOP_IDX, i, j)
        p_R = p_interface_values(BOTTOM_IDX, i, j + 1)

        call self%interface_state(rho=[rho_L, rho_R], &
                                  u=[u_L, u_R], &
                                  v=[v_L, v_R], &
                                  p=[p_L, p_R], &
                                  n=[n_x, n_y], &
                                  M_half=M_half, &
                                  p_half=p_half, &
                                  a_half=a_half, H=H)

        ! Edge flux values, via upwinding
        if(M_half > 0.0_rk) then
          mass_flux = M_half * a_half * rho_L
          self%jflux(1, i, j) = mass_flux
          self%jflux(2, i, j) = mass_flux * u_L + (n_x * p_half)
          self%jflux(3, i, j) = mass_flux * v_L + (n_y * p_half)
          self%jflux(4, i, j) = mass_flux * H(1)
        else
          mass_flux = M_half * a_half * rho_R
          self%jflux(1, i, j) = mass_flux
          self%jflux(2, i, j) = mass_flux * u_R + (n_x * p_half)
          self%jflux(3, i, j) = mass_flux * v_R + (n_y * p_half)
          self%jflux(4, i, j) = mass_flux * H(2)
        end if

      end do
    end do
    !$omp end do
    !$omp end parallel

    call self%flux_split_edges(grid, lbounds, d_rho_dt, d_rho_u_dt, d_rho_v_dt, d_rho_E_dt)

    deallocate(rho_interface_values)
    deallocate(u_interface_values)
    deallocate(v_interface_values)
    deallocate(p_interface_values)
    deallocate(edge_interpolator)

    deallocate(bc_plus_x)
    deallocate(bc_plus_y)
    deallocate(bc_minus_x)
    deallocate(bc_minus_y)

    if(allocated(self%iflux)) deallocate(self%iflux)
    if(allocated(self%jflux)) deallocate(self%jflux)
  end subroutine solve_ausm_plus

  pure function split_mach_deg_1(M, plus_or_minus) result(M_split)
    !< Implementation of the 1st order polynomial split Mach function M±(1). See Eq. 18 in Ref [1]
    real(rk), intent(in) :: M         !< Mach number
    character(len=1), intent(in) :: plus_or_minus !< which split? '+' or '-'
    real(rk) :: M_split !< Split Mach number

    if(plus_or_minus == '+') then
      M_split = 0.5_rk * (M + abs(M)) ! M+(1)
    else
      M_split = 0.5_rk * (M - abs(M)) ! M-(1)
    end if
  end function split_mach_deg_1

  pure function split_mach_deg_2(M, plus_or_minus) result(M_split)
    !< Implementation of the 2nd order polynomial split Mach function  M±(2). See Eq. 19 in Ref [1]
    real(rk), intent(in) :: M                     !< Mach number
    character(len=1), intent(in) :: plus_or_minus !< which split? '+' or '-'
    real(rk) :: M_split                           !< Split Mach number, M±(2)

    if(plus_or_minus == '+') then
      M_split = 0.25_rk * (M + 1.0_rk)**2  ! M+(2)
    else
      M_split = -0.25_rk * (M - 1.0_rk)**2 ! M-(2)
    end if

  end function split_mach_deg_2

  pure function split_mach_deg_4(M, beta, plus_or_minus) result(M_split)
    !< Implementation of the 4th order polynomial split Mach function M±(4).
    !< See Eq. 20 in Ref [1]

    real(rk), intent(in) :: M         !< Mach number
    real(rk), intent(in) :: beta      !< Mach number
    character(len=1), intent(in) :: plus_or_minus !< which split? '+' or '-'
    real(rk) :: M_split !< Split Mach number

    ! Locals
    real(rk) :: M_2_plus  !< 2nd order split Mach polynomial M+(2)
    real(rk) :: M_2_minus !< 2nd order split Mach polynomial M-(2)

    M_2_plus = split_mach_deg_2(M, '+')  ! M+(2)
    M_2_minus = split_mach_deg_2(M, '-') ! M-(2)

    if(plus_or_minus == '+') then
      if(abs(M) >= 1.0_rk) then
        M_split = split_mach_deg_1(M, '+') ! M+(1)
      else
        ! M+(4)
        M_split = M_2_plus * (1.0_rk - 16.0_rk * beta * M_2_minus)
      end if
    else ! '-'
      if(abs(M) >= 1.0_rk) then
        M_split = split_mach_deg_1(M, '-')
      else
        ! M-(4)
        M_split = M_2_minus * (1.0_rk + 16.0_rk * beta * M_2_plus)
      end if
    end if
  end function split_mach_deg_4

  pure function split_pressure_deg_5(M, plus_or_minus, alpha) result(P_split)
    !< Implementation of the 5th order polynomial split pressure function P±(5).
    !< See Eq. 24 in Ref [1]

    real(rk), intent(in) :: M                     !< Mach number
    real(rk), intent(in) :: alpha
    character(len=1), intent(in) :: plus_or_minus !< Which split? '+' or '-'
    real(rk) :: P_split                           !< Split pressure factor

    ! Locals
    real(rk) :: M_2_plus  !< 2nd order split Mach polynomial M+(2)
    real(rk) :: M_2_minus !< 2nd order split Mach polynomial M-(2)

    M_2_plus = split_mach_deg_2(M, '+')  ! M+(2)
    M_2_minus = split_mach_deg_2(M, '-') ! M-(2)

    if(plus_or_minus == '+') then
      ! P+(5)
      if(abs(M) >= 1.0_rk) then
        P_split = split_mach_deg_1(M, '+') / M
      else
        P_split = M_2_plus * ((2.0_rk - M) - 16.0_rk * alpha * M * M_2_minus)
      end if
    else
      ! P-(5)
      if(abs(M) >= 1.0_rk) then
        P_split = split_mach_deg_1(M, '-') / M
      else
        P_split = M_2_minus * ((-2.0_rk - M) + 16.0_rk * alpha * M * M_2_plus)
      end if
    end if
  end function split_pressure_deg_5

  pure real(rk) function scaling_factor(M_0) result(f_a)
    !< Scaling function, Eq. 72
    real(rk), intent(in) :: M_0 ! Reference Mach number
    f_a = M_0 * (2.0_rk - M_0)
  end function scaling_factor

  pure real(rk) function alpha(f_a)
    !< alpha parameter, Eq. 76
    real(rk), intent(in) :: f_a ! Scaling factor -> [0,1]
    alpha = (3.0_rk / 16.0_rk) * (-4.0_rk + 5.0_rk * f_a**2)
  end function alpha

  pure subroutine interface_state(self, n, rho, u, v, p, M_half, p_half, a_half, H)
    !< Interface Mach number, e.g. M_(1/2)
    class(ausm_plus_solver_t), intent(in) :: self
    real(rk), dimension(2), intent(in) :: n   !< (x,y); edge normal vector
    real(rk), dimension(2), intent(in) :: rho !< (L,R); density
    real(rk), dimension(2), intent(in) :: u   !< (L,R); x-velocity
    real(rk), dimension(2), intent(in) :: v   !< (L,R); y-velocity
    real(rk), dimension(2), intent(in) :: p   !< (L,R); pressure
    real(rk), intent(out) :: M_half           !< Mach of the interface
    real(rk), intent(out) :: p_half           !< pressure of the interface
    real(rk), intent(out) :: a_half         !< total enthalpy of the left side
    real(rk), dimension(2), intent(out) :: H  !< (L,R); Total enthalpy

    ! Locals
    real(rk) :: vel_L          !< total-velocity of the left side
    real(rk) :: vel_R          !< total-velocity of the right side
    real(rk) :: M_p            !< Pressure diffusion term
    real(rk) :: M_plus         !< Split Mach term M+, see Eq. 73
    real(rk) :: M_minus        !< Split Mach term M-, see Eq. 73
    real(rk) :: M_sum          !< M_plus + M_minus; used for round-off machine epsilon checks
    real(rk) :: M_bar_sq       !< Mean Mach number ^2
    real(rk) :: M_0            !< Reference Mach
    real(rk) :: M_0_sq         !< Reference Mach squared
    real(rk) :: M_L            !< Mach number of the left side of the interface
    real(rk) :: M_R            !< Mach number of the right side of the interface
    real(rk) :: rho_half       !< Final interface density
    real(rk) :: alpha_param    !< Alpha parameter
    real(rk) :: f_a            !< Scaling factor: f_a(M_0) = M0(2-M0) -> [0,1]
    real(rk) :: p_u            !< Velocity difference (diffusion) term
    real(rk) :: P_plus         !< Split pressure term P+, see Eq. 75
    real(rk) :: P_minus        !< Split pressure term P-, see Eq. 75
    ! real(rk) :: P_sum          !< P_plus + P_minus; used for round-off machine epsilon checks
    real(rk) :: a_circumflex_L !< â_L; used to find the left interface sound speed see Eq 28 in Ref [1]
    real(rk) :: a_circumflex_R !< â_R; used to find the rgight interface sound speed see Eq 28 in Ref [1]
    real(rk) :: a_crit_L       !< a*; critical sound speed of the left side
    real(rk) :: a_crit_R       !< a*; critical sound speed of the right state
    real(rk) :: gamma          !< EOS polytropic index
    real(rk) :: gamma_factor   !< Precomputed gamma factor = 2*(gamma-1)/(gamma+1)

    ! Precompute a few re-used scalars
    gamma = eos%get_gamma()
    gamma_factor = ((2.0_rk * (gamma - 1.0_rk)) / (gamma + 1.0_rk))
    f_a = 0.0_rk

    associate(n_x=>n(1), n_y=>n(2), &
              rho_L=>rho(1), rho_R=>rho(2), &
              u_L=>u(1), u_R=>u(2), &
              v_L=>v(1), v_R=>v(2), &
              p_L=>p(1), p_R=>p(2), &
              H_L=>H(1), H_R=>H(2))
      ! Dot-product to get total velocity based on the edge normals
      vel_L = u_L * n_x + v_L * n_y
      vel_R = u_R * n_x + v_R * n_y

      ! Find total enthalpy in order to find the critical sound speed of the interface
      call eos%total_enthalpy(rho=rho_L, u=u_L, v=v_L, p=p_L, H=H_L)
      call eos%total_enthalpy(rho=rho_R, u=u_R, v=v_R, p=p_R, H=H_R)

      ! Critical sound speed, Eq. 29
      a_crit_L = sqrt(gamma_factor * abs(H_L)) ! a*_L
      a_crit_R = sqrt(gamma_factor * abs(H_R)) ! a*_R
    end associate

    ! Eq. 30
    a_circumflex_L = a_crit_L**2 / max(a_crit_L, vel_L)  ! â_L
    a_circumflex_R = a_crit_R**2 / max(a_crit_R, -vel_R) ! â_R

    ! Interface sound speed, Eq. 28
    a_half = min(a_circumflex_L, a_circumflex_R)

    ! Left/Right Mach, Eq. 69
    M_L = vel_L / a_half
    M_R = vel_R / a_half

    if(self%ausm_plus_up_all_speed) then! AUSM+-up all-speed scheme
      ! Mean local Mach, Eq. 70
      M_bar_sq = (vel_L**2 + vel_R**2) / (2.0_rk * a_half**2)

      ! Reference Mach number, see Eq. 71
      M_0_sq = min(1.0_rk, max(M_bar_sq, self%M_inf_sq))
      M_0 = sqrt(M_0_sq)
      f_a = scaling_factor(M_0) ! Scaling function, Eq. 72
      alpha_param = alpha(f_a) ! alpha parameter, Eq. 76
    else
      ! Mean local Mach, Eq. 13
      M_bar_sq = 0.5_rk * (M_L + M_R)
      alpha_param = self%pressure_split_alpha ! alpha parameter, Eq. 76
    end if

    ! Interface Mach number, Eq. 73
    M_plus = split_mach_deg_4(M=M_L, plus_or_minus='+', beta=self%mach_split_beta)
    M_minus = split_mach_deg_4(M=M_R, plus_or_minus='-', beta=self%mach_split_beta)

    ! Interface density
    rho_half = 0.5_rk * (rho(2) + rho(1))

    ! The pressure diffusion term, Mp, as defined in Eq. 21 in Ref [1]
    ! This is the "p" in AUSM+-up
    associate(K_p=>self%pressure_diffusion_coeff, sigma=>self%sigma, &
              p_L=>p(1), p_R=>p(2))

      if(abs(p_R - p_L) < epsilon(1.0_rk)) then
        M_p = 0.0_rk
      else if(self%ausm_plus_up_basic) then ! AUSM+-up basic scheme
        M_p = -K_p * max(1.0_rk - sigma * M_bar_sq, 0.0_rk) * ((p_R - p_L) / (rho_half * a_half**2))
      else if(self%ausm_plus_up_all_speed) then! AUSM+-up all-speed scheme
        M_p = -(K_p / f_a) * max(1.0_rk - sigma * M_bar_sq, 0.0_rk) * ((p_R - p_L) / (rho_half * a_half**2))
      else
        M_p = 0.0_rk
      end if
    end associate

    ! Interface Mach number
    M_sum = M_plus + M_minus
    if(abs(M_sum) < 1e-13_rk) M_sum = 0.0_rk

    M_half = M_sum + M_p
    if(abs(M_half) < 1e-13_rk) M_half = 0.0_rk

    ! Pressure splitting functions in Eq. 75
    P_plus = split_pressure_deg_5(M=M_L, plus_or_minus='+', alpha=alpha_param)
    P_minus = split_pressure_deg_5(M=M_R, plus_or_minus='-', alpha=alpha_param)

    ! Velocity difference (diffusion) term p_u, NOT, Eq. 26, but rather the last term in Eq 75
    ! This is the "u" in the AUSM+-u and AUSM+-up schemes
    associate(K_u=>self%pressure_flux_coeff, rho_L=>rho(1), rho_R=>rho(2))

      if(abs(vel_R - vel_L) < epsilon(1.0_rk)) then
        p_u = 0.0_rk
      else if(self%ausm_plus_up_basic) then ! AUSM+-up basic scheme
        p_u = -K_u * P_plus * P_minus * (rho_L + rho_R) * a_half * (vel_R - vel_L)
      else if(self%ausm_plus_up_all_speed) then! AUSM+-up all-speed scheme
        p_u = -K_u * P_plus * P_minus * (rho_L + rho_R) * (f_a * a_half) * (vel_R - vel_L)
      else
        p_u = 0.0_rk
      end if
    end associate

    ! Interface pressure, Eq. 75
    associate(p_L=>p(1), p_R=>p(2))
      p_half = P_plus * p_L + P_minus * p_R + p_u
    end associate

    ! write(*, '(10(es16.6))') M_half, p_half, a_half, H
  end subroutine interface_state

  subroutine finalize(self)
    !< Class finalizer
    type(ausm_plus_solver_t), intent(inout) :: self

    call debug_print('Running ausm_plus_solver_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%iflux)) deallocate(self%iflux) ! these should already be deallocated
    if(allocated(self%jflux)) deallocate(self%jflux) ! these should already be deallocated

    ! if(allocated(self%reconstructor)) deallocate(self%reconstructor)
    ! if(allocated(self%bc_plus_x)) deallocate(self%bc_plus_x)
    ! if(allocated(self%bc_plus_y)) deallocate(self%bc_plus_y)
    ! if(allocated(self%bc_minus_x)) deallocate(self%bc_minus_x)
    ! if(allocated(self%bc_minus_y)) deallocate(self%bc_minus_y)
  end subroutine finalize
end module mod_ausm_plus_solver
