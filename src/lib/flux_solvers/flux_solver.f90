! MIT License
! Copyright (c) 2020 Sam Miller
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

module mod_flux_solver
  !< Summary: Provide a base Riemann solver class structure
  !< Date: 06/22/2020
  !< Author: Sam Miller

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_field, only: field_2d_t, field_2d
  use mod_globals, only: debug_print, enable_debug_print
  use mod_floating_point_utils, only: neumaier_sum_4
  use mod_grid_block_2d, only: grid_block_2d_t
  use mod_input, only: input_t
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_bc_factory, only: bc_factory

  implicit none

  private
  public :: flux_solver_t, edge_split_flux_solver_t

  type, abstract :: flux_solver_t
    type(input_t) :: input         !<
    character(len=32) :: name = '' !<
    integer(ik) :: iteration = 0   !<
    real(rk) :: time = 0.0_rk      !<
    real(rk) :: dt = 0.0_rk        !<
    integer(ik), dimension(2) :: lbounds = 0 !< (i,j); lower cell bounds
    integer(ik), dimension(2) :: ubounds = 0 !< (i,j); upper cell bounds
  contains
    procedure, public :: init_boundary_conditions
    procedure, public :: apply_primitive_bc
    ! Deferred methods
    procedure(initialize), deferred, public :: initialize
    procedure(solve), deferred, public :: solve

  endtype flux_solver_t

  type, abstract, extends(flux_solver_t) :: edge_split_flux_solver_t
    !< Directionally split flux solver class
    real(rk), dimension(:, :, :), allocatable :: iflux !< ((1:4), i, j) edge flux of the i-direction edges
    real(rk), dimension(:, :, :), allocatable :: jflux !< ((1:4), i, j) edge flux of the j-direction edges
  contains
    procedure, public :: flux_split_edges
  endtype edge_split_flux_solver_t

  abstract interface
    subroutine solve(self, dt, grid, rho, u, v, p, d_rho_dt, d_rho_u_dt, d_rho_v_dt, d_rho_E_dt)
      !< Solve and flux the edges
      import :: flux_solver_t, rk
      import :: grid_block_2d_t, field_2d_t
      class(flux_solver_t), intent(inout) :: self
      real(rk), intent(in) :: dt !< timestep delta t
      class(grid_block_2d_t), intent(in) :: grid
      class(field_2d_t), intent(inout) :: rho !< density
      class(field_2d_t), intent(inout) :: u   !< x-velocity
      class(field_2d_t), intent(inout) :: v   !< y-velocity
      class(field_2d_t), intent(inout) :: p   !< pressure

      real(rk), dimension(:, :), allocatable, intent(out) ::   d_rho_dt  !< d/dt of the density field
      real(rk), dimension(:, :), allocatable, intent(out) :: d_rho_u_dt  !< d/dt of the rhou field
      real(rk), dimension(:, :), allocatable, intent(out) :: d_rho_v_dt  !< d/dt of the rhov field
      real(rk), dimension(:, :), allocatable, intent(out) :: d_rho_E_dt  !< d/dt of the rhoE field

      ! type(field_2d_t), intent(out)   ::   d_rho_dt !< d/dt of the density field
      ! type(field_2d_t), intent(out)   :: d_rho_u_dt !< d/dt of the rhou field
      ! type(field_2d_t), intent(out)   :: d_rho_v_dt !< d/dt of the rhov field
      ! type(field_2d_t), intent(out)   :: d_rho_E_dt !< d/dt of the rhoE field
    endsubroutine

    subroutine initialize(self, input, time)
      import :: flux_solver_t, rk
      import :: input_t
      class(flux_solver_t), intent(inout) :: self
      class(input_t), intent(in) :: input
      real(rk), intent(in) :: time
    endsubroutine
  endinterface

contains

  subroutine init_boundary_conditions(self, grid, bc_plus_x, bc_minus_x, bc_plus_y, bc_minus_y)
    class(flux_solver_t), intent(inout) :: self
    class(grid_block_2d_t), intent(in) :: grid
    class(boundary_condition_t), pointer :: bc => null()

    class(boundary_condition_t), allocatable, intent(out):: bc_plus_x
    class(boundary_condition_t), allocatable, intent(out):: bc_plus_y
    class(boundary_condition_t), allocatable, intent(out):: bc_minus_x
    class(boundary_condition_t), allocatable, intent(out):: bc_minus_y

    ! Locals
    integer(ik) :: alloc_status

    if(enable_debug_print) call debug_print('Calling flux_solver_t%init_boundary_conditions()', __FILE__, __LINE__)

    ! Set boundary conditions
    bc => bc_factory(bc_type=self%input%plus_x_bc, location='+x', input=self%input, grid=grid, time=self%time)
    allocate(bc_plus_x, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate bc_plus_x"
    deallocate(bc)

    bc => bc_factory(bc_type=self%input%plus_y_bc, location='+y', input=self%input, grid=grid, time=self%time)
    allocate(bc_plus_y, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate bc_plus_y"
    deallocate(bc)

    bc => bc_factory(bc_type=self%input%minus_x_bc, location='-x', input=self%input, grid=grid, time=self%time)
    allocate(bc_minus_x, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate bc_minus_x"
    deallocate(bc)

    bc => bc_factory(bc_type=self%input%minus_y_bc, location='-y', input=self%input, grid=grid, time=self%time)
    allocate(bc_minus_y, source=bc, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate bc_minus_y"
    deallocate(bc)

  endsubroutine init_boundary_conditions

  subroutine apply_primitive_bc(self, rho, u, v, p, &
                                bc_plus_x, bc_minus_x, bc_plus_y, bc_minus_y)
    !< Apply the boundary conditions to the primitive variables (rho, u, v, and p). BC's are applied
    !< in order of their priority, with the highest going first
    class(flux_solver_t), intent(inout) :: self
    class(field_2d_t), intent(inout) :: rho
    class(field_2d_t), intent(inout) :: u
    class(field_2d_t), intent(inout) :: v
    class(field_2d_t), intent(inout) :: p
    class(boundary_condition_t), intent(inout):: bc_plus_x
    class(boundary_condition_t), intent(inout):: bc_plus_y
    class(boundary_condition_t), intent(inout):: bc_minus_x
    class(boundary_condition_t), intent(inout):: bc_minus_y

    integer(ik) :: priority
    integer(ik) :: max_priority_bc !< highest goes first

    call debug_print('Running flux_solver_t%apply_primitive_bc()', __FILE__, __LINE__)

    max_priority_bc = max(bc_plus_x%priority, bc_plus_y%priority, &
                          bc_minus_x%priority, bc_minus_y%priority)

    call bc_plus_x%set_time(time=self%time)
    call bc_minus_x%set_time(time=self%time)
    call bc_plus_y%set_time(time=self%time)
    call bc_minus_y%set_time(time=self%time)

    do priority = max_priority_bc, 0, -1

      if(bc_plus_x%priority == priority) then
        call bc_plus_x%apply(rho=rho, u=u, v=v, p=p)
      endif

      if(bc_plus_y%priority == priority) then
        call bc_plus_y%apply(rho=rho, u=u, v=v, p=p)
      endif

      if(bc_minus_x%priority == priority) then
        call bc_minus_x%apply(rho=rho, u=u, v=v, p=p)
      endif

      if(bc_minus_y%priority == priority) then
        call bc_minus_y%apply(rho=rho, u=u, v=v, p=p)
      endif

    enddo

  endsubroutine apply_primitive_bc

  subroutine flux_split_edges(self, grid, d_rho_dt, d_rho_u_dt, d_rho_v_dt, d_rho_E_dt)
    !< Flux the edges to get the residuals, e.g. 1/vol * d/dt U
    class(edge_split_flux_solver_t), intent(in) :: self
    class(grid_block_2d_t), intent(in) :: grid          !< grid topology class
    real(rk), dimension(:, :), allocatable, intent(out) ::   d_rho_dt  !< d/dt of the density field
    real(rk), dimension(:, :), allocatable, intent(out) :: d_rho_u_dt  !< d/dt of the rhou field
    real(rk), dimension(:, :), allocatable, intent(out) :: d_rho_v_dt  !< d/dt of the rhov field
    real(rk), dimension(:, :), allocatable, intent(out) :: d_rho_E_dt  !< d/dt of the rhoE field

    ! Locals
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    integer(ik) :: ilo_halo, ihi_halo, jlo_halo, jhi_halo
    real(rk), dimension(4) :: delta_l !< edge length
    real(rk), parameter :: EPS = 1e-13_rk
    real(rk), parameter :: FLUX_EPS = 1e-13_rk

    real(rk) :: rho_i, rho_j
    real(rk) :: rhou_i, rhou_j
    real(rk) :: rhov_i, rhov_j
    real(rk) :: rhoE_i, rhoE_j
    real(rk) :: rho_flux 
    real(rk) :: rhou_flux
    real(rk) :: rhov_flux
    real(rk) :: rhoE_flux

    call debug_print('Running flux_solver_t%flux_split_edges()', __FILE__, __LINE__)

    ! Block bounds
    ilo = grid%lbounds(1)
    ihi = grid%ubounds(1)
    jlo = grid%lbounds(2)
    jhi = grid%ubounds(2)

    ilo_halo = grid%lbounds_halo(1)
    ihi_halo = grid%ubounds_halo(1)
    jlo_halo = grid%lbounds_halo(2)
    jhi_halo = grid%ubounds_halo(2)

    allocate(d_rho_dt(ilo_halo:ihi_halo, jlo_halo:jhi_halo))
    allocate(d_rho_u_dt(ilo_halo:ihi_halo, jlo_halo:jhi_halo))
    allocate(d_rho_v_dt(ilo_halo:ihi_halo, jlo_halo:jhi_halo))
    allocate(d_rho_E_dt(ilo_halo:ihi_halo, jlo_halo:jhi_halo))

    !                                   /\
    !                  jflux(i,j)  'R'  |
    !                o--------------------o
    !                |                'L' |
    !            <---|                    |--->
    ! -iflux(i-1, j) |     cell (i,j)     | iflux(i, j)
    !                |                    |
    !                |                'L' | 'R'
    !                o--------------------o
    !                   jflux(i,j-1)   |
    !                                 \/
    !
    ! This is the numbering convention that this module uses

    ! ------------------------------------
    ! rho
    ! ------------------------------------
    do j = jlo, jhi
      do i = ilo, ihi
        delta_l = grid%edge_lengths(:, i, j)
        rho_i = (self%iflux(1, i, j) * delta_l(2)) + (-self%iflux(1, i - 1, j) * delta_l(4))
        rho_j = (self%jflux(1, i, j) * delta_l(3)) + (-self%jflux(1, i, j - 1) * delta_l(1))
        if(abs(rho_i) < EPS) rho_i = 0.0_rk
        if(abs(rho_j) < EPS) rho_j = 0.0_rk
        rho_flux = rho_i + rho_j
        if(abs(rho_flux) < FLUX_EPS) rho_flux = 0.0_rk
        d_rho_dt(i, j) = -rho_flux
      enddo
    enddo

    ! ------------------------------------
    ! rho u
    ! ------------------------------------
    do j = jlo, jhi
      do i = ilo, ihi
        delta_l = grid%edge_lengths(:, i, j)
        rhou_i = (self%iflux(2, i, j) * delta_l(2)) + (-self%iflux(2, i - 1, j) * delta_l(4))
        rhou_j = (self%jflux(2, i, j) * delta_l(3)) + (-self%jflux(2, i, j - 1) * delta_l(1))
        if(abs(rhou_i) < EPS) rhou_i = 0.0_rk
        if(abs(rhou_j) < EPS) rhou_j = 0.0_rk
        rhou_flux = rhou_i + rhou_j
        if(abs(rhou_flux) < FLUX_EPS) rhou_flux = 0.0_rk
        d_rho_u_dt(i, j) = -rhou_flux
      enddo
    enddo

    ! ------------------------------------
    ! rho v
    ! ------------------------------------
    do j = jlo, jhi
      do i = ilo, ihi
        delta_l = grid%edge_lengths(:, i, j)
        rhov_i = (self%iflux(3, i, j) * delta_l(2)) + (-self%iflux(3, i - 1, j) * delta_l(4))
        rhov_j = (self%jflux(3, i, j) * delta_l(3)) + (-self%jflux(3, i, j - 1) * delta_l(1))
        if(abs(rhov_i) < EPS) rhov_i = 0.0_rk
        if(abs(rhov_j) < EPS) rhov_j = 0.0_rk
        rhov_flux = rhov_i + rhov_j
        if(abs(rhov_flux) < FLUX_EPS) rhov_flux = 0.0_rk
        d_rho_v_dt(i, j) = -rhov_flux
      enddo
    enddo

    ! ------------------------------------
    ! rho E
    ! ------------------------------------
    do j = jlo, jhi
      do i = ilo, ihi
        delta_l = grid%edge_lengths(:, i, j)
        rhoE_i = (self%iflux(4, i, j) * delta_l(2)) + (-self%iflux(4, i - 1, j) * delta_l(4))
        rhoE_j = (self%jflux(4, i, j) * delta_l(3)) + (-self%jflux(4, i, j - 1) * delta_l(1))
        if(abs(rhoE_i) < EPS) rhoE_i = 0.0_rk
        if(abs(rhoE_j) < EPS) rhoE_j = 0.0_rk
        rhoE_flux = rhoE_i + rhoE_j
        if(abs(rhoE_flux) < FLUX_EPS) rhoE_flux = 0.0_rk
        d_rho_E_dt(i, j) = -rhoE_flux
      enddo
    enddo
      
    !$omp simd
    do j = jlo, jhi
      do i = ilo, ihi
        d_rho_dt(i, j) = d_rho_dt(i, j) / grid%volume(i, j)
      enddo
    enddo
    !$omp end simd

    !$omp simd
    do j = jlo, jhi
      do i = ilo, ihi
        d_rho_u_dt(i, j) = d_rho_u_dt(i, j) / grid%volume(i, j)
      enddo
    enddo
    !$omp end simd

    !$omp simd
    do j = jlo, jhi
      do i = ilo, ihi
        d_rho_v_dt(i, j) = d_rho_v_dt(i, j) / grid%volume(i, j)
      enddo
    enddo
    !$omp end simd
    
    !$omp simd
    do j = jlo, jhi
      do i = ilo, ihi
        d_rho_E_dt(i, j) = d_rho_E_dt(i, j) / grid%volume(i, j)
      enddo
    enddo
    !$omp end simd

    ! Zero out the halo layers
    ! call d_rho_dt%zero_out_halo()
    ! call d_rho_u_dt%zero_out_halo()
    ! call d_rho_v_dt%zero_out_halo()
    ! call d_rho_E_dt%zero_out_halo()

    associate(ilo_s => grid%lbounds_halo(1), ilo_e => grid%lbounds(1) - 1, &
              ihi_s => grid%ubounds(1) + 1, ihi_e => grid%ubounds_halo(1), &
              jlo_s => grid%lbounds_halo(2), jlo_e => grid%lbounds(2) - 1, &
              jhi_s => grid%ubounds(2) + 1, jhi_e => grid%ubounds_halo(2))

      d_rho_dt(ilo_s:ilo_e, :) = 0.0_rk ! lower i cells
      d_rho_dt(ihi_s:ihi_e, :) = 0.0_rk ! upper i cells
      d_rho_dt(:, jlo_s:jlo_e) = 0.0_rk ! lower j cells
      d_rho_dt(:, jhi_s:jhi_e) = 0.0_rk ! upper j cells

      d_rho_u_dt(ilo_s:ilo_e, :) = 0.0_rk ! lower i cells
      d_rho_u_dt(ihi_s:ihi_e, :) = 0.0_rk ! upper i cells
      d_rho_u_dt(:, jlo_s:jlo_e) = 0.0_rk ! lower j cells
      d_rho_u_dt(:, jhi_s:jhi_e) = 0.0_rk ! upper j cells

      d_rho_v_dt(ilo_s:ilo_e, :) = 0.0_rk ! lower i cells
      d_rho_v_dt(ihi_s:ihi_e, :) = 0.0_rk ! upper i cells
      d_rho_v_dt(:, jlo_s:jlo_e) = 0.0_rk ! lower j cells
      d_rho_v_dt(:, jhi_s:jhi_e) = 0.0_rk ! upper j cells

      d_rho_E_dt(ilo_s:ilo_e, :) = 0.0_rk ! lower i cells
      d_rho_E_dt(ihi_s:ihi_e, :) = 0.0_rk ! upper i cells
      d_rho_E_dt(:, jlo_s:jlo_e) = 0.0_rk ! lower j cells
      d_rho_E_dt(:, jhi_s:jhi_e) = 0.0_rk ! upper j cells
    endassociate

  endsubroutine flux_split_edges

endmodule mod_flux_solver
