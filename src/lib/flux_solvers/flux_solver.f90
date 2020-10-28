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
  use mod_field, only: field_2d_t
  use mod_globals, only: debug_print
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

  end type flux_solver_t

  type, abstract, extends(flux_solver_t) :: edge_split_flux_solver_t
    !< Directionally split flux solver class
    real(rk), dimension(:, :, :), allocatable :: iflux !< ((1:4), i, j) edge flux of the i-direction edges
    real(rk), dimension(:, :, :), allocatable :: jflux !< ((1:4), i, j) edge flux of the j-direction edges
  contains
    procedure, public :: flux_split_edges
  end type edge_split_flux_solver_t

  abstract interface
    subroutine solve(self, grid, rho, u, v, p, d_rho_dt, d_rho_u_dt, d_rho_v_dt, d_rho_E_dt)
      !< Solve and flux the edges
      import :: flux_solver_t
      import :: grid_block_2d_t, field_2d_t
      class(flux_solver_t), intent(inout) :: self
      class(grid_block_2d_t), intent(in) :: grid
      class(field_2d_t), intent(inout) :: rho !< density
      class(field_2d_t), intent(inout) :: u   !< x-velocity
      class(field_2d_t), intent(inout) :: v   !< y-velocity
      class(field_2d_t), intent(inout) :: p   !< pressure
      type(field_2d_t), intent(out)   ::   d_rho_dt !< d/dt of the density field
      type(field_2d_t), intent(out)   :: d_rho_u_dt !< d/dt of the rhou field
      type(field_2d_t), intent(out)   :: d_rho_v_dt !< d/dt of the rhov field
      type(field_2d_t), intent(out)   :: d_rho_E_dt !< d/dt of the rhoE field
    end subroutine

    subroutine initialize(self, input, time)
      import :: flux_solver_t, rk
      import :: input_t
      class(flux_solver_t), intent(inout) :: self
      class(input_t), intent(in) :: input
      real(rk), intent(in) :: time
    end subroutine
  end interface

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

  end subroutine init_boundary_conditions

  subroutine apply_primitive_bc(self, rho, u, v, p, &
                                bc_plus_x, bc_minus_x, bc_plus_y, bc_minus_y)
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

    call debug_print('Running flux_solver_t%apply()', __FILE__, __LINE__)

    max_priority_bc = max(bc_plus_x%priority, bc_plus_y%priority, &
                          bc_minus_x%priority, bc_minus_y%priority)

    call bc_plus_x%set_time(time=self%time)
    call bc_minus_x%set_time(time=self%time)
    call bc_plus_y%set_time(time=self%time)
    call bc_minus_y%set_time(time=self%time)

    do priority = max_priority_bc, 0, -1

      if(bc_plus_x%priority == priority) then
        call bc_plus_x%apply(rho=rho, u=u, v=v, p=p)
      end if

      if(bc_plus_y%priority == priority) then
        call bc_plus_y%apply(rho=rho, u=u, v=v, p=p)
      end if

      if(bc_minus_x%priority == priority) then
        call bc_minus_x%apply(rho=rho, u=u, v=v, p=p)
      end if

      if(bc_minus_y%priority == priority) then
        call bc_minus_y%apply(rho=rho, u=u, v=v, p=p)
      end if

    end do

  end subroutine apply_primitive_bc

  subroutine flux_split_edges(self, grid, d_rho_dt, d_rho_u_dt, d_rho_v_dt, d_rho_E_dt)
    !< Flux the edges to get the residuals, e.g. 1/vol * d/dt U
    class(edge_split_flux_solver_t), intent(in) :: self
    class(grid_block_2d_t), intent(in) :: grid          !< grid topology class
    type(field_2d_t), intent(out) ::   d_rho_dt  !< d/dt of the density field
    type(field_2d_t), intent(out) :: d_rho_u_dt  !< d/dt of the rhou field
    type(field_2d_t), intent(out) :: d_rho_v_dt  !< d/dt of the rhov field
    type(field_2d_t), intent(out) :: d_rho_E_dt  !< d/dt of the rhoE field

    ! Locals
    integer(ik) :: i, j, ilo, ihi, jlo, jhi
    real(rk), dimension(4) :: delta_l !< edge length
    real(rk) :: volume !< cell volume
    real(rk), parameter :: FLUX_EPS = 5e-13_rk
    real(rk), parameter :: REL_THRESHOLD = 1e-5_rk !< relative error threshold

    real(rk), dimension(4) :: rho_edge_fluxes
    real(rk), dimension(4) :: rhou_edge_fluxes
    real(rk), dimension(4) :: rhov_edge_fluxes
    real(rk), dimension(4) :: rhoE_edge_fluxes
    real(rk) :: ave_rho_edge_flux, rho_flux, rho_flux_threshold
    real(rk) :: ave_rhou_edge_flux, rhou_flux, rhou_flux_threshold
    real(rk) :: ave_rhov_edge_flux, rhov_flux, rhov_flux_threshold
    real(rk) :: ave_rhoE_edge_flux, rhoE_flux, rhoE_flux_threshold

    ! Block bounds
    ilo = self%lbounds(1)
    ihi = self%ubounds(1)
    jlo = self%lbounds(2)
    jhi = self%ubounds(2)

    ! d_rho_dt = field_2d(name='d_rho_dt  ', long_name='d/dt(rho)', descrip='', units='', &
    !                       dims=[grid%ni_cell, grid%nj_cell])
    ! d_rho_u_dt = field_2d(name='d_rho_u_dt', long_name='d/dt(rho u)', descrip='', units='', &
    !                         dims=[grid%ni_cell, grid%nj_cell])
    ! d_rho_v_dt = field_2d(name='d_rho_v_dt', long_name='d/dt(rho v)', descrip='', units='', &
    !                         dims=[grid%ni_cell, grid%nj_cell])
    ! d_rho_E_dt = field_2d(name='d_rho_E_dt', long_name='d/dt(rho E)', descrip='', units='', &
    !                         dims=[grid%ni_cell, grid%nj_cell])

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

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi), &
    !$omp shared(self, grid, d_rho_dt, d_rho_u_dt, d_rho_v_dt, d_rho_E_dt), &
    !$omp private(i, j, delta_l, volume) &
    !$omp private(ave_rho_edge_flux, rho_flux, rho_flux_threshold) &
    !$omp private(ave_rhou_edge_flux, rhou_flux, rhou_flux_threshold) &
    !$omp private(ave_rhov_edge_flux, rhov_flux, rhov_flux_threshold) &
    !$omp private(ave_rhoE_edge_flux, rhoE_flux, rhoE_flux_threshold) &
    !$omp private(rho_edge_fluxes, rhou_edge_fluxes, rhov_edge_fluxes, rhoE_edge_fluxes)
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi

        delta_l = grid%edge_lengths(:, i, j)
        volume = grid%volume(i, j)

        ! rho
        rho_edge_fluxes = [self%iflux(1, i, j) * delta_l(2), -self%iflux(1, i - 1, j) * delta_l(4), &
                           self%jflux(1, i, j) * delta_l(3), -self%jflux(1, i, j - 1) * delta_l(1)]
        ave_rho_edge_flux = 0.25_rk * sum(abs(rho_edge_fluxes))

        rho_flux = -neumaier_sum_4(rho_edge_fluxes)
        ! rho_flux = -sum(rho_edge_fluxes)

        rho_flux_threshold = abs(ave_rho_edge_flux) * REL_THRESHOLD
        if(abs(rho_flux) < rho_flux_threshold .or. abs(rho_flux) < epsilon(1.0_rk)) then
          rho_flux = 0.0_rk
        end if
        d_rho_dt%data(i, j) = rho_flux / volume

        ! rho u
        rhou_edge_fluxes = [self%iflux(2, i, j) * delta_l(2), -self%iflux(2, i - 1, j) * delta_l(4), &
                            self%jflux(2, i, j) * delta_l(3), -self%jflux(2, i, j - 1) * delta_l(1)]
        ave_rhou_edge_flux = 0.25_rk * sum(abs(rhou_edge_fluxes))

        rhou_flux = -neumaier_sum_4(rhou_edge_fluxes)
        ! rhou_flux = -sum(rhou_edge_fluxes)

        rhou_flux_threshold = abs(ave_rhou_edge_flux) * REL_THRESHOLD
        if(abs(rhou_flux) < rhou_flux_threshold .or. abs(rhou_flux) < epsilon(1.0_rk)) then
          rhou_flux = 0.0_rk
        end if
        d_rho_u_dt%data(i, j) = rhou_flux / volume

        ! rho v
        rhov_edge_fluxes = [self%iflux(3, i, j) * delta_l(2), -self%iflux(3, i - 1, j) * delta_l(4), &
                            self%jflux(3, i, j) * delta_l(3), -self%jflux(3, i, j - 1) * delta_l(1)]
        ave_rhov_edge_flux = 0.25_rk * sum(abs(rhov_edge_fluxes))

        rhov_flux = -neumaier_sum_4(rhov_edge_fluxes)
        ! rhov_flux = -sum(rhov_edge_fluxes)

        rhov_flux_threshold = abs(ave_rhov_edge_flux) * REL_THRESHOLD
        if(abs(rhov_flux) < rhov_flux_threshold .or. abs(rhov_flux) < epsilon(1.0_rk)) then
          rhov_flux = 0.0_rk
        end if
        d_rho_v_dt%data(i, j) = rhov_flux / volume

        ! rho E
        rhoE_edge_fluxes = [self%iflux(4, i, j) * delta_l(2), -self%iflux(4, i - 1, j) * delta_l(4), &
                            self%jflux(4, i, j) * delta_l(3), -self%jflux(4, i, j - 1) * delta_l(1)]
        ave_rhoE_edge_flux = 0.25_rk * sum(abs(rhoE_edge_fluxes))

        rhoE_flux = -neumaier_sum_4(rhoE_edge_fluxes)
        ! rhoE_flux = -sum(rhoE_edge_fluxes)

        rhoE_flux_threshold = abs(ave_rhoE_edge_flux) * REL_THRESHOLD
        if(abs(rhoE_flux) < rhoE_flux_threshold .or. abs(rhoE_flux) < epsilon(1.0_rk)) then
          rhoE_flux = 0.0_rk
        end if
        d_rho_E_dt%data(i, j) = rhoE_flux / volume
      end do
    end do
    !$omp end do
    !$omp end parallel

    ! Zero out the halo layers
    call d_rho_dt%zero_out_halo()
    call d_rho_u_dt%zero_out_halo()
    call d_rho_v_dt%zero_out_halo()
    call d_rho_E_dt%zero_out_halo()

  end subroutine flux_split_edges

end module mod_flux_solver
