module mod_fvleg_fluid
  !> Summary: Provide the basis for the fluid class that uses the FVLEG solver
  !> Date: 06/22/2020
  !> Author: Sam Miller
  !> Notes:
  !> References:
  !      [1]

  use mod_globals, only: debug_print, print_evolved_cell_data, print_recon_data
  use mod_nondimensionalization, only: scale_factors_set, rho_0, v_0, p_0, e_0, t_0
  use mod_base_fluid, only: base_fluid_t

  implicit none

  private
  public :: base_fluid_t, new_fluid

  logical, parameter :: filter_small_mach = .false.

  type, extends(base_fluid_t) :: fvleg_fluid_t
  contains
    procedure :: force_finalization
    procedure, public :: t => time_derivative
    procedure, nopass, private :: flux_edges
    final :: finalize
  end type fvleg_fluid_t

contains

  ! function time_derivative(self, fv, stage) result(d_dt)
  !   !< Implementation of dU/dt

  !   class(fvleg_fluid_t), intent(in) :: self
  !   class(finite_volume_scheme_t), intent(inout) :: fv
  !   integer(ik), intent(in) :: stage !< which stage in the time integration scheme are we in, e.g. RK2 stage 1

  !   ! Locals
  !   class(integrand_t), allocatable :: d_dt !< dU/dt (integrand_t to satisfy parent interface)
  !   type(fvleg_fluid_t), allocatable :: local_d_dt !< dU/dt
  !   integer(ik) :: error_code

  !   real(rk), dimension(:, :), allocatable :: evolved_corner_rho !< (i,j); Reconstructed rho at the corners
  !   real(rk), dimension(:, :), allocatable :: evolved_corner_u   !< (i,j); Reconstructed u at the corners
  !   real(rk), dimension(:, :), allocatable :: evolved_corner_v   !< (i,j); Reconstructed v at the corners
  !   real(rk), dimension(:, :), allocatable :: evolved_corner_p   !< (i,j); Reconstructed p at the corners
  !   real(rk), dimension(:, :), allocatable :: evolved_lr_mid_rho !< (i,j); Reconstructed rho at the left/right midpoints
  !   real(rk), dimension(:, :), allocatable :: evolved_lr_mid_u   !< (i,j); Reconstructed u at the left/right midpoints
  !   real(rk), dimension(:, :), allocatable :: evolved_lr_mid_v   !< (i,j); Reconstructed v at the left/right midpoints
  !   real(rk), dimension(:, :), allocatable :: evolved_lr_mid_p   !< (i,j); Reconstructed p at the left/right midpoints
  !   real(rk), dimension(:, :), allocatable :: evolved_du_mid_rho !< (i,j); Reconstructed rho at the down/up midpoints
  !   real(rk), dimension(:, :), allocatable :: evolved_du_mid_u   !< (i,j); Reconstructed u at the down/up midpoints
  !   real(rk), dimension(:, :), allocatable :: evolved_du_mid_v   !< (i,j); Reconstructed v at the down/up midpoints
  !   real(rk), dimension(:, :), allocatable :: evolved_du_mid_p   !< (i,j); Reconstructed p at the down/up midpoints

  !   real(rk), dimension(:, :, :), allocatable, target :: rho_recon_state
  !   !< ((corner1:midpoint4), i, j); reconstructed density, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint
  !   real(rk), dimension(:, :, :), allocatable, target :: u_recon_state
  !   !< ((corner1:midpoint4), i, j); reconstructed x-velocity, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint
  !   real(rk), dimension(:, :, :), allocatable, target :: v_recon_state
  !   !< ((corner1:midpoint4), i, j); reconstructed y-velocity, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint
  !   real(rk), dimension(:, :, :), allocatable, target :: p_recon_state
  !   !< ((corner1:midpoint4), i, j); reconstructed pressure, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint

  !   integer(ik), dimension(2) :: bounds

  !   call debug_print('Running fluid_t%time_derivative()', __FILE__, __LINE__)

  !   associate(imin=>fv%grid%ilo_bc_cell, imax=>fv%grid%ihi_bc_cell, &
  !             jmin=>fv%grid%jlo_bc_cell, jmax=>fv%grid%jhi_bc_cell)

  !     allocate(rho_recon_state(1:8, imin:imax, jmin:jmax))
  !     allocate(u_recon_state(1:8, imin:imax, jmin:jmax))
  !     allocate(v_recon_state(1:8, imin:imax, jmin:jmax))
  !     allocate(p_recon_state(1:8, imin:imax, jmin:jmax))
  !   end associate

  !   associate(imin_node=>fv%grid%ilo_node, imax_node=>fv%grid%ihi_node, &
  !             jmin_node=>fv%grid%jlo_node, jmax_node=>fv%grid%jhi_node, &
  !             imin_cell=>fv%grid%ilo_cell, imax_cell=>fv%grid%ihi_cell, &
  !             jmin_cell=>fv%grid%jlo_cell, jmax_cell=>fv%grid%jhi_cell)

  !     allocate(evolved_corner_rho(imin_node:imax_node, jmin_node:jmax_node))
  !     allocate(evolved_corner_u(imin_node:imax_node, jmin_node:jmax_node))
  !     allocate(evolved_corner_v(imin_node:imax_node, jmin_node:jmax_node))
  !     allocate(evolved_corner_p(imin_node:imax_node, jmin_node:jmax_node))

  !     allocate(evolved_lr_mid_rho(imin_cell:imax_cell, jmin_node:jmax_node))
  !     allocate(evolved_lr_mid_u(imin_cell:imax_cell, jmin_node:jmax_node))
  !     allocate(evolved_lr_mid_v(imin_cell:imax_cell, jmin_node:jmax_node))
  !     allocate(evolved_lr_mid_p(imin_cell:imax_cell, jmin_node:jmax_node))

  !     allocate(evolved_du_mid_rho(imin_node:imax_node, jmin_cell:jmax_cell))
  !     allocate(evolved_du_mid_u(imin_node:imax_node, jmin_cell:jmax_cell))
  !     allocate(evolved_du_mid_v(imin_node:imax_node, jmin_cell:jmax_cell))
  !     allocate(evolved_du_mid_p(imin_node:imax_node, jmin_cell:jmax_cell))

  !   end associate

  !   allocate(local_d_dt, source=self)

  !   if(.not. local_d_dt%prim_vars_updated) then
  !     call local_d_dt%calculate_derived_quantities()
  !   end if

  !   bounds = lbound(local_d_dt%rho)

  !   call fv%reconstruction_operator%set_cell_average_pointers(rho=local_d_dt%rho, &
  !                                                             p=local_d_dt%p, &
  !                                                             lbounds=bounds)
  !   call fv%apply_primitive_vars_bc(rho=local_d_dt%rho, &
  !                                   u=local_d_dt%u, &
  !                                   v=local_d_dt%v, &
  !                                   p=local_d_dt%p, lbounds=bounds)

  !   ! Now we can reconstruct the entire domain
  !   call debug_print('Reconstructing density', __FILE__, __LINE__)
  !   call fv%reconstruct(primitive_var=local_d_dt%rho, lbounds=bounds, &
  !                       reconstructed_var=rho_recon_state, name='rho', stage=stage)

  !   call debug_print('Reconstructing x-velocity', __FILE__, __LINE__)
  !   call fv%reconstruct(primitive_var=local_d_dt%u, lbounds=bounds, &
  !                       reconstructed_var=u_recon_state, name='u', stage=stage)

  !   call debug_print('Reconstructing y-velocity', __FILE__, __LINE__)
  !   call fv%reconstruct(primitive_var=local_d_dt%v, lbounds=bounds, &
  !                       reconstructed_var=v_recon_state, name='v', stage=stage)

  !   call debug_print('Reconstructing pressure', __FILE__, __LINE__)
  !   call fv%reconstruct(primitive_var=local_d_dt%p, lbounds=bounds, &
  !                       reconstructed_var=p_recon_state, name='p', stage=stage)

  !   ! The gradients have to be applied to the boundaries as well, since reconstructing
  !   ! at P'(x,y) requires the cell gradient
  !   call fv%apply_gradient_bc()

  !   ! Apply the reconstructed state to the ghost layers
  !   call fv%apply_reconstructed_state_bc(recon_rho=rho_recon_state, recon_u=u_recon_state, &
  !                                        recon_v=v_recon_state, recon_p=p_recon_state, lbounds=bounds)

  !   call fv%evolution_operator%set_reconstructed_state_pointers(rho=rho_recon_state, &
  !                                                               u=u_recon_state, &
  !                                                               v=v_recon_state, &
  !                                                               p=p_recon_state, &
  !                                                               lbounds=bounds)

  !   ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of down/up edge vectors
  !   call debug_print('Evolving down/up midpoints', __FILE__, __LINE__)
  !   bounds = lbound(evolved_du_mid_rho)
  !   call fv%evolution_operator%evolve(evolved_rho=evolved_du_mid_rho, evolved_u=evolved_du_mid_u, &
  !                                     evolved_v=evolved_du_mid_v, evolved_p=evolved_du_mid_p, &
  !                                     location='down/up midpoint', &
  !                                     lbounds=bounds, &
  !                                     error_code=error_code)

  !   ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of left/right edge vectors
  !   call debug_print('Evolving left/right midpoints', __FILE__, __LINE__)
  !   bounds = lbound(evolved_lr_mid_rho)
  !   call fv%evolution_operator%evolve(evolved_rho=evolved_lr_mid_rho, evolved_u=evolved_lr_mid_u, &
  !                                     evolved_v=evolved_lr_mid_v, evolved_p=evolved_lr_mid_p, &
  !                                     location='left/right midpoint', &
  !                                     lbounds=bounds, &
  !                                     error_code=error_code)

  !   ! Evolve, i.e. E0(R_omega), at all corner nodes
  !   call debug_print('Evolving corner nodes', __FILE__, __LINE__)
  !   bounds = lbound(evolved_corner_rho)
  !   call fv%evolution_operator%evolve(evolved_rho=evolved_corner_rho, evolved_u=evolved_corner_u, &
  !                                     evolved_v=evolved_corner_v, evolved_p=evolved_corner_p, &
  !                                     location='corner', &
  !                                     lbounds=bounds, &
  !                                     error_code=error_code)

  !   if(error_code /= 0) then
  !     fv%error_code = error_code
  !   end if

  !   nullify(fv%evolution_operator%reconstructed_rho)
  !   nullify(fv%evolution_operator%reconstructed_u)
  !   nullify(fv%evolution_operator%reconstructed_v)
  !   nullify(fv%evolution_operator%reconstructed_p)

  !   nullify(fv%reconstruction_operator%rho)
  !   nullify(fv%reconstruction_operator%p)

  !   call self%flux_edges(grid=fv%grid, &
  !                        evolved_corner_rho=evolved_corner_rho, &
  !                        evolved_corner_u=evolved_corner_u, &
  !                        evolved_corner_v=evolved_corner_v, &
  !                        evolved_corner_p=evolved_corner_p, &
  !                        evolved_lr_mid_rho=evolved_lr_mid_rho, &
  !                        evolved_lr_mid_u=evolved_lr_mid_u, &
  !                        evolved_lr_mid_v=evolved_lr_mid_v, &
  !                        evolved_lr_mid_p=evolved_lr_mid_p, &
  !                        evolved_du_mid_rho=evolved_du_mid_rho, &
  !                        evolved_du_mid_u=evolved_du_mid_u, &
  !                        evolved_du_mid_v=evolved_du_mid_v, &
  !                        evolved_du_mid_p=evolved_du_mid_p, &
  !                        d_rho_dt=local_d_dt%rho, &
  !                        d_rhou_dt=local_d_dt%rho_u, &
  !                        d_rhov_dt=local_d_dt%rho_v, &
  !                        d_rhoE_dt=local_d_dt%rho_E)

  !   call move_alloc(local_d_dt, d_dt)

  !   call d_dt%set_temp(calling_function='fluid_t%time_derivative (d_dt)', line=__LINE__)

  !   ! Now deallocate everything
  !   deallocate(evolved_corner_rho)
  !   deallocate(evolved_corner_u)
  !   deallocate(evolved_corner_v)
  !   deallocate(evolved_corner_p)

  !   deallocate(evolved_lr_mid_rho)
  !   deallocate(evolved_lr_mid_u)
  !   deallocate(evolved_lr_mid_v)
  !   deallocate(evolved_lr_mid_p)

  !   deallocate(evolved_du_mid_rho)
  !   deallocate(evolved_du_mid_u)
  !   deallocate(evolved_du_mid_v)
  !   deallocate(evolved_du_mid_p)

  !   deallocate(rho_recon_state)
  !   deallocate(u_recon_state)
  !   deallocate(v_recon_state)
  !   deallocate(p_recon_state)

  ! end function time_derivative

  ! subroutine flux_edges(grid, &
  !                       evolved_corner_rho, evolved_corner_u, evolved_corner_v, evolved_corner_p, &
  !                       evolved_lr_mid_rho, evolved_lr_mid_u, evolved_lr_mid_v, evolved_lr_mid_p, &
  !                       evolved_du_mid_rho, evolved_du_mid_u, evolved_du_mid_v, evolved_du_mid_p, &
  !                       d_rho_dt, d_rhou_dt, d_rhov_dt, d_rhoE_dt)
  !   !< Evaluate the fluxes along the edges. This is equation 13 in the paper
  !   ! FIXME: Move this into a separate class/module - this will allow easier unit testing
  !   class(grid_t), intent(in) :: grid

  !   real(rk), dimension(grid%ilo_node:, grid%jlo_node:), contiguous, intent(in) :: evolved_corner_rho !< (i,j); Reconstructed rho at the corners
  !   real(rk), dimension(grid%ilo_node:, grid%jlo_node:), contiguous, intent(in) :: evolved_corner_u   !< (i,j); Reconstructed u at the corners
  !   real(rk), dimension(grid%ilo_node:, grid%jlo_node:), contiguous, intent(in) :: evolved_corner_v   !< (i,j); Reconstructed v at the corners
  !   real(rk), dimension(grid%ilo_node:, grid%jlo_node:), contiguous, intent(in) :: evolved_corner_p   !< (i,j); Reconstructed p at the corners
  !   real(rk), dimension(grid%ilo_cell:, grid%jlo_node:), contiguous, intent(in) :: evolved_lr_mid_rho !< (i,j); Reconstructed rho at the left/right midpoints
  !   real(rk), dimension(grid%ilo_cell:, grid%jlo_node:), contiguous, intent(in) :: evolved_lr_mid_u   !< (i,j); Reconstructed u at the left/right midpoints
  !   real(rk), dimension(grid%ilo_cell:, grid%jlo_node:), contiguous, intent(in) :: evolved_lr_mid_v   !< (i,j); Reconstructed v at the left/right midpoints
  !   real(rk), dimension(grid%ilo_cell:, grid%jlo_node:), contiguous, intent(in) :: evolved_lr_mid_p   !< (i,j); Reconstructed p at the left/right midpoints
  !   real(rk), dimension(grid%ilo_node:, grid%jlo_cell:), contiguous, intent(in) :: evolved_du_mid_rho !< (i,j); Reconstructed rho at the down/up midpoints
  !   real(rk), dimension(grid%ilo_node:, grid%jlo_cell:), contiguous, intent(in) :: evolved_du_mid_u   !< (i,j); Reconstructed u at the down/up midpoints
  !   real(rk), dimension(grid%ilo_node:, grid%jlo_cell:), contiguous, intent(in) :: evolved_du_mid_v   !< (i,j); Reconstructed v at the down/up midpoints
  !   real(rk), dimension(grid%ilo_node:, grid%jlo_cell:), contiguous, intent(in) :: evolved_du_mid_p   !< (i,j); Reconstructed p at the down/up midpoints
  !   real(rk), dimension(grid%ilo_bc_cell:, grid%jlo_bc_cell:), contiguous, intent(inout) :: d_rho_dt
  !   real(rk), dimension(grid%ilo_bc_cell:, grid%jlo_bc_cell:), contiguous, intent(inout) :: d_rhou_dt
  !   real(rk), dimension(grid%ilo_bc_cell:, grid%jlo_bc_cell:), contiguous, intent(inout) :: d_rhov_dt
  !   real(rk), dimension(grid%ilo_bc_cell:, grid%jlo_bc_cell:), contiguous, intent(inout) :: d_rhoE_dt

  !   integer(ik) :: ilo, ihi, jlo, jhi
  !   integer(ik) :: i, j, k, edge, xy
  !   real(rk), dimension(2, 4) :: n_hat
  !   real(rk), dimension(4) :: delta_l
  !   integer(ik), dimension(2) :: bounds

  !   type(flux_array_t) :: corner_fluxes
  !   type(flux_array_t) :: downup_mid_fluxes
  !   type(flux_array_t) :: leftright_mid_fluxes

  !   real(rk) :: f_sum, g_sum
  !   real(rk) :: diff

  !   real(rk), dimension(4) :: bottom_flux
  !   real(rk), dimension(4) :: right_flux
  !   real(rk), dimension(4) :: top_flux
  !   real(rk), dimension(4) :: left_flux

  !   real(rk) :: rho_flux, ave_rho_flux
  !   real(rk) :: rhou_flux, ave_rhou_flux
  !   real(rk) :: rhov_flux, ave_rhov_flux
  !   real(rk) :: rhoE_flux, ave_rhoE_flux
  !   real(rk) :: threshold

  !   real(rk), parameter :: FLUX_EPS = 1e-13_rk !epsilon(1.0_rk)
  !   real(rk), parameter :: REL_THRESHOLD = 1e-5_rk

  !   call debug_print('Running fluid_t%flux_edges()', __FILE__, __LINE__)

  !   ! Get the flux arrays for each corner node or midpoint
  !   ilo = grid%ilo_node; ihi = grid%ihi_node
  !   jlo = grid%jlo_node; jhi = grid%jhi_node
  !   bounds = [ilo, jlo]
  !   corner_fluxes = get_fluxes(rho=evolved_corner_rho, u=evolved_corner_u, v=evolved_corner_v, &
  !                              p=evolved_corner_p, lbounds=bounds)

  !   ilo = grid%ilo_cell; ihi = grid%ihi_cell
  !   jlo = grid%jlo_node; jhi = grid%jhi_node
  !   bounds = [ilo, jlo]
  !   leftright_mid_fluxes = get_fluxes(rho=evolved_lr_mid_rho, u=evolved_lr_mid_u, v=evolved_lr_mid_v, &
  !                                     p=evolved_lr_mid_p, lbounds=bounds)

  !   ilo = grid%ilo_node; ihi = grid%ihi_node
  !   jlo = grid%jlo_cell; jhi = grid%jhi_cell
  !   bounds = [ilo, jlo]
  !   downup_mid_fluxes = get_fluxes(rho=evolved_du_mid_rho, u=evolved_du_mid_u, v=evolved_du_mid_v, &
  !                                  p=evolved_du_mid_p, lbounds=bounds)

  !   ilo = grid%ilo_cell; ihi = grid%ihi_cell
  !   jlo = grid%jlo_cell; jhi = grid%jhi_cell
  !   !$omp parallel default(none), &
  !   !$omp firstprivate(ilo, ihi, jlo, jhi) &
  !   !$omp private(i, j, k, delta_l, n_hat) &
  !   !$omp private(bottom_flux, right_flux, left_flux, top_flux) &
  !   !$omp private(f_sum, g_sum, diff) &
  !   !$omp shared(grid, corner_fluxes, leftright_mid_fluxes, downup_mid_fluxes) &
  !   !$omp shared(d_rho_dt, d_rhou_dt, d_rhov_dt, d_rhoE_dt) &
  !   !$omp private(rho_flux, ave_rho_flux, rhou_flux, ave_rhou_flux, rhov_flux, ave_rhov_flux, rhoE_flux, ave_rhoE_flux, threshold)
  !   !$omp do simd
  !   do j = jlo, jhi
  !     do i = ilo, ihi

  !       delta_l = grid%cell_edge_lengths(:, i, j)

  !       do edge = 1, 4
  !         do xy = 1, 2
  !           n_hat(xy, edge) = grid%cell_edge_norm_vectors(xy, edge, i, j)
  !         end do
  !       end do

  !       associate(F_c1=>corner_fluxes%F(:, i, j), &
  !                 G_c1=>corner_fluxes%G(:, i, j), &
  !                 F_c2=>corner_fluxes%F(:, i + 1, j), &
  !                 G_c2=>corner_fluxes%G(:, i + 1, j), &
  !                 F_c3=>corner_fluxes%F(:, i + 1, j + 1), &
  !                 G_c3=>corner_fluxes%G(:, i + 1, j + 1), &
  !                 F_c4=>corner_fluxes%F(:, i, j + 1), &
  !                 G_c4=>corner_fluxes%G(:, i, j + 1), &
  !                 F_m1=>leftright_mid_fluxes%F(:, i, j), &
  !                 G_m1=>leftright_mid_fluxes%G(:, i, j), &
  !                 F_m2=>downup_mid_fluxes%F(:, i + 1, j), &
  !                 G_m2=>downup_mid_fluxes%G(:, i + 1, j), &
  !                 F_m3=>leftright_mid_fluxes%F(:, i, j + 1), &
  !                 G_m3=>leftright_mid_fluxes%G(:, i, j + 1), &
  !                 G_m4=>downup_mid_fluxes%G(:, i, j), &
  !                 F_m4=>downup_mid_fluxes%F(:, i, j), &
  !                 n_hat_1=>grid%cell_edge_norm_vectors(:, 1, i, j), &
  !                 n_hat_2=>grid%cell_edge_norm_vectors(:, 2, i, j), &
  !                 n_hat_3=>grid%cell_edge_norm_vectors(:, 3, i, j), &
  !                 n_hat_4=>grid%cell_edge_norm_vectors(:, 4, i, j))

  !         !  Bottom
  !         ! bottom_flux (rho, rhou, rhov, rhoE)
  !         do k = 1, 4
  !           f_sum = sum([F_c1(k), 4.0_rk * F_m1(k), F_c2(k)]) * n_hat_1(1)
  !           g_sum = sum([G_c1(k), 4.0_rk * G_m1(k), G_c2(k)]) * n_hat_1(2)
  !           ! if(abs(f_sum + g_sum) < epsilon(1.0_rk)) then !FIXME: is this the right way to do it?
  !           !   bottom_flux(k) = 0.0_rk
  !           ! else
  !           bottom_flux(k) = (f_sum + g_sum) * (delta_l(1) / 6.0_rk)
  !           ! end if
  !         end do

  !         !  Right
  !         ! right_flux (rho, rhou, rhov, rhoE)
  !         do k = 1, 4
  !           f_sum = sum([F_c2(k), 4.0_rk * F_m2(k), F_c3(k)]) * n_hat_2(1)
  !           g_sum = sum([G_c2(k), 4.0_rk * G_m2(k), G_c3(k)]) * n_hat_2(2)
  !           ! if(abs(f_sum + g_sum) < epsilon(1.0_rk)) then
  !           !   right_flux(k) = 0.0_rk
  !           ! else
  !           right_flux(k) = (f_sum + g_sum) * (delta_l(2) / 6.0_rk)
  !           ! end if
  !         end do

  !         ! Top
  !         ! top_flux (rho, rhou, rhov, rhoE)
  !         do k = 1, 4
  !           f_sum = sum([F_c3(k), 4.0_rk * F_m3(k), F_c4(k)]) * n_hat_3(1)
  !           g_sum = sum([G_c3(k), 4.0_rk * G_m3(k), G_c4(k)]) * n_hat_3(2)
  !           ! if(abs(f_sum + g_sum) < epsilon(1.0_rk)) then
  !           !   top_flux(k) = 0.0_rk
  !           ! else
  !           top_flux(k) = (f_sum + g_sum) * (delta_l(3) / 6.0_rk)
  !           ! end if
  !         end do

  !         ! Left
  !         ! left_flux (rho, rhou, rhov, rhoE)
  !         do k = 1, 4
  !           f_sum = sum([F_c4(k), 4.0_rk * F_m4(k), F_c1(k)]) * n_hat_4(1)
  !           g_sum = sum([G_c4(k), 4.0_rk * G_m4(k), G_c1(k)]) * n_hat_4(2)
  !           ! if(abs(f_sum + g_sum) < epsilon(1.0_rk)) then
  !           ! left_flux(k) = 0.0_rk
  !           ! else
  !           left_flux(k) = (f_sum + g_sum) * (delta_l(4) / 6.0_rk)
  !           ! end if
  !         end do

  !         ave_rho_flux = sum([abs(left_flux(1)), abs(right_flux(1)), abs(top_flux(1)), abs(bottom_flux(1))]) / 4.0_rk
  !         ave_rhou_flux = sum([abs(left_flux(2)), abs(right_flux(2)), abs(top_flux(2)), abs(bottom_flux(2))]) / 4.0_rk
  !         ave_rhov_flux = sum([abs(left_flux(3)), abs(right_flux(3)), abs(top_flux(3)), abs(bottom_flux(3))]) / 4.0_rk
  !         ave_rhoE_flux = sum([abs(left_flux(4)), abs(right_flux(4)), abs(top_flux(4)), abs(bottom_flux(4))]) / 4.0_rk

  !         ! Now sum them all up together
  !         rho_flux = neumaier_sum_4([left_flux(1), right_flux(1), top_flux(1), bottom_flux(1)])
  !         rhou_flux = neumaier_sum_4([left_flux(2), right_flux(2), top_flux(2), bottom_flux(2)])
  !         rhov_flux = neumaier_sum_4([left_flux(3), right_flux(3), top_flux(3), bottom_flux(3)])
  !         rhoE_flux = neumaier_sum_4([left_flux(4), right_flux(4), top_flux(4), bottom_flux(4)])

  !         threshold = abs(ave_rho_flux) * REL_THRESHOLD
  !         if(abs(rho_flux) < threshold .or. abs(rho_flux) < epsilon(1.0_rk)) then
  !           rho_flux = 0.0_rk
  !         end if

  !         threshold = abs(ave_rhou_flux) * REL_THRESHOLD
  !         if(abs(rhou_flux) < threshold .or. abs(rhou_flux) < epsilon(1.0_rk)) then
  !           rhou_flux = 0.0_rk
  !         end if

  !         threshold = abs(ave_rhov_flux) * REL_THRESHOLD
  !         if(abs(rhov_flux) < threshold .or. abs(rhov_flux) < epsilon(1.0_rk)) then
  !           rhov_flux = 0.0_rk
  !         end if

  !         threshold = abs(ave_rhoE_flux) * REL_THRESHOLD
  !         if(abs(rhoE_flux) < threshold .or. abs(rhoE_flux) < epsilon(1.0_rk)) then
  !           rhoE_flux = 0.0_rk
  !         end if

  !         d_rho_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * rho_flux
  !         d_rhou_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * rhou_flux
  !         d_rhov_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * rhov_flux
  !         d_rhoE_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * rhoE_flux

  !         ! if(abs(rhov_flux) > 0.0_rk) then
  !         !   print *, 'i,j', i, j
  !         !   write(*, '(a, es16.6)') 'rhov_flux:                        : ', rhov_flux
  !         !   ! write(*, '(a, es16.6)') 'FLUX_EPS                          : ', FLUX_EPS
  !         !   ! write(*, '(a, es16.6)') 'ave_rhov_flux                     : ', ave_rhov_flux
  !         !   ! write(*, '(a, es16.6)') 'REL_THRESHOLD                     : ', REL_THRESHOLD
  !         !   ! write(*, '(a, es16.6)') 'abs(ave_rhov_flux) * REL_THRESHOLD: ', abs(ave_rhov_flux) * REL_THRESHOLD
  !         !   ! write(*, '(a, 4(es16.6, 1x))') ' -> ', left_flux(3), right_flux(3), top_flux(3), bottom_flux(3)

  !         !   print*, 'bottom flux'
  !         !   print*, 'n_hat_1: ', n_hat_1
  !         !   write(*, '(a, 4(es16.6))') "F_c1:          ", F_c1
  !         !   write(*, '(a, 4(es16.6))') "4.0_rk * F_m1: ", 4.0_rk * F_m1
  !         !   write(*, '(a, 4(es16.6))') "F_c2:          ", F_c2
  !         !   write(*, '(a, 4(es16.6))') "G_c1:          ", G_c1
  !         !   write(*, '(a, 4(es16.6))') "4.0_rk * G_m1: ", 4.0_rk * G_m1
  !         !   write(*, '(a, 4(es16.6))') "G_c2:          ", G_c2
  !         !   print*, '(delta_l(1) / 6.0_rk) ', (delta_l(1) / 6.0_rk)
  !         !   do k = 1, 4
  !         !     f_sum = sum([F_c1(k), 4.0_rk * F_m1(k), F_c2(k)]) * n_hat_1(1)
  !         !     g_sum = sum([G_c1(k), 4.0_rk * G_m1(k), G_c2(k)]) * n_hat_1(2)
  !         !     print*, 'f_sum', k, f_sum
  !         !     print*, 'g_sum', k, g_sum
  !         !   end do

  !         !   error stop
  !         ! end if

  !       end associate

  !     end do ! i
  !   end do ! j
  !   !$omp end do simd
  !   !$omp end parallel

  !   ilo = grid%ilo_bc_cell; ihi = grid%ihi_bc_cell
  !   jlo = grid%jlo_bc_cell; jhi = grid%jhi_bc_cell
  !   d_rho_dt(ilo, :) = 0.0_rk
  !   d_rho_dt(ihi, :) = 0.0_rk
  !   d_rho_dt(:, jlo) = 0.0_rk
  !   d_rho_dt(:, jhi) = 0.0_rk

  !   d_rhou_dt(ilo, :) = 0.0_rk
  !   d_rhou_dt(ihi, :) = 0.0_rk
  !   d_rhou_dt(:, jlo) = 0.0_rk
  !   d_rhou_dt(:, jhi) = 0.0_rk

  !   d_rhov_dt(ilo, :) = 0.0_rk
  !   d_rhov_dt(ihi, :) = 0.0_rk
  !   d_rhov_dt(:, jlo) = 0.0_rk
  !   d_rhov_dt(:, jhi) = 0.0_rk

  !   d_rhoE_dt(ilo, :) = 0.0_rk
  !   d_rhoE_dt(ihi, :) = 0.0_rk
  !   d_rhoE_dt(:, jlo) = 0.0_rk
  !   d_rhoE_dt(:, jhi) = 0.0_rk

  !   ! print *, 'fluxed vars differences'
  !   ! write(*, '(a, 6(es16.6))') "d/dt rho   : ", maxval(d_rho_dt(204, jlo + 1:jhi - 1)) - minval(d_rho_dt(204, jlo + 1:jhi - 1))
  !   ! write(*, '(a, 6(es16.6))') "d/dt rho u : ", maxval(d_rhou_dt(204, jlo + 1:jhi - 1)) - minval(d_rhou_dt(204, jlo + 1:jhi - 1))
  !   ! write(*, '(a, 6(es16.6))') "d/dt rho v : ", maxval(d_rhov_dt(204, jlo + 1:jhi - 1)) - minval(d_rhov_dt(204, jlo + 1:jhi - 1))
  !   ! write(*, '(a, 6(es16.6))') "d/dt rho E : ", maxval(d_rhoE_dt(204, jlo + 1:jhi - 1)) - minval(d_rhoE_dt(204, jlo + 1:jhi - 1))
  !   ! print *

  ! end subroutine flux_edges

  subroutine force_finalization(self)
    type(fvleg_fluid_t), intent(inout) :: self
    call self%finalize()
  end subroutine force_finalization

  subroutine finalize(self)
    type(fvleg_fluid_t), intent(inout) :: self

    call debug_print('Running fluid_t%finalize()', __FILE__, __LINE__)
    if(allocated(self%rho)) deallocate(self%rho)
    if(allocated(self%u)) deallocate(self%u)
    if(allocated(self%v)) deallocate(self%v)
    if(allocated(self%p)) deallocate(self%p)
    if(allocated(self%rho_u)) deallocate(self%rho_u)
    if(allocated(self%rho_v)) deallocate(self%rho_v)
    if(allocated(self%rho_E)) deallocate(self%rho_E)
    if(allocated(self%cs)) deallocate(self%cs)
    if(allocated(self%mach_u)) deallocate(self%mach_u)
    if(allocated(self%mach_v)) deallocate(self%mach_v)
    if(allocated(self%time_integrator)) deallocate(self%time_integrator)
  end subroutine finalize

end module mod_fvleg_fluid
