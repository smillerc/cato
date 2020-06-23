module mod_fvleg_solver
  !> Summary: Provide a FVLEG solver class structure
  !> Date: 06/22/2020
  !> Author: Sam Miller
  !> Notes:
  !> References:
  !      [1]

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_riemann_solver, only: riemann_solver_t
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_finite_volume_schemes, only: finite_volume_scheme_t
  use mod_grid, only: grid_t
  use mod_input, only: input_t

  implicit none

  private
  public :: riemann_solver_t

  type, extends(riemann_solver_t) :: fvleg_solver_t
    class(abstract_evo_operator_t), allocatable :: evolution_operator
  contains
    procedure, public :: initialize => initialize_fvleg
    procedure, public :: solve => solve_fvleg
    final :: finalize
  end type fvleg_solver_t

contains
  subroutine initialize_fvleg(self, grid, input)
    class(fvleg_solver_t), intent(in) :: self
    class(grid_t), intent(in) :: grid
    class(input_t), intent(in) :: input
  end subroutine initialize_fvleg

  subroutine solve_fvleg(self, time, grid, lbounds, rho, u, v, p, rho_u, rho_v, rho_E)
    class(fvleg_solver_t), intent(in) :: self
    real(rk), intent(in) :: time
    class(grid_t), intent(in) :: grid
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(inout) :: rho
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in) :: u
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in) :: v
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in) :: p
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: rho_u
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: rho_v
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: rho_E

    ! class(fluid_t), intent(in) :: self
    ! class(finite_volume_scheme_t), intent(inout) :: fv
    ! integer(ik), intent(in) :: stage !< which stage in the time integration scheme are we in, e.g. RK2 stage 1

    ! ! Locals
    ! class(integrand_t), allocatable :: d_dt !< dU/dt (integrand_t to satisfy parent interface)
    ! type(fluid_t), allocatable :: local_d_dt !< dU/dt
    ! integer(ik) :: error_code

    real(rk), dimension(:, :), allocatable :: evolved_corner_rho !< (i,j); Reconstructed rho at the corners
    real(rk), dimension(:, :), allocatable :: evolved_corner_u   !< (i,j); Reconstructed u at the corners
    real(rk), dimension(:, :), allocatable :: evolved_corner_v   !< (i,j); Reconstructed v at the corners
    real(rk), dimension(:, :), allocatable :: evolved_corner_p   !< (i,j); Reconstructed p at the corners
    real(rk), dimension(:, :), allocatable :: evolved_lr_mid_rho !< (i,j); Reconstructed rho at the left/right midpoints
    real(rk), dimension(:, :), allocatable :: evolved_lr_mid_u   !< (i,j); Reconstructed u at the left/right midpoints
    real(rk), dimension(:, :), allocatable :: evolved_lr_mid_v   !< (i,j); Reconstructed v at the left/right midpoints
    real(rk), dimension(:, :), allocatable :: evolved_lr_mid_p   !< (i,j); Reconstructed p at the left/right midpoints
    real(rk), dimension(:, :), allocatable :: evolved_du_mid_rho !< (i,j); Reconstructed rho at the down/up midpoints
    real(rk), dimension(:, :), allocatable :: evolved_du_mid_u   !< (i,j); Reconstructed u at the down/up midpoints
    real(rk), dimension(:, :), allocatable :: evolved_du_mid_v   !< (i,j); Reconstructed v at the down/up midpoints
    real(rk), dimension(:, :), allocatable :: evolved_du_mid_p   !< (i,j); Reconstructed p at the down/up midpoints

    real(rk), dimension(:, :, :), allocatable, target :: rho_recon_state
    !< ((corner1:midpoint4), i, j); reconstructed density, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint
    real(rk), dimension(:, :, :), allocatable, target :: u_recon_state
    !< ((corner1:midpoint4), i, j); reconstructed x-velocity, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint
    real(rk), dimension(:, :, :), allocatable, target :: v_recon_state
    !< ((corner1:midpoint4), i, j); reconstructed y-velocity, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint
    real(rk), dimension(:, :, :), allocatable, target :: p_recon_state
    !< ((corner1:midpoint4), i, j); reconstructed pressure, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint

    integer(ik), dimension(2) :: bounds

    call debug_print('Running fluid_t%time_derivative()', __FILE__, __LINE__)

    ! associate(imin=>fv%grid%ilo_bc_cell, imax=>fv%grid%ihi_bc_cell, &
    !           jmin=>fv%grid%jlo_bc_cell, jmax=>fv%grid%jhi_bc_cell)

    !   allocate(rho_recon_state(1:8, imin:imax, jmin:jmax))
    !   allocate(u_recon_state(1:8, imin:imax, jmin:jmax))
    !   allocate(v_recon_state(1:8, imin:imax, jmin:jmax))
    !   allocate(p_recon_state(1:8, imin:imax, jmin:jmax))
    ! end associate

    ! associate(imin_node=>fv%grid%ilo_node, imax_node=>fv%grid%ihi_node, &
    !           jmin_node=>fv%grid%jlo_node, jmax_node=>fv%grid%jhi_node, &
    !           imin_cell=>fv%grid%ilo_cell, imax_cell=>fv%grid%ihi_cell, &
    !           jmin_cell=>fv%grid%jlo_cell, jmax_cell=>fv%grid%jhi_cell)

    !   allocate(evolved_corner_rho(imin_node:imax_node, jmin_node:jmax_node))
    !   allocate(evolved_corner_u(imin_node:imax_node, jmin_node:jmax_node))
    !   allocate(evolved_corner_v(imin_node:imax_node, jmin_node:jmax_node))
    !   allocate(evolved_corner_p(imin_node:imax_node, jmin_node:jmax_node))

    !   allocate(evolved_lr_mid_rho(imin_cell:imax_cell, jmin_node:jmax_node))
    !   allocate(evolved_lr_mid_u(imin_cell:imax_cell, jmin_node:jmax_node))
    !   allocate(evolved_lr_mid_v(imin_cell:imax_cell, jmin_node:jmax_node))
    !   allocate(evolved_lr_mid_p(imin_cell:imax_cell, jmin_node:jmax_node))

    !   allocate(evolved_du_mid_rho(imin_node:imax_node, jmin_cell:jmax_cell))
    !   allocate(evolved_du_mid_u(imin_node:imax_node, jmin_cell:jmax_cell))
    !   allocate(evolved_du_mid_v(imin_node:imax_node, jmin_cell:jmax_cell))
    !   allocate(evolved_du_mid_p(imin_node:imax_node, jmin_cell:jmax_cell))

    ! end associate

    ! allocate(local_d_dt, source=self)

    ! if(.not. local_d_dt%prim_vars_updated) then
    !   call local_d_dt%calculate_derived_quantities()
    ! end if

    ! bounds = lbound(local_d_dt%rho)

    ! call fv%reconstruction_operator%set_cell_average_pointers(rho=local_d_dt%rho, &
    !                                                           p=local_d_dt%p, &
    !                                                           lbounds=bounds)
    ! call fv%apply_primitive_vars_bc(rho=local_d_dt%rho, &
    !                                 u=local_d_dt%u, &
    !                                 v=local_d_dt%v, &
    !                                 p=local_d_dt%p, lbounds=bounds)

    ! ! Now we can reconstruct the entire domain
    ! call debug_print('Reconstructing density', __FILE__, __LINE__)
    ! call fv%reconstruct(primitive_var=local_d_dt%rho, lbounds=bounds, &
    !                     reconstructed_var=rho_recon_state, name='rho', stage=stage)

    ! call debug_print('Reconstructing x-velocity', __FILE__, __LINE__)
    ! call fv%reconstruct(primitive_var=local_d_dt%u, lbounds=bounds, &
    !                     reconstructed_var=u_recon_state, name='u', stage=stage)

    ! call debug_print('Reconstructing y-velocity', __FILE__, __LINE__)
    ! call fv%reconstruct(primitive_var=local_d_dt%v, lbounds=bounds, &
    !                     reconstructed_var=v_recon_state, name='v', stage=stage)

    ! call debug_print('Reconstructing pressure', __FILE__, __LINE__)
    ! call fv%reconstruct(primitive_var=local_d_dt%p, lbounds=bounds, &
    !                     reconstructed_var=p_recon_state, name='p', stage=stage)

    ! ! The gradients have to be applied to the boundaries as well, since reconstructing
    ! ! at P'(x,y) requires the cell gradient
    ! call fv%apply_gradient_bc()

    ! ! Apply the reconstructed state to the ghost layers
    ! call fv%apply_reconstructed_state_bc(recon_rho=rho_recon_state, recon_u=u_recon_state, &
    !                                       recon_v=v_recon_state, recon_p=p_recon_state, lbounds=bounds)

    ! call fv%evolution_operator%set_reconstructed_state_pointers(rho=rho_recon_state, &
    !                                                             u=u_recon_state, &
    !                                                             v=v_recon_state, &
    !                                                             p=p_recon_state, &
    !                                                             lbounds=bounds)

    ! ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of down/up edge vectors
    ! call debug_print('Evolving down/up midpoints', __FILE__, __LINE__)
    ! bounds = lbound(evolved_du_mid_rho)
    ! call fv%evolution_operator%evolve(evolved_rho=evolved_du_mid_rho, evolved_u=evolved_du_mid_u, &
    !                                   evolved_v=evolved_du_mid_v, evolved_p=evolved_du_mid_p, &
    !                                   location='down/up midpoint', &
    !                                   lbounds=bounds, &
    !                                   error_code=error_code)

    ! ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of left/right edge vectors
    ! call debug_print('Evolving left/right midpoints', __FILE__, __LINE__)
    ! bounds = lbound(evolved_lr_mid_rho)
    ! call fv%evolution_operator%evolve(evolved_rho=evolved_lr_mid_rho, evolved_u=evolved_lr_mid_u, &
    !                                   evolved_v=evolved_lr_mid_v, evolved_p=evolved_lr_mid_p, &
    !                                   location='left/right midpoint', &
    !                                   lbounds=bounds, &
    !                                   error_code=error_code)

    ! ! Evolve, i.e. E0(R_omega), at all corner nodes
    ! call debug_print('Evolving corner nodes', __FILE__, __LINE__)
    ! bounds = lbound(evolved_corner_rho)
    ! call fv%evolution_operator%evolve(evolved_rho=evolved_corner_rho, evolved_u=evolved_corner_u, &
    !                                   evolved_v=evolved_corner_v, evolved_p=evolved_corner_p, &
    !                                   location='corner', &
    !                                   lbounds=bounds, &
    !                                   error_code=error_code)

    ! if(error_code /= 0) then
    !   fv%error_code = error_code
    ! end if

    ! nullify(fv%evolution_operator%reconstructed_rho)
    ! nullify(fv%evolution_operator%reconstructed_u)
    ! nullify(fv%evolution_operator%reconstructed_v)
    ! nullify(fv%evolution_operator%reconstructed_p)

    ! nullify(fv%reconstruction_operator%rho)
    ! nullify(fv%reconstruction_operator%p)

    ! call self%flux_edges(grid=fv%grid, &
    !                       evolved_corner_rho=evolved_corner_rho, &
    !                       evolved_corner_u=evolved_corner_u, &
    !                       evolved_corner_v=evolved_corner_v, &
    !                       evolved_corner_p=evolved_corner_p, &
    !                       evolved_lr_mid_rho=evolved_lr_mid_rho, &
    !                       evolved_lr_mid_u=evolved_lr_mid_u, &
    !                       evolved_lr_mid_v=evolved_lr_mid_v, &
    !                       evolved_lr_mid_p=evolved_lr_mid_p, &
    !                       evolved_du_mid_rho=evolved_du_mid_rho, &
    !                       evolved_du_mid_u=evolved_du_mid_u, &
    !                       evolved_du_mid_v=evolved_du_mid_v, &
    !                       evolved_du_mid_p=evolved_du_mid_p, &
    !                       d_rho_dt=local_d_dt%rho, &
    !                       d_rhou_dt=local_d_dt%rho_u, &
    !                       d_rhov_dt=local_d_dt%rho_v, &
    !                       d_rhoE_dt=local_d_dt%rho_E)

    ! call move_alloc(local_d_dt, d_dt)

    ! call d_dt%set_temp(calling_function='fluid_t%time_derivative (d_dt)', line=__LINE__)

    ! Now deallocate everything
    deallocate(evolved_corner_rho)
    deallocate(evolved_corner_u)
    deallocate(evolved_corner_v)
    deallocate(evolved_corner_p)

    deallocate(evolved_lr_mid_rho)
    deallocate(evolved_lr_mid_u)
    deallocate(evolved_lr_mid_v)
    deallocate(evolved_lr_mid_p)

    deallocate(evolved_du_mid_rho)
    deallocate(evolved_du_mid_u)
    deallocate(evolved_du_mid_v)
    deallocate(evolved_du_mid_p)

    deallocate(rho_recon_state)
    deallocate(u_recon_state)
    deallocate(v_recon_state)
    deallocate(p_recon_state)

  end subroutine solve_fvleg

  subroutine finalize(self)
    !< Class finalizer
    type(fvleg_solver_t), intent(inout) :: self
    if(allocated(self%reconstructor)) deallocate(self%reconstructor)
    if(allocated(self%evolution_operator)) deallocate(self%evolution_operator)
    if(allocated(self%bc_plus_x)) deallocate(self%bc_plus_x)
    if(allocated(self%bc_plus_y)) deallocate(self%bc_plus_y)
    if(allocated(self%bc_minus_x)) deallocate(self%bc_minus_x)
    if(allocated(self%bc_minus_y)) deallocate(self%bc_minus_y)
  end subroutine finalize
end module mod_fvleg_solver
