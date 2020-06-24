module mod_fvleg_solver
  !> Summary: Provide a FVLEG solver class structure
  !> Date: 06/22/2020
  !> Author: Sam Miller
  !> Notes:
  !> References:
  !      [1]

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_flux_solver, only: flux_solver_t
  use mod_floating_point_utils, only: neumaier_sum_4
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_local_evo_operator, only: local_evo_operator_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_reconstruction_factory, only: reconstruction_factory
  use mod_flux_array, only: get_fluxes, flux_array_t
  use mod_grid, only: grid_t
  use mod_input, only: input_t

  implicit none

  private
  public :: fvleg_solver_t

  type, extends(flux_solver_t) :: fvleg_solver_t
    !< Solver class that uses the FVLEG method to reconstruct, evolve, and flux quantities
    !< at the cell edge midpoints and corners
    ! class(abstract_evo_operator_t), allocatable :: evolution_operator !< E0 in the paper

    logical :: reconstruct_p_prime = .false.
    !< for the Mach cones, reconstruct the solution at P' or just use P0 from the cell that contains P'

  contains
    ! Public methods
    procedure, public :: initialize => initialize_fvleg
    procedure, public :: solve => solve_fvleg
    procedure, public, pass(lhs) :: copy => copy_fvleg

    ! Private methods
    procedure, private, nopass :: flux_edges
    procedure, private :: apply_primitive_bc
    procedure, private :: apply_gradient_bc
    procedure, private :: apply_reconstructed_bc
    final :: finalize

    ! Operators
    generic :: assignment(=) => copy
  end type fvleg_solver_t

contains
  subroutine initialize_fvleg(self, grid, input)
    !< Construct the FVLEG solver

    class(fvleg_solver_t), intent(inout) :: self
    class(grid_t), intent(in), target :: grid
    class(input_t), intent(in) :: input

    ! Locals

    integer(ik) :: alloc_status
    alloc_status = 0

    call debug_print('Running fvleg_solver_t%initialize_fvleg()', __FILE__, __LINE__)

    self%name = 'fvleg'
    self%input = input

  end subroutine initialize_fvleg

  subroutine copy_fvleg(lhs, rhs)
    !< Implement LHS = RHS
    class(fvleg_solver_t), intent(inout) :: lhs
    type(fvleg_solver_t), intent(in) :: rhs

    call debug_print('Running fvleg_solver_t%copy()', __FILE__, __LINE__)

    lhs%iteration = rhs%iteration
    lhs%time = rhs%time
    lhs%dt = rhs%dt
    lhs%reconstruct_p_prime = rhs%reconstruct_p_prime

  end subroutine copy_fvleg

  subroutine solve_fvleg(self, dt, grid, lbounds, rho, u, v, p, &
                         d_rho_dt, d_rho_u_dt, d_rho_v_dt, d_rho_E_dt)
    !< Solve and flux the edges
    class(fvleg_solver_t), intent(inout) :: self
    class(grid_t), intent(in) :: grid
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), intent(in) :: dt
    real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: rho !< density
    real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: u   !< x-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: v   !< y-velocity
    real(rk), dimension(lbounds(1):, lbounds(2):), target, intent(inout) :: p   !< pressure
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) ::   d_rho_dt    !< d/dt of the density field
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: d_rho_u_dt    !< d/dt of the rhou field
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: d_rho_v_dt    !< d/dt of the rhov field
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: d_rho_E_dt    !< d/dt of the rhoE field

    ! Locals
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

    integer(ik), dimension(2) :: evo_bounds
    integer(ik) :: error_code
    integer(ik) :: stage = 0
    character(len=50) :: stage_name = ''

    class(abstract_reconstruction_t), pointer :: reconstructor => null()
    type(local_evo_operator_t) :: E0 !< Local Evolution Operator, E0

    class(boundary_condition_t), allocatable:: bc_plus_x
    class(boundary_condition_t), allocatable:: bc_plus_y
    class(boundary_condition_t), allocatable:: bc_minus_x
    class(boundary_condition_t), allocatable:: bc_minus_y

    error_code = 0
    evo_bounds = 0

    if(dt < tiny(1.0_rk)) error stop "Error in fvleg_solver_t%solve_fvleg(), the timestep dt is < tiny(1.0_rk)"
    self%time = self%time + dt
    self%dt = dt
    self%iteration = self%iteration + 1

    reconstructor => reconstruction_factory(input=self%input, grid_target=grid)
    call reconstructor%set_grid_pointer(grid)

    call E0%initialize(grid_target=grid, dt=dt, recon_operator_target=reconstructor)

    call self%init_boundary_conditions(grid, &
                                       bc_plus_x=bc_plus_x, bc_minus_x=bc_minus_x, &
                                       bc_plus_y=bc_plus_y, bc_minus_y=bc_minus_y)

    call debug_print('Running fvleg_t%solve_fvleg()', __FILE__, __LINE__)

    associate(imin=>grid%ilo_bc_cell, imax=>grid%ihi_bc_cell, &
              jmin=>grid%jlo_bc_cell, jmax=>grid%jhi_bc_cell)

      allocate(rho_recon_state(1:8, imin:imax, jmin:jmax))
      allocate(u_recon_state(1:8, imin:imax, jmin:jmax))
      allocate(v_recon_state(1:8, imin:imax, jmin:jmax))
      allocate(p_recon_state(1:8, imin:imax, jmin:jmax))
    end associate

    associate(imin_node=>grid%ilo_node, imax_node=>grid%ihi_node, &
              jmin_node=>grid%jlo_node, jmax_node=>grid%jhi_node, &
              imin_cell=>grid%ilo_cell, imax_cell=>grid%ihi_cell, &
              jmin_cell=>grid%jlo_cell, jmax_cell=>grid%jhi_cell)

      allocate(evolved_corner_rho(imin_node:imax_node, jmin_node:jmax_node))
      allocate(evolved_corner_u(imin_node:imax_node, jmin_node:jmax_node))
      allocate(evolved_corner_v(imin_node:imax_node, jmin_node:jmax_node))
      allocate(evolved_corner_p(imin_node:imax_node, jmin_node:jmax_node))

      allocate(evolved_lr_mid_rho(imin_cell:imax_cell, jmin_node:jmax_node))
      allocate(evolved_lr_mid_u(imin_cell:imax_cell, jmin_node:jmax_node))
      allocate(evolved_lr_mid_v(imin_cell:imax_cell, jmin_node:jmax_node))
      allocate(evolved_lr_mid_p(imin_cell:imax_cell, jmin_node:jmax_node))

      allocate(evolved_du_mid_rho(imin_node:imax_node, jmin_cell:jmax_cell))
      allocate(evolved_du_mid_u(imin_node:imax_node, jmin_cell:jmax_cell))
      allocate(evolved_du_mid_v(imin_node:imax_node, jmin_cell:jmax_cell))
      allocate(evolved_du_mid_p(imin_node:imax_node, jmin_cell:jmax_cell))

    end associate

    write(stage_name, '(2(a, i0))') 'iter_', self%iteration, 'stage_', stage

    call reconstructor%set_cell_average_pointers(rho=rho, p=p, lbounds=lbounds)
    call self%apply_primitive_bc(rho=rho, u=u, v=v, p=p, lbounds=lbounds, &
                                 bc_plus_x=bc_plus_x, bc_minus_x=bc_minus_x, &
                                 bc_plus_y=bc_plus_y, bc_minus_y=bc_minus_y)

    ! Now we can reconstruct the entire domain
    call debug_print('Reconstructing density', __FILE__, __LINE__)
    call reconstructor%reconstruct(primitive_var=rho, lbounds=lbounds, &
                                   reconstructed_var=rho_recon_state, name='rho', stage_name=stage_name)

    call debug_print('Reconstructing x-velocity', __FILE__, __LINE__)
    call reconstructor%reconstruct(primitive_var=u, lbounds=lbounds, &
                                   reconstructed_var=u_recon_state, name='u', stage_name=stage_name)

    call debug_print('Reconstructing y-velocity', __FILE__, __LINE__)
    call reconstructor%reconstruct(primitive_var=v, lbounds=lbounds, &
                                   reconstructed_var=v_recon_state, name='v', stage_name=stage_name)

    call debug_print('Reconstructing pressure', __FILE__, __LINE__)
    call reconstructor%reconstruct(primitive_var=p, lbounds=lbounds, &
                                   reconstructed_var=p_recon_state, name='p', stage_name=stage_name)

    ! The gradients have to be applied to the boundaries as well, since reconstructing
    ! at P'(x,y) requires the cell gradient
    call self%apply_gradient_bc(reconstructor, &
                                bc_plus_x=bc_plus_x, bc_minus_x=bc_minus_x, &
                                bc_plus_y=bc_plus_y, bc_minus_y=bc_minus_y)

    ! Apply the reconstructed state to the ghost layers
    call self%apply_reconstructed_bc(recon_rho=rho_recon_state, recon_u=u_recon_state, &
                                     recon_v=v_recon_state, recon_p=p_recon_state, lbounds=lbounds, &
                                     bc_plus_x=bc_plus_x, bc_minus_x=bc_minus_x, &
                                     bc_plus_y=bc_plus_y, bc_minus_y=bc_minus_y)

    E0%reconstructed_rho => rho_recon_state
    E0%reconstructed_u => u_recon_state
    E0%reconstructed_v => v_recon_state
    E0%reconstructed_p => p_recon_state
    ! call E0%set_reconstructed_state_pointers(rho=rho_recon_state, &
    !                                                               u=u_recon_state, &
    !                                                               v=v_recon_state, &
    !                                                               p=p_recon_state, &
    !                                                               lbounds=lbounds)

    ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of down/up edge vectors
    call debug_print('Evolving down/up midpoints', __FILE__, __LINE__)
    evo_bounds = lbound(evolved_du_mid_rho)
    call E0%evolve(dt=dt, evolved_rho=evolved_du_mid_rho, evolved_u=evolved_du_mid_u, &
                   evolved_v=evolved_du_mid_v, evolved_p=evolved_du_mid_p, &
                   location='down/up midpoint', &
                   lbounds=evo_bounds, &
                   error_code=error_code)

    ! Evolve, i.e. E0(R_omega), at all midpoint nodes that are composed of left/right edge vectors
    call debug_print('Evolving left/right midpoints', __FILE__, __LINE__)
    evo_bounds = lbound(evolved_lr_mid_rho)
    call E0%evolve(dt=dt, evolved_rho=evolved_lr_mid_rho, evolved_u=evolved_lr_mid_u, &
                   evolved_v=evolved_lr_mid_v, evolved_p=evolved_lr_mid_p, &
                   location='left/right midpoint', &
                   lbounds=evo_bounds, &
                   error_code=error_code)

    ! Evolve, i.e. E0(R_omega), at all corner nodes
    call debug_print('Evolving corner nodes', __FILE__, __LINE__)
    evo_bounds = lbound(evolved_corner_rho)
    call E0%evolve(dt=dt, evolved_rho=evolved_corner_rho, evolved_u=evolved_corner_u, &
                   evolved_v=evolved_corner_v, evolved_p=evolved_corner_p, &
                   location='corner', &
                   lbounds=evo_bounds, &
                   error_code=error_code)

    call self%flux_edges(grid=grid, &
                         evolved_corner_rho=evolved_corner_rho, &
                         evolved_corner_u=evolved_corner_u, &
                         evolved_corner_v=evolved_corner_v, &
                         evolved_corner_p=evolved_corner_p, &
                         evolved_lr_mid_rho=evolved_lr_mid_rho, &
                         evolved_lr_mid_u=evolved_lr_mid_u, &
                         evolved_lr_mid_v=evolved_lr_mid_v, &
                         evolved_lr_mid_p=evolved_lr_mid_p, &
                         evolved_du_mid_rho=evolved_du_mid_rho, &
                         evolved_du_mid_u=evolved_du_mid_u, &
                         evolved_du_mid_v=evolved_du_mid_v, &
                         evolved_du_mid_p=evolved_du_mid_p, &
                         d_rho_dt=d_rho_dt, &
                         d_rhou_dt=d_rho_u_dt, &
                         d_rhov_dt=d_rho_v_dt, &
                         d_rhoE_dt=d_rho_E_dt)

    ! Now deallocate everything

    deallocate(bc_plus_x)
    deallocate(bc_plus_y)
    deallocate(bc_minus_x)
    deallocate(bc_minus_y)

    deallocate(reconstructor)

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

  subroutine apply_primitive_bc(self, lbounds, rho, u, v, p, &
                                bc_plus_x, bc_minus_x, bc_plus_y, bc_minus_y)
    class(fvleg_solver_t), intent(inout) :: self
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

  subroutine apply_reconstructed_bc(self, lbounds, recon_rho, recon_u, recon_v, recon_p, &
                                    bc_plus_x, bc_minus_x, bc_plus_y, bc_minus_y)
    class(fvleg_solver_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_rho
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_u
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_v
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(inout) :: recon_p
    class(boundary_condition_t), intent(inout):: bc_plus_x
    class(boundary_condition_t), intent(inout):: bc_plus_y
    class(boundary_condition_t), intent(inout):: bc_minus_x
    class(boundary_condition_t), intent(inout):: bc_minus_y

    integer(ik) :: priority
    integer(ik) :: max_priority_bc !< highest goes first

    call debug_print('Running fvleg_solver_t%apply_reconstructed_state_bc()', __FILE__, __LINE__)

    max_priority_bc = max(bc_plus_x%priority, bc_plus_y%priority, &
                          bc_minus_x%priority, bc_minus_y%priority)

    do priority = max_priority_bc, 0, -1

      if(bc_plus_x%priority == priority) then
        call bc_plus_x%apply_reconstructed_state_bc(recon_rho=recon_rho, &
                                                    recon_u=recon_u, &
                                                    recon_v=recon_v, &
                                                    recon_p=recon_p, lbounds=lbounds)
      end if

      if(bc_plus_y%priority == priority) then
        call bc_plus_y%apply_reconstructed_state_bc(recon_rho=recon_rho, &
                                                    recon_u=recon_u, &
                                                    recon_v=recon_v, &
                                                    recon_p=recon_p, lbounds=lbounds)
      end if

      if(bc_minus_x%priority == priority) then
        call bc_minus_x%apply_reconstructed_state_bc(recon_rho=recon_rho, &
                                                     recon_u=recon_u, &
                                                     recon_v=recon_v, &
                                                     recon_p=recon_p, lbounds=lbounds)
      end if

      if(bc_minus_y%priority == priority) then
        call bc_minus_y%apply_reconstructed_state_bc(recon_rho=recon_rho, &
                                                     recon_u=recon_u, &
                                                     recon_v=recon_v, &
                                                     recon_p=recon_p, lbounds=lbounds)
      end if

    end do

  end subroutine apply_reconstructed_bc

  subroutine apply_gradient_bc(self, reconstructor, bc_plus_x, bc_minus_x, bc_plus_y, bc_minus_y)

    class(fvleg_solver_t), intent(inout) :: self
    class(abstract_reconstruction_t), intent(inout) :: reconstructor
    class(boundary_condition_t), intent(inout):: bc_plus_x
    class(boundary_condition_t), intent(inout):: bc_plus_y
    class(boundary_condition_t), intent(inout):: bc_minus_x
    class(boundary_condition_t), intent(inout):: bc_minus_y

    integer(ik) :: priority
    integer(ik) :: max_priority_bc !< highest goes first
    integer(ik), dimension(2) :: lbounds

    call debug_print('Running fvleg_solver_t%apply_gradient_bc()', __FILE__, __LINE__)

    lbounds = lbound(reconstructor%grad_x_rho)
    max_priority_bc = max(bc_plus_x%priority, bc_plus_y%priority, &
                          bc_minus_x%priority, bc_minus_y%priority)

    do priority = max_priority_bc, 0, -1

      if(bc_plus_x%priority == priority) then
        call bc_plus_x%apply_gradient_bc(grad_x=reconstructor%grad_x_rho, &
                                         grad_y=reconstructor%grad_y_rho, &
                                         lbounds=lbounds)
        call bc_plus_x%apply_gradient_bc(grad_x=reconstructor%grad_x_p, &
                                         grad_y=reconstructor%grad_y_p, &
                                         lbounds=lbounds)
      end if

      if(bc_plus_y%priority == priority) then
        call bc_plus_y%apply_gradient_bc(grad_x=reconstructor%grad_x_rho, &
                                         grad_y=reconstructor%grad_y_rho, &
                                         lbounds=lbounds)
        call bc_plus_y%apply_gradient_bc(grad_x=reconstructor%grad_x_p, &
                                         grad_y=reconstructor%grad_y_p, &
                                         lbounds=lbounds)
      end if

      if(bc_minus_x%priority == priority) then
        call bc_minus_x%apply_gradient_bc(grad_x=reconstructor%grad_x_rho, &
                                          grad_y=reconstructor%grad_y_rho, &
                                          lbounds=lbounds)
        call bc_minus_x%apply_gradient_bc(grad_x=reconstructor%grad_x_p, &
                                          grad_y=reconstructor%grad_y_p, &
                                          lbounds=lbounds)
      end if

      if(bc_minus_y%priority == priority) then
        call bc_minus_y%apply_gradient_bc(grad_x=reconstructor%grad_x_rho, &
                                          grad_y=reconstructor%grad_y_rho, &
                                          lbounds=lbounds)
        call bc_minus_y%apply_gradient_bc(grad_x=reconstructor%grad_x_p, &
                                          grad_y=reconstructor%grad_y_p, &
                                          lbounds=lbounds)
      end if

    end do

  end subroutine apply_gradient_bc

  subroutine flux_edges(grid, &
                        evolved_corner_rho, evolved_corner_u, evolved_corner_v, evolved_corner_p, &
                        evolved_lr_mid_rho, evolved_lr_mid_u, evolved_lr_mid_v, evolved_lr_mid_p, &
                        evolved_du_mid_rho, evolved_du_mid_u, evolved_du_mid_v, evolved_du_mid_p, &
                        d_rho_dt, d_rhou_dt, d_rhov_dt, d_rhoE_dt)
    !< Evaluate the fluxes along the edges. This is equation 13 in the paper
    class(grid_t), intent(in) :: grid

    real(rk), dimension(grid%ilo_node:, grid%jlo_node:), contiguous, intent(in) :: evolved_corner_rho !< (i,j); Reconstructed rho at the corners
    real(rk), dimension(grid%ilo_node:, grid%jlo_node:), contiguous, intent(in) :: evolved_corner_u   !< (i,j); Reconstructed u at the corners
    real(rk), dimension(grid%ilo_node:, grid%jlo_node:), contiguous, intent(in) :: evolved_corner_v   !< (i,j); Reconstructed v at the corners
    real(rk), dimension(grid%ilo_node:, grid%jlo_node:), contiguous, intent(in) :: evolved_corner_p   !< (i,j); Reconstructed p at the corners
    real(rk), dimension(grid%ilo_cell:, grid%jlo_node:), contiguous, intent(in) :: evolved_lr_mid_rho !< (i,j); Reconstructed rho at the left/right midpoints
    real(rk), dimension(grid%ilo_cell:, grid%jlo_node:), contiguous, intent(in) :: evolved_lr_mid_u   !< (i,j); Reconstructed u at the left/right midpoints
    real(rk), dimension(grid%ilo_cell:, grid%jlo_node:), contiguous, intent(in) :: evolved_lr_mid_v   !< (i,j); Reconstructed v at the left/right midpoints
    real(rk), dimension(grid%ilo_cell:, grid%jlo_node:), contiguous, intent(in) :: evolved_lr_mid_p   !< (i,j); Reconstructed p at the left/right midpoints
    real(rk), dimension(grid%ilo_node:, grid%jlo_cell:), contiguous, intent(in) :: evolved_du_mid_rho !< (i,j); Reconstructed rho at the down/up midpoints
    real(rk), dimension(grid%ilo_node:, grid%jlo_cell:), contiguous, intent(in) :: evolved_du_mid_u   !< (i,j); Reconstructed u at the down/up midpoints
    real(rk), dimension(grid%ilo_node:, grid%jlo_cell:), contiguous, intent(in) :: evolved_du_mid_v   !< (i,j); Reconstructed v at the down/up midpoints
    real(rk), dimension(grid%ilo_node:, grid%jlo_cell:), contiguous, intent(in) :: evolved_du_mid_p   !< (i,j); Reconstructed p at the down/up midpoints
    real(rk), dimension(grid%ilo_bc_cell:, grid%jlo_bc_cell:), contiguous, intent(inout) :: d_rho_dt
    real(rk), dimension(grid%ilo_bc_cell:, grid%jlo_bc_cell:), contiguous, intent(inout) :: d_rhou_dt
    real(rk), dimension(grid%ilo_bc_cell:, grid%jlo_bc_cell:), contiguous, intent(inout) :: d_rhov_dt
    real(rk), dimension(grid%ilo_bc_cell:, grid%jlo_bc_cell:), contiguous, intent(inout) :: d_rhoE_dt

    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: i, j, k, edge, xy
    real(rk), dimension(2, 4) :: n_hat
    real(rk), dimension(4) :: delta_l
    integer(ik), dimension(2) :: bounds

    type(flux_array_t) :: corner_fluxes
    type(flux_array_t) :: downup_mid_fluxes
    type(flux_array_t) :: leftright_mid_fluxes

    real(rk) :: f_sum, g_sum

    real(rk), dimension(4) :: bottom_flux
    real(rk), dimension(4) :: right_flux
    real(rk), dimension(4) :: top_flux
    real(rk), dimension(4) :: left_flux

    real(rk) :: rho_flux, ave_rho_flux
    real(rk) :: rhou_flux, ave_rhou_flux
    real(rk) :: rhov_flux, ave_rhov_flux
    real(rk) :: rhoE_flux, ave_rhoE_flux
    real(rk) :: threshold

    real(rk), parameter :: REL_THRESHOLD = 1e-5_rk

    call debug_print('Running fvleg_solver_t%flux_edges()', __FILE__, __LINE__)

    ! Get the flux arrays for each corner node or midpoint
    ilo = grid%ilo_node; ihi = grid%ihi_node
    jlo = grid%jlo_node; jhi = grid%jhi_node
    bounds = [ilo, jlo]
    corner_fluxes = get_fluxes(rho=evolved_corner_rho, u=evolved_corner_u, v=evolved_corner_v, &
                               p=evolved_corner_p, lbounds=bounds)

    ilo = grid%ilo_cell; ihi = grid%ihi_cell
    jlo = grid%jlo_node; jhi = grid%jhi_node
    bounds = [ilo, jlo]
    leftright_mid_fluxes = get_fluxes(rho=evolved_lr_mid_rho, u=evolved_lr_mid_u, v=evolved_lr_mid_v, &
                                      p=evolved_lr_mid_p, lbounds=bounds)

    ilo = grid%ilo_node; ihi = grid%ihi_node
    jlo = grid%jlo_cell; jhi = grid%jhi_cell
    bounds = [ilo, jlo]
    downup_mid_fluxes = get_fluxes(rho=evolved_du_mid_rho, u=evolved_du_mid_u, v=evolved_du_mid_v, &
                                   p=evolved_du_mid_p, lbounds=bounds)

    ilo = grid%ilo_cell; ihi = grid%ihi_cell
    jlo = grid%jlo_cell; jhi = grid%jhi_cell
    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j, k, delta_l, n_hat) &
    !$omp private(bottom_flux, right_flux, left_flux, top_flux) &
    !$omp private(f_sum, g_sum) &
    !$omp shared(grid, corner_fluxes, leftright_mid_fluxes, downup_mid_fluxes) &
    !$omp shared(d_rho_dt, d_rhou_dt, d_rhov_dt, d_rhoE_dt) &
    !$omp private(rho_flux, ave_rho_flux, rhou_flux, ave_rhou_flux, rhov_flux, ave_rhov_flux, rhoE_flux, ave_rhoE_flux, threshold)
    !$omp do simd
    do j = jlo, jhi
      do i = ilo, ihi

        delta_l = grid%cell_edge_lengths(:, i, j)

        do edge = 1, 4
          do xy = 1, 2
            n_hat(xy, edge) = grid%cell_edge_norm_vectors(xy, edge, i, j)
          end do
        end do

        associate(F_c1=>corner_fluxes%F(:, i, j), &
                  G_c1=>corner_fluxes%G(:, i, j), &
                  F_c2=>corner_fluxes%F(:, i + 1, j), &
                  G_c2=>corner_fluxes%G(:, i + 1, j), &
                  F_c3=>corner_fluxes%F(:, i + 1, j + 1), &
                  G_c3=>corner_fluxes%G(:, i + 1, j + 1), &
                  F_c4=>corner_fluxes%F(:, i, j + 1), &
                  G_c4=>corner_fluxes%G(:, i, j + 1), &
                  F_m1=>leftright_mid_fluxes%F(:, i, j), &
                  G_m1=>leftright_mid_fluxes%G(:, i, j), &
                  F_m2=>downup_mid_fluxes%F(:, i + 1, j), &
                  G_m2=>downup_mid_fluxes%G(:, i + 1, j), &
                  F_m3=>leftright_mid_fluxes%F(:, i, j + 1), &
                  G_m3=>leftright_mid_fluxes%G(:, i, j + 1), &
                  G_m4=>downup_mid_fluxes%G(:, i, j), &
                  F_m4=>downup_mid_fluxes%F(:, i, j), &
                  n_hat_1=>grid%cell_edge_norm_vectors(:, 1, i, j), &
                  n_hat_2=>grid%cell_edge_norm_vectors(:, 2, i, j), &
                  n_hat_3=>grid%cell_edge_norm_vectors(:, 3, i, j), &
                  n_hat_4=>grid%cell_edge_norm_vectors(:, 4, i, j))

          ! Bottom
          ! bottom_flux (rho, rhou, rhov, rhoE)
          do k = 1, 4
            f_sum = sum([F_c1(k), 4.0_rk * F_m1(k), F_c2(k)]) * n_hat_1(1)
            g_sum = sum([G_c1(k), 4.0_rk * G_m1(k), G_c2(k)]) * n_hat_1(2)
            bottom_flux(k) = (f_sum + g_sum) * (delta_l(1) / 6.0_rk)
          end do

          ! Right
          ! right_flux (rho, rhou, rhov, rhoE)
          do k = 1, 4
            f_sum = sum([F_c2(k), 4.0_rk * F_m2(k), F_c3(k)]) * n_hat_2(1)
            g_sum = sum([G_c2(k), 4.0_rk * G_m2(k), G_c3(k)]) * n_hat_2(2)
            right_flux(k) = (f_sum + g_sum) * (delta_l(2) / 6.0_rk)
          end do

          ! Top
          ! top_flux (rho, rhou, rhov, rhoE)
          do k = 1, 4
            f_sum = sum([F_c3(k), 4.0_rk * F_m3(k), F_c4(k)]) * n_hat_3(1)
            g_sum = sum([G_c3(k), 4.0_rk * G_m3(k), G_c4(k)]) * n_hat_3(2)
            top_flux(k) = (f_sum + g_sum) * (delta_l(3) / 6.0_rk)
          end do

          ! Left
          ! left_flux (rho, rhou, rhov, rhoE)
          do k = 1, 4
            f_sum = sum([F_c4(k), 4.0_rk * F_m4(k), F_c1(k)]) * n_hat_4(1)
            g_sum = sum([G_c4(k), 4.0_rk * G_m4(k), G_c1(k)]) * n_hat_4(2)
            left_flux(k) = (f_sum + g_sum) * (delta_l(4) / 6.0_rk)
          end do

          ave_rho_flux = sum([abs(left_flux(1)), abs(right_flux(1)), abs(top_flux(1)), abs(bottom_flux(1))]) / 4.0_rk
          ave_rhou_flux = sum([abs(left_flux(2)), abs(right_flux(2)), abs(top_flux(2)), abs(bottom_flux(2))]) / 4.0_rk
          ave_rhov_flux = sum([abs(left_flux(3)), abs(right_flux(3)), abs(top_flux(3)), abs(bottom_flux(3))]) / 4.0_rk
          ave_rhoE_flux = sum([abs(left_flux(4)), abs(right_flux(4)), abs(top_flux(4)), abs(bottom_flux(4))]) / 4.0_rk

          ! Now sum them all up together
          rho_flux = neumaier_sum_4([left_flux(1), right_flux(1), top_flux(1), bottom_flux(1)])
          rhou_flux = neumaier_sum_4([left_flux(2), right_flux(2), top_flux(2), bottom_flux(2)])
          rhov_flux = neumaier_sum_4([left_flux(3), right_flux(3), top_flux(3), bottom_flux(3)])
          rhoE_flux = neumaier_sum_4([left_flux(4), right_flux(4), top_flux(4), bottom_flux(4)])

          threshold = abs(ave_rho_flux) * REL_THRESHOLD
          if(abs(rho_flux) < threshold .or. abs(rho_flux) < epsilon(1.0_rk)) then
            rho_flux = 0.0_rk
          end if

          threshold = abs(ave_rhou_flux) * REL_THRESHOLD
          if(abs(rhou_flux) < threshold .or. abs(rhou_flux) < epsilon(1.0_rk)) then
            rhou_flux = 0.0_rk
          end if

          threshold = abs(ave_rhov_flux) * REL_THRESHOLD
          if(abs(rhov_flux) < threshold .or. abs(rhov_flux) < epsilon(1.0_rk)) then
            rhov_flux = 0.0_rk
          end if

          threshold = abs(ave_rhoE_flux) * REL_THRESHOLD
          if(abs(rhoE_flux) < threshold .or. abs(rhoE_flux) < epsilon(1.0_rk)) then
            rhoE_flux = 0.0_rk
          end if

          d_rho_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * rho_flux
          d_rhou_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * rhou_flux
          d_rhov_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * rhov_flux
          d_rhoE_dt(i, j) = (-1.0_rk / grid%cell_volume(i, j)) * rhoE_flux
        end associate

      end do ! i
    end do ! j
    !$omp end do simd
    !$omp end parallel

    ilo = grid%ilo_bc_cell; ihi = grid%ihi_bc_cell
    jlo = grid%jlo_bc_cell; jhi = grid%jhi_bc_cell
    d_rho_dt(ilo, :) = 0.0_rk
    d_rho_dt(ihi, :) = 0.0_rk
    d_rho_dt(:, jlo) = 0.0_rk
    d_rho_dt(:, jhi) = 0.0_rk

    d_rhou_dt(ilo, :) = 0.0_rk
    d_rhou_dt(ihi, :) = 0.0_rk
    d_rhou_dt(:, jlo) = 0.0_rk
    d_rhou_dt(:, jhi) = 0.0_rk

    d_rhov_dt(ilo, :) = 0.0_rk
    d_rhov_dt(ihi, :) = 0.0_rk
    d_rhov_dt(:, jlo) = 0.0_rk
    d_rhov_dt(:, jhi) = 0.0_rk

    d_rhoE_dt(ilo, :) = 0.0_rk
    d_rhoE_dt(ihi, :) = 0.0_rk
    d_rhoE_dt(:, jlo) = 0.0_rk
    d_rhoE_dt(:, jhi) = 0.0_rk
  end subroutine flux_edges

  subroutine finalize(self)
    !< Class finalizer
    type(fvleg_solver_t), intent(inout) :: self

    call debug_print('Running fvleg_solver_t%finalize()', __FILE__, __LINE__)
  end subroutine finalize
end module mod_fvleg_solver
