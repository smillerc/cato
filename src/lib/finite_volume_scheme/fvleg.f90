module mod_fvleg

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_surrogate, only: surrogate
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t
  use mod_input, only: input_t
  use mod_finite_volume_schemes, only: finite_volume_scheme_t
  use mod_reconstruction_factory, only: reconstruction_factory_t
  use mod_local_evo_operator, only: local_evo_operator_t
  use mod_second_order_reconstruction, only: second_order_reconstruction_t
  use mod_flux_tensor, only: H => flux_tensor_t
  ! use mod_first_order_reconstruction, only: first_order_reconstruction_t

  implicit none
  private
  public :: fvleg_t

  type, extends(finite_volume_scheme_t) :: fvleg_t
    !< Implementation of the finite volume local evolution Galerkin (fvleg) scheme type

  contains
    procedure, public :: reconstruct => reconstruct_fvleg
    ! procedure, public :: eval_fluxes => eval_fvleg_fluxes
    procedure, public :: t => time_derivative

    procedure, public :: add => add_fvleg
    procedure, public :: multiply => multiply_fvleg
    procedure, public :: assign => assign_fvleg
    final :: finalize
  end type

  interface fvleg_t
    module procedure :: constructor
  end interface

contains

  function constructor(input) result(scheme)
    !< Construct the finite volume local evolution Galerkin (fvleg) scheme
    class(input_t), intent(in) :: input
    class(fvleg_t), allocatable :: scheme
    type(reconstruction_factory_t) :: recon_factory
    ! type(grid_factory_t) :: grid_factory
    integer(ik) :: alloc_status

    allocate(fvleg_t :: scheme)

    ! grid_factory = grid_factory_t(input)
    ! scheme%grid = grid_factory%create_grid()

    ! recon_factory = reconstruction_factory_t(input)
    ! scheme%reconstruction_operator = recon_factory%create_reconstruction(input)

    allocate(scheme%conserved_vars(4, input%ni - 1, input%nj - 1), stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%conserved_vars"

    allocate(scheme%evolved_corner_state(4, input%ni, input%nj), stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%evolved_corner_state"
    !< ((rho,u,v,p), i, j); Reconstructed U at each corner

    allocate(scheme%corner_reference_state(4, input%ni, input%nj), stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%corner_reference_state"
    !< ((rho, u ,v, p), i, j); Reference state (tilde) at each corner

    allocate(scheme%evolved_downup_midpoints_state(4, input%ni, input%nj), stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%evolved_downup_midpoints_state"
    !< ((rho,u,v,p), i, j); Reconstructed U at each midpoint on the up/down edges (edges 2 and 4)

    allocate(scheme%evolved_leftright_midpoints_state(4, input%ni, input%nj), stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%evolved_leftright_midpoints_state"
    !< ((rho,u,v,p), i, j); Reconstructed U at each midpoint on the left/right edges (edges 1 and 3)

    allocate(scheme%downup_midpoints_reference_state(4, input%ni, input%nj), stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%downup_midpoints_reference_state"
    !< ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the up/down edges (edges 2 and 4)

    allocate(scheme%leftright_midpoints_reference_state(4, input%ni, input%nj), stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%leftright_midpoints_reference_state"
    !< ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the left/right edges (edges 1 and 3)

    allocate(scheme%reconstructed_state(4, 4, 2, input%ni - 1, input%nj - 1), stat=alloc_status)
    !< ((rho, u ,v, p), point, node/midpoint, i, j); this is a cell-based value, so imax=ni-1, etc
    if(alloc_status /= 0) error stop "Unable to allocate fvleg_t%reconstructed_state"

  end function

  subroutine finalize(self)
    !< Implementation of the class cleanup
    type(fvleg_t), intent(inout) :: self
    integer(ik) :: alloc_status

    if(allocated(self%grid)) deallocate(self%grid, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%grid"

    if(allocated(self%reconstruction_operator)) deallocate(self%reconstruction_operator, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%reconstruction_operator"

    if(allocated(self%conserved_vars)) deallocate(self%conserved_vars, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%conserved_vars"

    if(allocated(self%evolved_corner_state)) deallocate(self%evolved_corner_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%evolved_corner_state"

    if(allocated(self%corner_reference_state)) deallocate(self%corner_reference_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%corner_reference_state"

    if(allocated(self%evolved_downup_midpoints_state)) deallocate(self%evolved_downup_midpoints_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%evolved_downup_midpoints_state"

    if(allocated(self%evolved_leftright_midpoints_state)) deallocate(self%evolved_leftright_midpoints_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%evolved_leftright_midpoints_state"

    if(allocated(self%downup_midpoints_reference_state)) deallocate(self%downup_midpoints_reference_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%downup_midpoints_reference_state"

    if(allocated(self%leftright_midpoints_reference_state)) deallocate(self%leftright_midpoints_reference_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%leftright_midpoints_reference_state"

    if(allocated(self%reconstructed_state)) deallocate(self%reconstructed_state, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate fvleg_t%reconstructed_state"

  end subroutine

  subroutine reconstruct_fvleg(self)
    !< Implementation of the FVLEG reconstruction.
    !< This reconstructs the entire grid at all the nodes/midpoints
    class(fvleg_t), intent(inout) :: self

    call self%reconstruction_operator%reconstruct_domain(self%reconstructed_state)

    ! do j = self%grid%jlo, self%grid%jhi
    !   do i = self%grid%ilo, self%grid%ihi

    !     ! The cell average (and gradient if higher order) is resused for each cell
    !     call self%reconstruction_operator%select_and_find_gradient(i,j)

    !     ! First do corners, then to midpoints
    !     do n = 1, ubound(self%grid%reconstructed_state, dim=3)

    !       ! Loop through each point (N1-N4, and M1-M4)
    !       do p = 1, ubound(self%grid%reconstructed_state, dim=2)

    !         associate(V_bar => self%reconstructed_state, &
    !                   R_omega => self%reconstruction_operator%reconstruct, &
    !                   xy => self%grid%cell_node_xy(:, p, n, i, j))

    !           ! reconstructed_state(rho:p, point, node/midpoint, i, j)
    !           V_bar(:, p, n, i, j) = R_omega(xy)
    !         end associate

    !       end do
    !     end do
    !   end do
    ! end do

  end subroutine

  function time_derivative(self) result(dU_dt)
    !< Implementation of dU/dt = 1/Omega_ij Sum(1,k) Integral(H . n dl)
    class(fvleg_t), intent(in) :: self
    class(integrand_t), allocatable :: dU_dt
    type(fvleg_t), allocatable :: local_dU_dt
    integer(ik) :: alloc_status

    allocate(local_dU_dt, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate local_dU_dt"

    local_dU_dt%grid = self%grid
    local_dU_dt%evolution_operator = self%evolution_operator
    local_dU_dt%reconstruction_operator = self%reconstruction_operator

    ! // TODO: fix this
    ! allocate(local_dU_dt%reconstructed_state(size(self%reconstructed_state, dim=1), &
    !                                          size(self%reconstructed_state, dim=2), &
    !                                          size(self%reconstructed_state, dim=3), &
    !                                          size(self%reconstructed_state, dim=4), &
    !                                          size(self%reconstructed_state, dim=5)), &
    !          stat=alloc_status)
    ! if(alloc_status /= 0) error stop "Unable to allocate local_dU_dt%reconstructed_state"

    ! allocate(local_dU_dt%reference_state(size(self%reference_state, dim=1), &
    !                                      size(self%reference_state, dim=2), &
    !                                      size(self%reference_state, dim=3), &
    !                                      size(self%reference_state, dim=4), &
    !                                      size(self%reference_state, dim=5)), &
    !          stat=alloc_status)
    ! if(alloc_status /= 0) error stop "Unable to allocate local_dU_dt%reference_state"

    ! allocate(local_dU_dt%conserved_vars(size(self%conserved_vars, dim=1), &
    !                                     size(self%conserved_vars, dim=2), &
    !                                     size(self%conserved_vars, dim=3)), &
    !          stat=alloc_status)
    ! if(alloc_status /= 0) error stop "Unable to allocate local_dU_dt%conserved_vars"

    ! Set the RHS of Eq. 3 in the main text
    ! local_dU_dt%conserved_vars = self%integrate_fluxes() //TODO: fix

    call move_alloc(local_dU_dt, dU_dt)

  end function time_derivative

  pure function integrate_fluxes(self) result(rhs)
    !< Evaluate the fluxes along the edges. This is equation 13 in the paper
    use mod_flux_tensor, only: operator(.dot.)
    class(fvleg_t), intent(in) :: self
    real(rk), dimension(4) :: edge_flux
    real(rk), dimension(size(self%conserved_vars, dim=1), &
                        size(self%conserved_vars, dim=2), &
                        size(self%conserved_vars, dim=3)) :: rhs !< RHS of Eq 3 in the main text
    ! type(flux_tensor_t) :: H1, H2, H3
    integer(ik) :: i, j, k

    do j = self%grid%jlo, self%grid%jhi
      do i = self%grid%ilo, self%grid%ihi

        edge_flux = 0.0_rk

        do k = 1, 4  ! edge
          associate(E0=>self%evolution_operator, &
                    r_omega_u_bar=>self%reconstructed_state, &
                    n_hat=>self%grid%cell_edge_norm_vectors(:, k, i, j), &
                    delta_l=>self%grid%cell_edge_lengths(k, i, j))

            edge_flux = edge_flux + &
                        ( &
                        ((H(E0%evolve(r_omega_u_bar, i, j, edge=k, loc='k,1')) + &
                          H(E0%evolve(r_omega_u_bar, i, j, edge=k, loc='k,c')) * 4.0_rk + &
                          H(E0%evolve(r_omega_u_bar, i, j, edge=k, loc='k,2'))) .dot.n_hat) * (delta_l / 6.0_rk))

          end associate
        end do ! k

        rhs(:, i, j) = (-1.0_rk / self%grid%cell_volume(i, j)) * edge_flux

      end do ! i
    end do ! j

  end function integrate_fluxes

  function add_fvleg(lhs, rhs) result(operator_result)
    class(fvleg_t), intent(in) :: lhs
    class(integrand_t), intent(in) :: rhs
    class(integrand_t), allocatable :: operator_result
  end function add_fvleg

  function multiply_fvleg(lhs, rhs) result(operator_result)
    class(fvleg_t), intent(in) :: lhs
    class(integrand_t), allocatable :: operator_result
    real(rk), intent(in) :: rhs
  end function multiply_fvleg

  subroutine assign_fvleg(lhs, rhs)
    class(fvleg_t), intent(inout) :: lhs
    class(integrand_t), intent(in) :: rhs
  end subroutine assign_fvleg

end module mod_fvleg
