module mod_finite_volume_schemes

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_surrogate, only: surrogate
  use mod_boundary_conditions, only: boundary_condition_t
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t
  use mod_grid, only: grid_t
  use mod_abstract_evo_operator, only: abstract_evo_operator_t

  implicit none
  private
  public :: finite_volume_scheme_t

  type, abstract, extends(integrand_t) :: finite_volume_scheme_t
    !< Abstract representation of the finite volume scheme. This is essentially a puppeteer that
    !< manages the reconstruction, grid, and evolution operator to ultimately calculate the conserved
    !< state variables of each finite cell. The reconstruction, grid, and evolution implementations are passed
    !< on to decendents like the FVLEG scheme.

    class(abstract_reconstruction_t), allocatable :: reconstruction_operator
    !< R_Omega reconstruction operator used to reconstruct the corners/midpoints based on the cell
    !< average (and gradient if high(er) order reconstruction used)

    real(rk), dimension(:, :, :), allocatable :: conserved_vars
    !< ((rho, u ,v, p), i, j); Conserved variables for each cell

    class(abstract_evo_operator_t), allocatable :: evolution_operator
    !< Evolution operator to construct (rho, u, v, p) at each corner and midpoint

    class(grid_t), allocatable :: grid
    !< Grid class to hold geometric information (edge lengths, volumes, etc.)

    class(boundary_condition_t), allocatable :: bc_plus_x
    class(boundary_condition_t), allocatable :: bc_plus_y
    class(boundary_condition_t), allocatable :: bc_minus_x
    class(boundary_condition_t), allocatable :: bc_minus_y

    ! Corner/midpoint index convention
    ! --------------------------------
    !
    !   C----M----C----M----C
    !   |         |         |
    !   O    x    O    x    O
    !   |         |         |
    !   C----M----C----M----C
    !   |         |         |
    !   O    x    O    x    O
    !   |         |         |
    !   C----M----C----M----C

    ! C: corner, M: left/right midpoint, O: up/down midpoint, x: cell

    ! Since the reconstructed state at the corners (C) and midpoints (M) are reused by multiple cells,
    ! the datastructures are set up for maximum reuse.
    ! If they were indexed via cell, each cell would duplicate information since they share corners and midpoints

    real(rk), dimension(:, :, :), allocatable :: evolved_corner_state
    !< ((rho,u,v,p), i, j); Reconstructed U at each corner

    real(rk), dimension(:, :, :), allocatable :: corner_reference_state
    !< ((rho, u ,v, p), i, j); Reference state (tilde) at each corner

    ! Indexing the midpoints is a pain, so they're split by the up/down edges and left/right edges

    real(rk), dimension(:, :, :), allocatable :: evolved_downup_midpoints_state
    !< ((rho,u,v,p), i, j); Reconstructed U at each midpoint on the up/down edges (edges 2 and 4)

    real(rk), dimension(:, :, :), allocatable :: evolved_leftright_midpoints_state
    !< ((rho,u,v,p), i, j); Reconstructed U at each midpoint on the left/right edges (edges 1 and 3)

    real(rk), dimension(:, :, :), allocatable :: downup_midpoints_reference_state
    !< ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the down/up edges (edges 2 and 4)

    real(rk), dimension(:, :, :), allocatable :: leftright_midpoints_reference_state
    !< ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the left/right edges (edges 1 and 3)

    real(rk), dimension(:, :, :, :, :), allocatable :: reconstructed_state
    !< ((rho, u ,v, p), point, node/midpoint, i, j); The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints. Note, this DOES repeat nodes, since corners and midpoints are
    !< shared by neighboring cells, but each point has its own reconstructed value based on the parent cell's state

  contains
    procedure(reconstruct), public, deferred :: reconstruct
    ! procedure(eval_fluxes), public, deferred :: eval_fluxes
    procedure, public :: calculate_reference_state
  end type

  abstract interface
    subroutine reconstruct(self)
      import :: finite_volume_scheme_t
      class(finite_volume_scheme_t), intent(inout) :: self
    end subroutine

    pure function integrate_fluxes(self) result(rhs)
      import :: finite_volume_scheme_t
      import :: rk

      class(finite_volume_scheme_t), intent(in) :: self
      real(rk), dimension(4) :: edge_flux
      real(rk), dimension(size(self%conserved_vars, dim=1), &
                          size(self%conserved_vars, dim=2), &
                          size(self%conserved_vars, dim=3)) :: rhs !< RHS of Eq 3 in the main text
    end function

  end interface
contains

  subroutine calculate_reference_state(self)
    !< Calculate the reference state at each corner/midpoint. This is just an average of
    !< the neighboring cells
    class(finite_volume_scheme_t), intent(inout) :: self
    integer(ik) :: i, j

    ! left/right midpoints -> needs to average cells above and below
    do j = self%grid%jlo_cell, self%grid%jhi_cell
      do i = self%grid%ilo_cell, self%grid%ihi_cell
        associate(U_tilde=>self%leftright_midpoints_reference_state, U=>self%conserved_vars)
          U_tilde(:, i, j) = 0.5_rk * (U(:, i, j) + U(:, i, j - 1))
        end associate
      end do
    end do

    ! up/down midpoints -> needs to average cells right and left
    do j = self%grid%jlo_cell, self%grid%jhi_cell
      do i = self%grid%ilo_cell, self%grid%ihi_cell
        associate(U_tilde=>self%downup_midpoints_reference_state, U=>self%conserved_vars)
          U_tilde(:, i, j) = 0.5_rk * (U(:, i - 1, j) + U(:, i, j))
        end associate
      end do
    end do

    ! Corners
    do j = self%grid%jlo_cell, self%grid%jhi_cell
      do i = self%grid%ilo_cell, self%grid%ihi_cell
        associate(U_tilde=>self%corner_reference_state, U=>self%conserved_vars)
          U_tilde(:, i, j) = 0.25_rk * (U(:, i, j) + U(:, i - 1, j) + &
                                        U(:, i, j - 1) + U(:, i - 1, j - 1))
        end associate
      end do
    end do

  end subroutine

end module mod_finite_volume_schemes
