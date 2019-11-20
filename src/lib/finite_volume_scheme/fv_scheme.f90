module mod_finite_volume_schemes

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_surrogate, only: surrogate
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t
  use mod_conserved_vars, only: conserved_vars_t
  use mod_grid, only: grid_t
  use mod_abstract_evo_operator, only: abstract_evo_operator_t

  implicit none
  private
  public :: finite_volume_scheme_t

  type, abstract, extends(integrand_t) :: finite_volume_scheme_t
    class(grid_t), allocatable :: grid
    class(abstract_reconstruction_t), allocatable :: reconstruction_operator

    real(rk), dimension(:, :, :), allocatable :: conserved_vars
    !< Conserved variables U = [rho, u ,v, p]; indexed via ([rho, u ,v, p], i, j)

    real(rk), dimension(:, :, :, :, :), allocatable :: reconstructed_state
    !< ([rho, u ,v, p], point, node/midpoint, i, j); The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints

    real(rk), dimension(:, :, :, :, :), allocatable :: reference_state
    !< ([rho, u ,v, p], point, node/midpoint, i, j)

    class(abstract_evo_operator_t), allocatable :: evolution_operator

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

    ! Midpoints
    do j = self%grid%jlo, self%grid%jhi
      do i = self%grid%ilo, self%grid%ihi
        associate(U_tilde=>self%reference_state, U=>self%conserved_vars)
          ! Midpoint 1 (bottom) ! = [rho,u,v,p]
          U_tilde(:, 1, 2, i, j) = 0.5_rk * (U(:, i, j) + U(:, i, j - 1))

          ! Midpoint 2 (right)
          U_tilde(:, 2, 2, i, j) = 0.5_rk * (U(:, i, j) + U(:, i + 1, j))

          ! Midpoint 3 (top)
          U_tilde(:, 3, 2, i, j) = 0.5_rk * (U(:, i, j) + U(:, i, j + 1))

          ! Midpoint 4 (left)
          U_tilde(:, 4, 2, i, j) = 0.5_rk * (U(:, i, j) + U(:, i - 1, j))
        end associate
      end do
    end do

    ! Corners
    do j = self%grid%jlo, self%grid%jhi
      do i = self%grid%ilo, self%grid%ihi
        associate(U_tilde=>self%reference_state, U=>self%conserved_vars)
          ! Corner 1 (bottom-left)
          U_tilde(:, 1, 1, i, j) = 0.25_rk * (U(:, i, j) + U(:, i, j - 1) + &
                                              U(:, i - 1, j) + U(:, i - 1, j - 1))

          ! Corner 2 (bottom-right)
          U_tilde(:, 1, 1, i, j) = 0.25_rk * (U(:, i, j) + U(:, i, j - 1) + &
                                              U(:, i + 1, j) + U(:, i + 1, j - 1))

          ! Corner 3 (top-right)
          U_tilde(:, 1, 1, i, j) = 0.25_rk * (U(:, i, j) + U(:, i, j + 1) + &
                                              U(:, i + 1, j) + U(:, i + 1, j + 1))

          ! Corner 4 (top-left)
          U_tilde(:, 1, 1, i, j) = 0.25_rk * (U(:, i, j) + U(:, i, j + 1) + &
                                              U(:, i - 1, j) + U(:, i - 1, j + 1))
        end associate
      end do
    end do

  end subroutine

end module mod_finite_volume_schemes
