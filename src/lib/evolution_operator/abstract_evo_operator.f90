module mod_abstract_evo_operator
  !< Provide a base type for the evolution operators

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_mach_cone_geometry, only: mach_cone_geometry_t
  use mod_grid, only: grid_t

  implicit none

  private
  public :: abstract_evo_operator_t

  type, abstract :: abstract_evo_operator_t
    !< This class provides the framework for any evolution operators. For example, the local evolution
    !< operator will inherit from this class. Most of the attributes are the same, just that the
    !< implementation varies between those who inherit this class

    class(grid_t), pointer :: grid
    !< pointer to the grid object

    real(rk), dimension(:, :, :), pointer :: conserved_vars
    !< pointer to U ((rho, u, v, p), i, j)

    real(rk), dimension(:, :, :, :, :), pointer :: reference_state
    !< pointer to the reference state at each point ((rho, u ,v, p), point, node/midpoint, i, j)

    class(abstract_reconstruction_t), pointer :: reconstruction_operator
    !< pointer to the R_Omega operator used to provide values at the P' location

    character(:), allocatable :: name  !< Name of the evolution operator
    real(rk) :: tau !< time increment to evolve
  contains
    ! procedure(evolve), deferred, public :: evolve
    procedure(evolve_location), public, deferred :: evolve_leftright_midpoints
    procedure(evolve_location), public, deferred :: evolve_downup_midpoints
    procedure(evolve_location), public, deferred :: evolve_corners
    procedure, public, non_overridable :: set_tau
  end type abstract_evo_operator_t

  abstract interface
    ! pure function evolve(self, reconstructed_state, i, j, edge, loc) result(U_bar)
    !   !< Evolve the interface (edge) quantities based on their local reconstruction and neighboring cells
    !   import :: abstract_evo_operator_t
    !   import :: abstract_reconstruction_t
    !   import :: rk, ik

    !   real(rk), dimension(4) :: U_bar !< Conserved variable values specified location
    !   class(abstract_evo_operator_t), intent(in) :: self
    !   real(rk), dimension(:, :, :, :, :), intent(in) :: reconstructed_state
    !   !< ((rho, u ,v, p), point, node/midpoint, i, j);
    !   !< The node/midpoint dimension just selects which set of points,
    !   !< e.g. 1 - all corners, 2 - all midpoints

    !   integer(ik), intent(in) :: i, j !< Cell indices
    !   integer(ik), intent(in) :: edge !< Cell edge index
    !   character(len=3), intent(in) :: loc !< where is this evolution taking place? Needs to be one of the following ['k,1','k,c','k,2']
    ! end function evolve
    pure subroutine evolve_location(self, reconstructed_state, reference_state, evolved_state)
      import :: abstract_evo_operator_t
      import :: rk, ik
      class(abstract_evo_operator_t), intent(in) :: self
      real(rk), dimension(:, :, :, :, :), intent(in) :: reconstructed_state
      !< ((rho, u ,v, p), point, node/midpoint, i, j); The reconstructed state of each point P with respect to its parent cell
      real(rk), dimension(:, :, :), intent(in) :: reference_state
      !< ((rho,u,v,p), i, j); Reference state (tilde) at each location
      real(rk), dimension(:, :, :), intent(out) :: evolved_state
      !< ((rho,u,v,p), i, j); Reconstructed U at each location
    end subroutine
  end interface

contains
  subroutine set_tau(self, tau)
    !< Public interface to set the time increment
    class(abstract_evo_operator_t), intent(inout) :: self
    real(rk), intent(in) :: tau
    self%tau = abs(tau)
  end subroutine
end module mod_abstract_evo_operator
