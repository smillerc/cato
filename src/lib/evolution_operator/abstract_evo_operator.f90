module mod_abstract_evo_operator
  !< Provide a base type for the evolution operators

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_mach_cone_geometry, only: mach_cone_geometry_t
  use mod_input, only: input_t
  use mod_grid, only: grid_t

  implicit none

  private
  public :: abstract_evo_operator_t

  type, abstract :: abstract_evo_operator_t
    !< This class provides the framework for any evolution operators. For example, the local evolution
    !< operator will inherit from this class. Most of the attributes are the same, just that the
    !< implementation varies between those who inherit this class

    class(grid_t), pointer :: grid => null()
    !< pointer to the grid object

    ! real(rk), dimension(:, :, :), pointer :: conserved_vars => null()
    !< pointer to U ((rho, u, v, p), i, j)

    real(rk), dimension(:, :, :, :, :), pointer :: reconstructed_state => null()
    !< ((rho, u ,v, p), point, node/midpoint, i, j);
    !< The reconstructed state of each point P with respect to its parent cell

    ! real(rk), dimension(:, :, :, :, :), pointer :: reference_state => null()
    ! !< pointer to the reference state at each point ((rho, u ,v, p), point, node/midpoint, i, j)

    class(abstract_reconstruction_t), pointer :: reconstruction_operator => null()
    !< pointer to the R_Omega operator used to provide values at the P' location

    character(:), allocatable :: name  !< Name of the evolution operator
    real(rk) :: tau !< time increment to evolve
  contains
    procedure(initialize), public, deferred :: initialize
    procedure(evolve_location), public, deferred :: evolve_leftright_midpoints
    procedure(evolve_location), public, deferred :: evolve_downup_midpoints
    procedure(evolve_location), public, deferred :: evolve_corners
    procedure, public, non_overridable :: set_grid_pointer
    procedure, public, non_overridable :: set_reconstructed_state_pointer
    procedure, public, non_overridable :: set_reconstruction_operator_pointer
    procedure, public, non_overridable :: set_tau
    procedure, public, non_overridable :: nullify_pointer_members
    procedure(copy_evo), public, deferred :: copy
    generic :: assignment(=) => copy
  end type abstract_evo_operator_t

  abstract interface
    subroutine initialize(self, input, grid_target, recon_operator_target, reconstructed_state_target, lbounds)
      import :: abstract_evo_operator_t
      import :: input_t
      import :: grid_t
      import :: abstract_reconstruction_t
      import :: ik, rk
      class(abstract_evo_operator_t), intent(inout) :: self
      class(input_t), intent(in) :: input
      class(grid_t), intent(in), target :: grid_target
      class(abstract_reconstruction_t), intent(in), target :: recon_operator_target
      integer(ik), dimension(5), intent(in) :: lbounds
      real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):, &
                          lbounds(4):, lbounds(5):), intent(in), target :: reconstructed_state_target

    end subroutine

    subroutine copy_evo(out_evo, in_evo)
      import :: abstract_evo_operator_t
      class(abstract_evo_operator_t), intent(in) :: in_evo
      class(abstract_evo_operator_t), intent(inout) :: out_evo
    end subroutine

    subroutine evolve_location(self, reference_state, evolved_state)
      import :: abstract_evo_operator_t
      import :: rk, ik
      class(abstract_evo_operator_t), intent(in) :: self
      ! !< ((rho, u ,v, p), point, node/midpoint, i, j); The reconstructed state of each point P with respect to its parent cell
      real(rk), dimension(:, :, :), intent(in) :: reference_state
      ! !< ((rho,u,v,p), i, j); Reference state (tilde) at each location
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

  subroutine set_grid_pointer(self, grid_target)
    class(abstract_evo_operator_t), intent(inout) :: self
    class(grid_t), intent(in), target :: grid_target
    self%grid => grid_target
  end subroutine

  subroutine set_reconstructed_state_pointer(self, reconstructed_state_target, lbounds)
    class(abstract_evo_operator_t), intent(inout) :: self
    integer(ik), dimension(5), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):, &
                        lbounds(4):, lbounds(5):), intent(in), target :: reconstructed_state_target
    self%reconstructed_state => reconstructed_state_target
  end subroutine

  subroutine set_reconstruction_operator_pointer(self, operator_target)
    class(abstract_evo_operator_t), intent(inout) :: self
    class(abstract_reconstruction_t), intent(in), target :: operator_target
    self%reconstruction_operator => operator_target
  end subroutine

  subroutine nullify_pointer_members(self)
    class(abstract_evo_operator_t), intent(inout) :: self
    if(associated(self%grid)) nullify(self%grid)
    if(associated(self%reconstructed_state)) nullify(self%reconstructed_state)
    if(associated(self%reconstruction_operator)) nullify(self%reconstruction_operator)
  end subroutine nullify_pointer_members

end module mod_abstract_evo_operator
