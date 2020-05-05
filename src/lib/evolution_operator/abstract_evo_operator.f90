module mod_abstract_evo_operator
  !< Provide a base type for the evolution operators

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_input, only: input_t
  use mod_grid, only: grid_t
  use mod_mach_cone_collection, only: mach_cone_collection_t

  implicit none

  private
  public :: abstract_evo_operator_t

  type, abstract :: abstract_evo_operator_t
    !< This class provides the framework for any evolution operators. For example, the local evolution
    !< operator will inherit from this class. Most of the attributes are the same, just that the
    !< implementation varies between those who inherit this class

    class(grid_t), pointer :: grid => null()
    !< pointer to the grid object

    real(rk), dimension(:, :, :), pointer :: reconstructed_rho => null()
    !< (point 1:8, i, j); pointer to reconstructed density
    real(rk), dimension(:, :, :), pointer :: reconstructed_u => null()
    !< (point 1:8, i, j); pointer to reconstructed x-velocity
    real(rk), dimension(:, :, :), pointer :: reconstructed_v => null()
    !< (point 1:8, i, j); pointer to reconstructed y-velocity
    real(rk), dimension(:, :, :), pointer :: reconstructed_p => null()
    !< (point 1:8, i, j); pointer to reconstructed pressure

    class(abstract_reconstruction_t), pointer :: reconstruction_operator => null()
    !< pointer to the R_Omega operator used to provide values at the P' location

    integer(ik), dimension(:, :, :, :), allocatable :: leftright_midpoint_neighbors
    !< ((i,j), (cell 1:2), i, j); (i,j) indices of the cells on either side of the left/right vector midpoint

    integer(ik), dimension(:, :, :, :), allocatable :: downup_midpoint_neighbors
    !< ((i,j), (cell 1:2), i, j); (i,j) indices of the cells on either side of the down/up vector midpoint

    integer(ik), dimension(:, :, :, :), allocatable :: corner_neighbors
    !< ((i,j), (cell 1:4), i, j); (i,j) indices of the cells surrounding the corner vector set

    type(mach_cone_collection_t) :: leftright_midpoint_mach_cones
    type(mach_cone_collection_t) :: downup_midpoint_mach_cones
    type(mach_cone_collection_t) :: corner_mach_cones

    character(:), allocatable :: name  !< Name of the evolution operator
    real(rk) :: tau !< time increment to evolve

    integer(ik) :: iteration
    real(rk) :: time
    real(rk) :: time_step
  contains
    procedure(initialize), public, deferred :: initialize
    procedure(evolve_location), public, deferred :: evolve
    procedure, public, non_overridable :: set_grid_pointer
    procedure, public, non_overridable :: set_reconstructed_state_pointers
    procedure, public, non_overridable :: set_reconstruction_operator_pointer
    procedure, public, non_overridable :: set_tau
    procedure, public, non_overridable :: nullify_pointer_members
    procedure(copy_evo), public, deferred :: copy
    generic :: assignment(=) => copy
  end type abstract_evo_operator_t

  abstract interface
    subroutine initialize(self, input, grid_target, recon_operator_target)
      import :: abstract_evo_operator_t
      import :: input_t
      import :: grid_t
      import :: abstract_reconstruction_t
      import :: ik, rk
      class(abstract_evo_operator_t), intent(inout) :: self
      class(input_t), intent(in) :: input
      class(grid_t), intent(in), target :: grid_target
      class(abstract_reconstruction_t), intent(in), target :: recon_operator_target
    end subroutine

    subroutine copy_evo(out_evo, in_evo)
      import :: abstract_evo_operator_t
      class(abstract_evo_operator_t), intent(in) :: in_evo
      class(abstract_evo_operator_t), intent(inout) :: out_evo
    end subroutine

    subroutine evolve_location(self, location, evolved_rho, evolved_u, evolved_v, evolved_p, lbounds, error_code)
      import :: abstract_evo_operator_t
      import :: rk, ik
      class(abstract_evo_operator_t), intent(inout) :: self

      integer(ik), dimension(2), intent(in) :: lbounds
      character(len=*), intent(in) :: location !< Mach cone location ['corner', 'left/right midpoint', or 'down/up midpoint']
      integer(ik), intent(out) :: error_code
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: evolved_rho
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: evolved_u
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: evolved_v
      real(rk), dimension(lbounds(1):, lbounds(2):), intent(out) :: evolved_p
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

  subroutine set_reconstructed_state_pointers(self, rho, u, v, p, lbounds)
    class(abstract_evo_operator_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(in), target :: rho
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(in), target :: u
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(in), target :: v
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(in), target :: p
    self%reconstructed_rho => rho
    self%reconstructed_u => u
    self%reconstructed_v => v
    self%reconstructed_p => p
  end subroutine

  subroutine set_reconstruction_operator_pointer(self, operator_target)
    class(abstract_evo_operator_t), intent(inout) :: self
    class(abstract_reconstruction_t), intent(in), target :: operator_target
    self%reconstruction_operator => operator_target
  end subroutine

  subroutine nullify_pointer_members(self)
    class(abstract_evo_operator_t), intent(inout) :: self

    if(associated(self%grid)) nullify(self%grid)

    if(associated(self%reconstructed_rho)) nullify(self%reconstructed_rho)
    if(associated(self%reconstructed_u)) nullify(self%reconstructed_u)
    if(associated(self%reconstructed_v)) nullify(self%reconstructed_v)
    if(associated(self%reconstructed_p)) nullify(self%reconstructed_p)
    if(associated(self%reconstruction_operator)) nullify(self%reconstruction_operator)

  end subroutine nullify_pointer_members

end module mod_abstract_evo_operator
