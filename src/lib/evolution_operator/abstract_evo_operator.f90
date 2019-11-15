module mod_abstract_evo_operator
  !< Provide a base type for the evolution operators

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t

  implicit none

  private
  public :: abstract_evo_operator_t

  type, abstract :: abstract_evo_operator_t
    character(:), allocatable :: name
    real(rk), private :: tau !< time increment to evolve
  contains
    procedure(evolve), deferred, public :: evolve
    procedure, public, non_overridable :: set_tau
  end type abstract_evo_operator_t

  abstract interface
    function evolve(self, i, j, reconstruction_operator, edge, loc) result(U_bar)
      !< Evolve the interface (edge) quantities based on their local reconstruction and neighboring cells
      import :: abstract_evo_operator_t
      import :: abstract_reconstruction_t
      import :: rk, ik

      class(abstract_evo_operator_t), intent(in) :: self
      integer(ik), intent(in) :: i, j !< Cell indices
      class(abstract_reconstruction_t), pointer, intent(in) :: reconstruction_operator

      ! Optional arguments (depending on the order of the reconstruction)
      real(rk), dimension(4) :: U_bar !< Conserved variable values at the interface (edge)
      integer(ik), intent(in), optional :: edge !< Cell edge index
      character(len=3), intent(in), optional :: loc !< where is this evolution taking place? Needs to be one of the following ['k,1','k,c','k,2']
    end function evolve
  end interface

contains
  subroutine set_tau(self, tau)
    !< Public interface to set the time increment
    class(abstract_evo_operator_t), intent(inout) :: self
    real(rk), intent(in) :: tau
    self%tau = abs(tau)
  end subroutine
end module mod_abstract_evo_operator
