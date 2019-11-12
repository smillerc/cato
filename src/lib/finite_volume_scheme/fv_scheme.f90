module mod_finite_volume_schemes

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_surrogate, only: surrogate
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand

  implicit none
  private
  public :: finite_volume_scheme_t

  type, abstract, extends(integrand) :: finite_volume_scheme_t
    class(abstract_reconstruction_t), allocatable :: reconstruction_operator
  contains
    procedure(reconstruct), public, deferred :: reconstruct
    procedure(eval_fluxes), public, deferred :: eval_fluxes
  end type

  abstract interface
    subroutine reconstruct(self)
      import :: finite_volume_scheme_t
      class(finite_volume_scheme_t), intent(inout) :: self
    end subroutine

    subroutine eval_fluxes(self)
      import :: finite_volume_scheme_t
      class(finite_volume_scheme_t), intent(inout) :: self
    end subroutine

  end interface
contains

end module mod_finite_volume_schemes
