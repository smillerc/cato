module mod_finite_volume_schemes

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_surrogate, only: surrogate
  use mod_strategy, only: strategy
  use mod_integrand, only: integrand_t
  use mod_conserved_vars, only: conserved_vars_t
  use mod_grid, only: grid_t

  implicit none
  private
  public :: finite_volume_scheme_t

  type, abstract, extends(integrand_t) :: finite_volume_scheme_t
    class(grid_t), allocatable :: grid
    class(conserved_vars_t), allocatable :: conserved_vars
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

    pure function eval_fluxes(self) result(H)
      import :: finite_volume_scheme_t
      real(rk), dimension(2, 4) :: H !< Flux tensor
      class(finite_volume_scheme_t), intent(inout) :: self
    end function

  end interface
contains

end module mod_finite_volume_schemes
