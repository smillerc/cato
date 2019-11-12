module mod_abstract_reconstruction

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_conserved_vars, only: conserved_vars_t
  use mod_regular_2d_grid, only: regular_2d_grid_t
  use mod_slope_limiter, only: slope_limiter_t
  ! use slope_limiter, only: limit

  implicit none

  private
  public :: abstract_reconstruction_t

  type, abstract :: abstract_reconstruction_t
    integer(ik) :: order
    character(:), allocatable :: name
    class(slope_limiter_t), allocatable :: limiter
  contains
    procedure(initialize), public, deferred :: initialize
    procedure(reconstruct), public, deferred :: reconstruct
  end type abstract_reconstruction_t

  abstract interface
    subroutine initialize(self)
      import :: abstract_reconstruction_t
      class(abstract_reconstruction_t), intent(inout) :: self
    end subroutine initialize

    pure subroutine reconstruct(self, U, grid, i, j, U_bar)
      import :: abstract_reconstruction_t
      import :: conserved_vars_t
      import :: regular_2d_grid_t
      import :: ik, rk
      class(abstract_reconstruction_t), intent(in) :: self
      class(regular_2d_grid_t), intent(in) :: grid
      class(conserved_vars_t), intent(in) :: U
      integer(ik), intent(in) :: i, j
      real(rk), dimension(4), intent(out) :: U_bar  !< U_bar = [rho, u, v, p]
    end subroutine reconstruct
  end interface

contains

  ! pure function reconstruct_piecewise() result(W)
  !   !< Follows the reconstruction scheme of DOI: 10.1016/j.jcp.2009.04.001, although it
  !   !< inherits it from DOI: 10.1016/j.jcp.2006.03.018, where the explaination is a bit clearer,
  !   !< and without typos
  ! end function

end module mod_abstract_reconstruction
