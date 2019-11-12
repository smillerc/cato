module mod_first_order_reconstruction
  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_conserved_vars, only: conserved_vars_t
  use mod_regular_2d_grid, only: regular_2d_grid_t

  implicit none

  private
  public :: first_order_reconstruction_t

  type, extends(abstract_reconstruction_t) :: first_order_reconstruction_t
  contains
    procedure, public :: initialize => init_first_order
    procedure, public :: reconstruct => reconstruct_first_order
  end type

contains
  subroutine init_first_order(self)
    class(first_order_reconstruction_t), intent(inout) :: self

  end subroutine init_first_order

  pure subroutine reconstruct_first_order(self, U, grid, i, j, U_bar)
    class(first_order_reconstruction_t), intent(in) :: self
    class(regular_2d_grid_t), intent(in) :: grid
    class(conserved_vars_t), intent(in) :: U
    integer(ik), intent(in) :: i, j
    real(rk), dimension(4), intent(out) :: U_bar  !< Approximate reconstructed [rho, u, v, p]

  end subroutine
end module mod_first_order_reconstruction
