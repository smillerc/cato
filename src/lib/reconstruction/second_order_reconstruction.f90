module mod_second_order_reconstruction
  use iso_fortran_env, only: int32, real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_conserved_vars, only: conserved_vars_t

  implicit none

  private
  public :: second_order_reconstruction_t

  type, extends(abstract_reconstruction_t) :: second_order_reconstruction_t
    class(slope_limiter_t), allocatable :: limiter
  contains
    procedure, public :: initialize => init_second_order
    procedure, public :: reconstruct => reconstruct_second_order
  end type

contains

  subroutine init_second_order(self)
    class(second_order_reconstruction_t), intent(inout) :: self

  end subroutine init_second_order

  pure subroutine reconstruct_second_order(self, U, i, j, U_bar)
    class(second_order_reconstruction_t), intent(in) :: self
    class(conserved_vars_t), intent(in) :: U
    integer(int32), intent(in) :: i, j
    real(real64), dimension(4), intent(out) :: U_bar  !< Approximate reconstructed [rho, u, v, p]

  end subroutine

end module mod_second_order_reconstruction
