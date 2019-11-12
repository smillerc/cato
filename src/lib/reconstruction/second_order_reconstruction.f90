module mod_second_order_reconstruction
  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_conserved_vars, only: conserved_vars_t
  use mod_regular_2d_grid, only: regular_2d_grid_t

  implicit none

  private
  public :: second_order_reconstruction_t

  type, extends(abstract_reconstruction_t) :: second_order_reconstruction_t

  contains
    procedure, public :: initialize => init_second_order
    procedure, public :: reconstruct => reconstruct_second_order
  end type

contains

  subroutine init_second_order(self)
    class(second_order_reconstruction_t), intent(inout) :: self

  end subroutine init_second_order

  pure subroutine reconstruct_second_order(self, U, grid, i, j, U_bar)
    class(second_order_reconstruction_t), intent(in) :: self
    class(regular_2d_grid_t), intent(in) :: grid
    class(conserved_vars_t), intent(in) :: U
    integer(ik), intent(in) :: i, j
    real(rk), dimension(4), intent(out) :: U_bar  !< Approximate reconstructed [rho, u, v, p]

  end subroutine

  !<  Numbering convention for a 2D quadrilateral cell
  !<
  !<                     F3
  !<                     M3
  !<              N4-----o-----N3
  !<              |            |
  !<      F4   M4 o      C     o M2   F2
  !<              |            |
  !<              N1-----o-----N2
  !<                     M1
  !<                     F1
  !< N: node or vertex
  !< F: face or edge
  !< M: midpoint of the edge (o)
  !< C: cell or control volume (in finite-volume lingo)

  ! pure function estimate_gradient(self,grid,v,i,j) result(grad_v)
  !   class(second_order_reconstruction_t), intent(in) :: self
  !   class(regular_2d_grid_t), intent(in) :: grid
  !   real(rk), dimension(2,2) :: grad_v

  !   associate(limit => self%limiter%limit, &
  !             volume => grid%elements(i,j)%volume, &
  !             n1 => grid%elements(i,j)%edge_norm_vectors(1,:,:), &
  !             n2 => grid%elements(i,j)%edge_norm_vectors(2,:,:), &
  !             n3 => grid%elements(i,j)%edge_norm_vectors(3,:,:), &
  !             n4 => grid%elements(i,j)%edge_norm_vectors(4,:,:), &
  !             delta_l1 => grid%elements(i,j)%edge_lengths(1), &
  !             delta_l2 => grid%elements(i,j)%edge_lengths(2), &
  !             delta_l3 => grid%elements(i,j)%edge_lengths(3), &
  !             delta_l4 => grid%elements(i,j)%edge_lengths(4) )

  !     grad_v = (1._rk / (2*volume)) * &
  !               ( limit(v(i+1,j) - v(i,j), v(i,j) - v(i-1,j)) * (n2*delta_l2 - n4*delta_l4) + &
  !                 limit(v(i,j+1) - v(i,j), v(i,j) - v(i,j-1)) * (n3*delta_l3 - n1*delta_l1))
  !   end associate

  ! end function

end module mod_second_order_reconstruction
