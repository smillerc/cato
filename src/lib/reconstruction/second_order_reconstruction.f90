module mod_second_order_reconstruction
  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_conserved_vars, only: conserved_vars_t
  use mod_grid, only: grid_t
  use mod_slope_limiter, only: slope_limiter_t

  implicit none

  private
  public :: second_order_reconstruction_t

  type, extends(abstract_reconstruction_t) :: second_order_reconstruction_t
    private
    !< 2nd order (piecewise linear) reconstruction operator
    class(grid_t), pointer :: grid !< grid type to hold control volume info

    real(rk), dimension(4, 2) :: cell_gradient = 0.0_rk !< approximated cell gradient (unique to each cell)
  contains
    procedure, public :: reconstruct => reconstruct_second_order
    procedure, public :: select_cell_to_reconstruct
    procedure, private :: estimate_gradients
    procedure, private :: estimate_single_gradient
    final :: finalize
  end type

  interface second_order_reconstruction_t
    module procedure :: constructor
  end interface

contains

  pure function constructor(grid, conserved_vars, slope_limiter) result(operator)
    !< Construct the second_order_reconstruction_t type
    class(grid_t), intent(in), target :: grid !< grid type to hold control volume info
    class(conserved_vars_t), intent(in), target :: conserved_vars !< aka U, the vector of conserved variables
    class(slope_limiter_t), intent(in), target :: slope_limiter !< The sloper limiter to be used in the reconstruction process
    class(second_order_reconstruction_t) :: operator

    operator%grid => grid
    operator%conserved_vars => conserved_vars
    operator%limiter => slope_limiter

  end function constructor

  pure function finalize(self)
    type(second_order_reconstruction_t), intent(in) :: self
    nullify(self%grid)
    nullify(self%conserved_vars)
    nullify(self%limiter)
    self%cell_is_selected = .false.
    self%cell_gradient = 0.0_rk
    self%cell_average = 0.0_rk
  end function finalize

  subroutine select_cell_to_reconstruct(self, i, j)
    !< This sets the interface to select the cell to reconstruct and calculate the gradient. Because the
    !< gradient is reused, this is stored in the type to be used later by the reconstruct function
    class(second_order_reconstruction_t), intent(inout) :: self

    self%cell_average = self%find_cell_average(i, j)
    self%cell_gradient = self%estimate_gradients(i, j)
    self%cell_is_selected = .true.
  end subroutine

  pure function reconstruct(self, x, y, i, j) result(U_bar)
    !< Reconstruct the value of the conserved variables (U) in a cell (i,j) at location (x,y) based on the
    !> cell average and gradient.

    class(second_order_reconstruction_t), intent(in) :: self
    real(rk), dimension(4) :: U_bar  !< U_bar = [rho, u, v, p]
    real(rk), intent(in) :: x, y !< where should U_bar be reconstructed at?
    integer(ik), intent(in) :: i, j

    real(rk) :: grad_rho, grad_u, grad_v, grad_p

   if(.not. self%cell_is_selected) error stop "Cell i,j was never selected in the reconstruction procedure! (this shouldn't happen)"

    associate(dV_dx=>self%cell_gradient(1, :), dV_dy=>self%cell_gradient(2, :), &
              x_ij=>self%grid%cell_centroid_xy(i, j, 1), y_ij=>self%grid%cell_centroid_xy(i, j, 2))

      U_bar = self%cell_average + dV_dx * (x - x_ij) + dV_dy * (y - y_ij)
    end associate
  end function reconstruct

  pure function estimate_gradients(self, i, j) result(gradients)
    class(second_order_reconstruction_t), intent(in) :: self
    real(rk), dimension(4, 2) :: gradients
    integer(ik), intent(in) :: i, j

    associate(rho=>self%conserved_vars%density, p=>self%conserved_vars%pressure, &
              u=>self%conserved_vars%x_velocity, v=>self%conserved_vars%y_velocity)

      ! density
      gradients(1, :, :) = self%estimate_single_gradient(rho, i, j)

      ! x velocity
      gradients(2, :, :) = self%estimate_single_gradient(u, i, j)

      ! y velocity
      gradients(3, :, :) = self%estimate_single_gradient(v, i, j)

      ! pressure
      gradients(4, :, :) = self%estimate_single_gradient(p, i, j)
    end associate

  end function

  pure function estimate_single_gradient(self, v, i, j) result(grad_v)
    !< Find the gradient of a variable (v) within a cell at indices (i,j) based on the neighbor information.
    !< See Eq. 9 in https://doi.org/10.1016/j.jcp.2006.03.018. The slope limiter is set via the constructor
    !< of this derived type.

    class(second_order_reconstruction_t), intent(in) :: self
    real(rk), dimension(:, :) :: v !< variable to estimate the gradient
    integer(ik), intent(in) :: i, j
    real(rk), dimension(2, 2) :: grad_v

    associate(limit=>self%limiter%limit, &
              volume=>, &
              n1=>self%grid%cell_edge_norm_vectors(i, j, 1, :), &
              n2=>self%grid%cell_edge_norm_vectors(i, j, 2, :), &
              n3=>self%grid%cell_edge_norm_vectors(i, j, 3, :), &
              n4=>self%grid%cell_edge_norm_vectors(i, j, 4, :), &
              delta_l1=>self%grid%edge_lengths(i, j, 1), &
              delta_l2=>self%grid%edge_lengths(i, j, 2), &
              delta_l3=>self%grid%edge_lengths(i, j, 3), &
              delta_l4=>self%grid%edge_lengths(i, j, 4))

      grad_v = (1._rk / (2.0_rk * volume)) * &
               (limit(v(i + 1, j) - v(i, j), v(i, j) - v(i - 1, j)) * (n2 * delta_l2 - n4 * delta_l4) + &
                limit(v(i, j + 1) - v(i, j), v(i, j) - v(i, j - 1)) * (n3 * delta_l3 - n1 * delta_l1))
    end associate

  end function

end module mod_second_order_reconstruction
