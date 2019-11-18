module mod_second_order_reconstruction
  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_conserved_vars, only: conserved_vars_t
  use mod_grid, only: grid_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t

  implicit none

  private
  public :: second_order_reconstruction_t

  type, extends(abstract_reconstruction_t) :: second_order_reconstruction_t
    private
    !< 2nd order (piecewise linear) reconstruction operator
    ! class(grid_t), pointer :: grid !< grid type to hold control volume info
    real(rk), dimension(4, 2) :: cell_gradient = 0.0_rk !< approximated cell gradient (unique to each cell)
  contains
    procedure, public :: initialize => init_second_order
    procedure, public :: reconstruct => reconstruct_second_order
    procedure, public :: select_and_find_gradient
    procedure, private :: estimate_gradients
    procedure, private :: estimate_single_gradient
    ! procedure, private :: finalize
    ! final :: finalize_second_order
  end type

  ! interface second_order_reconstruction_t
  !   module procedure :: constructor
  ! end interface

contains

  subroutine init_second_order(self, input)
    !< Construct the second_order_reconstruction_t type
    class(second_order_reconstruction_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    ! self%cell_is_selected = .false.
    ! self%current_cell_ij = [0,0]
    ! self%cell_gradient = 0.0_rk
    ! self%cell_average = 0.0_rk
    self%order = 2
    self%name = 'piecewise_linear_reconstruction'
    call self%set_slope_limiter(name=input%slope_limiter)

  end subroutine init_second_order

  ! subroutine finalize(self)
  !   type(second_order_reconstruction_t), intent(inout) :: self
  !   call self%finalize_second_order
  ! end subroutine finalize

  ! subroutine finalize_second_order(self)
  !   type(second_order_reconstruction_t), intent(inout) :: self
  !   ! nullify(self%limiter)
  ! end subroutine finalize_second_order

  subroutine select_and_find_gradient(self, i, j, grid, conserved_vars)
    !< This sets the interface to select the cell to reconstruct and calculate the gradient. Because the
    !< gradient is reused, this is stored in the type to be used later by the reconstruct function
    class(second_order_reconstruction_t), intent(inout) :: self
    class(grid_t), intent(in) :: grid
    class(conserved_vars_t), intent(in) :: conserved_vars
    integer(ik) :: i, j

    call self%find_cell_average(i, j, conserved_vars)
    self%cell_gradient = self%estimate_gradients(i, j, conserved_vars, grid)
    self%cell_is_selected = .true.
  end subroutine select_and_find_gradient

  pure function reconstruct_second_order(self, x, y, grid) result(U_bar)
    !< Reconstruct the value of the conserved variables (U) in a cell (i,j) at location (x,y) based on the
    !> cell average and gradient.

    class(second_order_reconstruction_t), intent(in) :: self
    real(rk), intent(in) :: x, y !< where should U_bar be reconstructed at?
    class(grid_t), intent(in) :: grid
    real(rk), dimension(4) :: U_bar  !< U_bar = [rho, u, v, p]

    real(rk), dimension(2) :: centroid_xy
    real(rk) :: grad_rho, grad_u, grad_v, grad_p

   if(.not. self%cell_is_selected) error stop "Cell i,j was never selected in the reconstruction procedure! (this shouldn't happen)"

    centroid_xy = grid%get_cell_centroid_xy(i=self%current_cell_ij(1), j=self%current_cell_ij(1))

    associate(dV_dx=>self%cell_gradient(1, :), dV_dy=>self%cell_gradient(2, :), &
              x_ij=>centroid_xy(1), y_ij=>centroid_xy(2))

      U_bar = self%cell_average + dV_dx * (x - x_ij) + dV_dy * (y - y_ij)
    end associate
  end function reconstruct_second_order

  pure function estimate_gradients(self, i, j, conserved_vars, grid) result(gradients)
    class(second_order_reconstruction_t), intent(in) :: self
    class(grid_t), intent(in) :: grid
    class(conserved_vars_t), intent(in) :: conserved_vars
    real(rk), dimension(4, 2) :: gradients
    integer(ik), intent(in) :: i, j

    associate(rho=>conserved_vars%density, p=>conserved_vars%pressure, &
              u=>conserved_vars%x_velocity, v=>conserved_vars%y_velocity)

      ! density
      gradients(1, :) = self%estimate_single_gradient(rho, i, j, grid)

      ! x velocity
      gradients(2, :) = self%estimate_single_gradient(u, i, j, grid)

      ! y velocity
      gradients(3, :) = self%estimate_single_gradient(v, i, j, grid)

      ! pressure
      gradients(4, :) = self%estimate_single_gradient(p, i, j, grid)
    end associate

  end function estimate_gradients

  pure function estimate_single_gradient(self, v, i, j, grid) result(grad_v)
    !< Find the gradient of a variable (v) within a cell at indices (i,j) based on the neighbor information.
    !< See Eq. 9 in https://doi.org/10.1016/j.jcp.2006.03.018. The slope limiter is set via the constructor
    !< of this derived type.

    class(second_order_reconstruction_t), intent(in) :: self
    class(grid_t), intent(in) :: grid
    real(rk), dimension(:, :), intent(in) :: v !< variable to estimate the gradient
    integer(ik), intent(in) :: i, j
    real(rk), dimension(2) :: grad_v

    ! associate(limit => self%limiter%limit, &
    !           volume => grid%get_cell_volumes(i,j), &
    !           n1 => grid%cell_edge_norm_vectors(i, j, 1, :), &
    !           n2 => grid%cell_edge_norm_vectors(i, j, 2, :), &
    !           n3 => grid%cell_edge_norm_vectors(i, j, 3, :), &
    !           n4 => grid%cell_edge_norm_vectors(i, j, 4, :), &
    !           delta_l1 => grid%edge_lengths(i, j, 1), &
    !           delta_l2 => grid%edge_lengths(i, j, 2), &
    !           delta_l3 => grid%edge_lengths(i, j, 3), &
    !           delta_l4 => grid%edge_lengths(i, j, 4))

    !   grad_v = (1._rk / (2.0_rk * volume)) * &
    !            (limit(v(i + 1, j) - v(i, j), v(i, j) - v(i - 1, j)) * (n2 * delta_l2 - n4 * delta_l4) + &
    !             limit(v(i, j + 1) - v(i, j), v(i, j) - v(i, j - 1)) * (n3 * delta_l3 - n1 * delta_l1))
    ! end associate

  end function estimate_single_gradient

end module mod_second_order_reconstruction
