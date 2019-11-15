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
    logical :: cell_is_selected = .false.
    class(slope_limiter_t), allocatable :: limiter
    class(conserved_vars_t), pointer :: conserved_vars !< aka U, the vector of conserved variables
    real(rk), dimension(4) :: cell_average = 0.0_rk !< average values of the cell conserved variables
  contains
    procedure, private :: find_cell_average
    procedure(reconstruct), public, deferred :: reconstruct
  end type abstract_reconstruction_t

  abstract interface

    pure function reconstruct(self, x, y, i, j) result(U_bar)
      !< Reconstruct the value of the conserved variables (U) in a cell (i,j) at location (x,y) based on the
      !> cell average and gradient.
      import :: abstract_reconstruction_t
      import :: ik, rk

      class(abstract_reconstruction_t), intent(in) :: self
      real(rk), dimension(4) :: U_bar  !< U_bar = [rho, u, v, p]
      real(rk), intent(in) :: x, y !< where should U_bar be reconstructed at?
      integer(ik), intent(in) :: i, j
    end function reconstruct
  end interface

contains

  pure function find_cell_average(self, i, j) result(U_bar)
    class(abstract_reconstruction_t), intent(in) :: self
    real(rk), dimension(4) :: U_bar
    integer(ik), intent(in) :: i, j

    associate(rho=>self%conserved_vars%density, p=>self%conserved_vars%pressure, &
              u=>self%conserved_vars%x_velocity, v=>self%conserved_vars%y_velocity)

      ! density
      U_bar(1) = sum(rho(i - 1:i + 1, j - 1:j + 1)) / 5.0_rk

      ! x velocity
      U_bar(2) = sum(u(i - 1:i + 1, j - 1:j + 1)) / 5.0_rk

      ! y velocity
      U_bar(3) = sum(v(i - 1:i + 1, j - 1:j + 1)) / 5.0_rk

      ! pressure
      U_bar(4) = sum(p(i - 1:i + 1, j - 1:j + 1)) / 5.0_rk
    end associate

  end function

end module mod_abstract_reconstruction
