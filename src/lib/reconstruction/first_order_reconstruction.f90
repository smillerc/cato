module mod_first_order_reconstruction
!   use iso_fortran_env, only: ik => int32, rk => real64
!   use mod_abstract_reconstruction, only: abstract_reconstruction_t
!   use mod_conserved_vars, only: conserved_vars_t
!   use mod_regular_2d_grid, only: regular_2d_grid_t
!   use mod_input, only: input_t
!   use mod_grid, only: grid_t

!   implicit none

!   private
!   public :: first_order_reconstruction_t

!   type, extends(abstract_reconstruction_t) :: first_order_reconstruction_t
!   contains
!     procedure, public :: initialize => init_first_order
!     procedure, public :: reconstruct => reconstruct_first_order
!     procedure, public :: select_and_find_gradient
!     ! procedure, private :: finalize
!     ! final :: finalize_first_order
!   end type

! contains
!   subroutine init_first_order(self, input)
!     class(first_order_reconstruction_t), intent(inout) :: self
!     class(input_t), intent(in) :: input

!     ! self%cell_is_selected = .false.
!     ! self%current_cell_ij = [0,0]
!     ! self%cell_gradient = 0.0_rk
!     ! self%cell_average = 0.0_rk
!     self%order = 1
!     self%name = 'cell_average'

!   end subroutine init_first_order

!   subroutine select_and_find_gradient(self, i, j, grid, conserved_vars)
!     !< This sets the interface to select the cell to reconstruct and calculate the gradient. Because the
!     !< gradient is reused, this is stored in the type to be used later by the reconstruct function
!     class(first_order_reconstruction_t), intent(inout) :: self
!     class(grid_t), intent(in) :: grid
!     class(conserved_vars_t), intent(in) :: conserved_vars
!     integer(ik) :: i, j

!     call self%find_cell_average(i, j, conserved_vars)
!     self%cell_is_selected = .true.
!   end subroutine select_and_find_gradient

!   pure function reconstruct_first_order(self, x, y, grid) result(U_bar)
!     !< Reconstruct the value of the conserved variables (U) in a cell (i,j) at location (x,y) based on the
!     !> cell average.

!     class(first_order_reconstruction_t), intent(in) :: self
!     real(rk), intent(in) :: x, y !< where should U_bar be reconstructed at?
!     class(grid_t), intent(in) :: grid

!     real(rk), dimension(4) :: U_bar  !< U_bar = [rho, u, v, p]

!     real(rk) :: grad_rho, grad_u, grad_v, grad_p

!    if(.not. self%cell_is_selected) error stop "Cell i,j was never selected in the reconstruction procedure! (this shouldn't happen)"

!     U_bar = self%cell_average
!   end function reconstruct_first_order

!   ! subroutine finalize(self)
!   !   type(first_order_reconstruction_t), intent(inout) :: self
!   !   call self%finalize_first_order()
!   ! end subroutine

!   ! subroutine finalize_first_order(self)
!   !   type(first_order_reconstruction_t), intent(inout) :: self
!   !   integer(ik) :: alloc_stat

!   !   self%cell_is_selected = .false.
!   !   ! if (allocated(self%slope_limiter)) then
!   !   !   deallocate(self%slope_limiter, stat=alloc_stat)
!   !   !   if (alloc_stat /= 0) error stop 'Unable to deallocate first_order_reconstruction_t%slope_limiter'
!   !   ! end if
!   ! end subroutine finalize_first_order
end module mod_first_order_reconstruction
