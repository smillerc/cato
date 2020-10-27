! MIT License
! Copyright (c) 2019 Sam Miller
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module mod_piecewise_constant_reconstruction
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit, std_out => output_unit
  use mod_globals, only: debug_print, n_ghost_layers
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_grid_block, only: grid_block_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t

  implicit none

  private
  public :: piecewise_constant_reconstruction_t

  type, extends(abstract_reconstruction_t) :: piecewise_constant_reconstruction_t
  contains
    procedure, public :: initialize => init_first_order
    procedure, public :: reconstruct
    procedure, public :: reconstruct_at_point
    final :: finalize
  end type

contains

  subroutine init_first_order(self, input, grid_target)
    !< Construct the piecewise_constant_reconstruction_t type

    class(piecewise_constant_reconstruction_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid_target

    integer(ik) :: alloc_status

    call debug_print('Initializing second_order_reconstruction_t', __FILE__, __LINE__)

    self%order = 1
    self%name = 'cell_average_reconstruction'

    self%grid => grid_target

  end subroutine init_first_order

  subroutine finalize(self)
    !< Finalize the piecewise_constant_reconstruction_t type
    type(piecewise_constant_reconstruction_t), intent(inout) :: self
    integer(ik) :: alloc_status

    call debug_print('Running piecewise_constant_reconstruction_t%finalize()', __FILE__, __LINE__)

    if(associated(self%grid)) nullify(self%grid)
    if(associated(self%rho)) nullify(self%rho)
    if(associated(self%p)) nullify(self%p)

    if(allocated(self%grad_x_rho)) deallocate(self%grad_x_rho)
    if(allocated(self%grad_x_p)) deallocate(self%grad_x_p)
    if(allocated(self%grad_y_rho)) deallocate(self%grad_y_rho)
    if(allocated(self%grad_y_p)) deallocate(self%grad_y_p)

  end subroutine finalize

  pure real(rk) function reconstruct_at_point(self, i, j, x, y, var) result(q)
    class(piecewise_constant_reconstruction_t), intent(in) :: self
    real(rk), intent(in) :: x, y  !< location within cell
    integer(ik), intent(in) :: i, j !< cell indices
    character(len=*), intent(in) :: var !< variable to reconstruct ('rho', or 'p')

    select case(trim(var))
    case('rho')
      if(.not. associated(self%rho)) then
        error stop "Error in piecewise_linear_reconstruction_t%reconstruct_at_point(), self%rho isn't associated!"
      end if
      q = self%rho(i, j)

    case('p')
      if(.not. associated(self%p)) then
        error stop "Error in piecewise_linear_reconstruction_t%reconstruct_at_point(), self%p isn't associated!"
      end if
      q = self%p(i, j)
    case default
      error stop "Error in piecewise_linear_reconstruction_t%reconstruct_at_point(), var must be 'p' or 'rho'"
    end select
  end function

  subroutine reconstruct(self, primitive_var, reconstructed_var, lbounds, name, stage_name)
    !< Reconstruct the entire domain. Rather than do it a point at a time, this reuses some
    !< of the data necessary, like the cell average and gradient

    class(piecewise_constant_reconstruction_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    character(len=*), intent(in) :: stage_name
    character(len=*), intent(in) :: name
    real(rk), dimension(lbounds(1):, lbounds(2):), intent(in), contiguous :: primitive_var !< (i,j); cell primitive variable to reconstruct
    real(rk), dimension(:, lbounds(1):, lbounds(2):), intent(out), contiguous :: reconstructed_var
    !< ((corner1:midpoint4), i, j); reconstructed variable, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint

    integer(ik) :: i, j, p
    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: ilo_bc, ihi_bc, jlo_bc, jhi_bc

    if(.not. associated(self%grid)) error stop "Grid not associated"

    ! Bounds do not include ghost cells. Ghost cells get their
    ! reconstructed values and gradients from the boundary conditions
    ilo_bc = lbound(primitive_var, dim=1)
    ihi_bc = ubound(primitive_var, dim=1)
    jlo_bc = lbound(primitive_var, dim=2)
    jhi_bc = ubound(primitive_var, dim=2)

    ilo = ilo_bc + n_ghost_layers
    ihi = ihi_bc - n_ghost_layers
    jlo = jlo_bc + n_ghost_layers
    jhi = jhi_bc - n_ghost_layers

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j) &
    !$omp shared(reconstructed_var, primitive_var)
    !$omp do simd
    do j = jlo, jhi
      do i = ilo, ihi
        do p = 1, 8
          reconstructed_var(p, i, j) = primitive_var(i, j)
        end do
      end do
    end do
    !$omp end do simd
    !$omp end parallel

  end subroutine reconstruct

end module mod_piecewise_constant_reconstruction
