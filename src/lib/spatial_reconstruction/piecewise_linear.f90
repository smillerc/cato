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

module mod_piecewise_linear_reconstruction

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print, n_ghost_layers
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_grid_block, only: grid_block_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use mod_floating_point_utils, only: equal
  use mod_gradients, only: green_gauss_gradient
  use mod_edge_interpolator, only: edge_iterpolator_t
  use mod_edge_interpolator_factory, only: edge_interpolator_factory

  implicit none

  private
  public :: piecewise_linear_reconstruction_t

  type, extends(abstract_reconstruction_t) :: piecewise_linear_reconstruction_t
    !< Implementation of a 2nd order piecewise-linear reconstruction operator
    class(edge_iterpolator_t), allocatable :: edge_interpolator !< how are the edges interpolated?, e.g. TVD2, MLP3, etc.
  contains
    procedure, public :: initialize
    procedure, public :: reconstruct
    procedure, public :: reconstruct_at_point
    final :: finalize
  end type

contains

  subroutine initialize(self, input, grid_target)
    !< Construct the piecewise_linear_reconstruction_t type

    class(piecewise_linear_reconstruction_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid_target

    class(edge_iterpolator_t), pointer :: edge_interpolator => null()

    call debug_print('Initializing piecewise_linear_reconstruction_t', __FILE__, __LINE__)

    self%order = 2
    self%name = 'piecewise_linear'
    self%grid => grid_target
    edge_interpolator => edge_interpolator_factory(input)
    allocate(self%edge_interpolator, source=edge_interpolator)
    deallocate(edge_interpolator)

  end subroutine initialize

  subroutine finalize(self)
    !< Finalize the piecewise_linear_reconstruction_t type
    type(piecewise_linear_reconstruction_t), intent(inout) :: self

    call debug_print('Running piecewise_linear_reconstruction_t%finalize()', __FILE__, __LINE__)

    if(allocated(self%edge_interpolator)) deallocate(self%edge_interpolator)

    if(associated(self%grid)) nullify(self%grid)
    if(associated(self%rho)) nullify(self%rho)
    if(associated(self%p)) nullify(self%p)

    if(allocated(self%grad_x_rho)) deallocate(self%grad_x_rho)
    if(allocated(self%grad_x_p)) deallocate(self%grad_x_p)
    if(allocated(self%grad_y_rho)) deallocate(self%grad_y_rho)
    if(allocated(self%grad_y_p)) deallocate(self%grad_y_p)

  end subroutine finalize

  real(rk) function reconstruct_at_point(self, i, j, x, y, var) result(q)
    class(piecewise_linear_reconstruction_t), intent(in) :: self
    real(rk), intent(in) :: x, y  !< location within cell
    integer(ik), intent(in) :: i, j !< cell indices
    character(len=*), intent(in) :: var !< variable to reconstruct ('rho', or 'p')

    real(rk) :: x_ij !< cell centroid x location
    real(rk) :: y_ij !< cell centroid y location

    real(rk) :: q_x, q_y
    real(rk) :: q_bar !< cell average quantity
    real(rk) :: grad_x, grad_y

    x_ij = self%grid%cell_centroid_x(i, j)
    y_ij = self%grid%cell_centroid_y(i, j)

    select case(trim(var))
    case('rho')
      if(.not. associated(self%rho)) then
        error stop "Error in piecewise_linear_reconstruction_t%reconstruct_at_point(), self%rho isn't associated!"
      end if

      if(.not. allocated(self%grad_x_rho) .or. .not. allocated(self%grad_y_rho)) then
        error stop "Error in piecewise_linear_reconstruction_t%reconstruct_at_point(), grad_rho isn't allocated!"
      end if
      q_bar = self%rho(i, j)
      grad_x = self%grad_x_rho(i, j)
      grad_y = self%grad_y_rho(i, j)

    case('p')
      if(.not. associated(self%p)) then
        error stop "Error in piecewise_linear_reconstruction_t%reconstruct_at_point(), self%p isn't associated!"
      end if

      if(.not. allocated(self%grad_x_p) .or. .not. allocated(self%grad_y_p)) then
        error stop "Error in piecewise_linear_reconstruction_t%reconstruct_at_point(), grad_p isn't allocated!"
      end if
      q_bar = self%p(i, j)
      grad_x = self%grad_x_p(i, j)
      grad_y = self%grad_y_p(i, j)

      ! write(*, '(3(a,es16.6))') 'ave p: ', q_bar, ' d/dx: ', grad_x, ' d/dy: ', grad_y
    case default
      error stop "Error in piecewise_linear_reconstruction_t%reconstruct_at_point(), var must be 'p' or 'rho'"
    end select

    q_x = grad_x * (x - x_ij)
    q_y = grad_y * (y - y_ij)

    if(abs(q_x) < 1e-9_rk) q_x = 0.0_rk
    if(abs(q_y) < 1e-9_rk) q_y = 0.0_rk

    q = q_bar + (q_x + q_y)
    if(abs(q - q_bar) < 1e-12_rk) q = q_bar
  end function reconstruct_at_point

  subroutine reconstruct(self, primitive_var, reconstructed_var, lbounds, name, stage_name)
    !< Reconstruct each corner/midpoint. This converts the cell centered conserved
    !< quantities [rho, rho u, rho v, e] to reconstructed primitive variables [rho, u, v, p]
    !< based on the chosen reconstruction order, e.g. using a piecewise-linear function based on the
    !< selected cell and it's neighbors. Rather than do it a point at a time, this reuses some
    !< of the data necessary, like the cell average and gradient
    class(piecewise_linear_reconstruction_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    character(len=*), intent(in) :: stage_name
    character(len=*), intent(in) :: name
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: primitive_var !< (i,j); cell primitive variable to reconstruct
    real(rk), dimension(:, lbounds(1):, lbounds(2):), contiguous, intent(out) :: reconstructed_var
    !< ((corner1:midpoint4), i, j); reconstructed variable, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint

    real(rk), dimension(:, :, :), allocatable :: edge_values !< ((edge 1:n), i, j); interpolated edge/cell interface values
    real(rk), dimension(:, :), allocatable :: grad_x !< (i,j); x-gradient of the primitive variable
    real(rk), dimension(:, :), allocatable :: grad_y !< (i,j); y-gradient of the primitive variable

    integer(ik) :: i, j, p  !< cell i,j index
    integer(ik) :: ilo, ihi, jlo, jhi !< grid bounds w/o ghost layers
    integer(ik) :: ilo_bc, ihi_bc, jlo_bc, jhi_bc !< grid bounds w/ ghost layers
    real(rk) :: x, y  !< reconstruction location
    real(rk) :: dx, dy
    real(rk) :: recon_var
    real(rk) :: x_ij, y_ij !< cell centroin location
    integer(ik), dimension(3) :: edge_lbounds

    call debug_print('Running piecewise_linear_reconstruction_t%reconstruct()', __FILE__, __LINE__)

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

    ! Reconstruct the values at the cell interfaces
    call self%edge_interpolator%interpolate_edge_values(q=primitive_var, &
                                                        lbounds=lbounds, &
                                                        edge_values=edge_values)

    ! Now find the cell gradient
    edge_lbounds = lbound(edge_values)
    allocate(grad_x(ilo_bc:ihi_bc, jlo_bc:jhi_bc))
    allocate(grad_y(ilo_bc:ihi_bc, jlo_bc:jhi_bc))
    grad_x = 0.0_rk
    grad_y = 0.0_rk
    call green_gauss_gradient(edge_vars=edge_values, lbounds=edge_lbounds, grid=self%grid, &
                              grad_x=grad_x, grad_y=grad_y, name=name, stage_name=stage_name)

    select case(trim(name))
    case('rho')
      if(.not. allocated(self%grad_x_rho)) allocate(self%grad_x_rho, mold=grad_x)
      if(.not. allocated(self%grad_y_rho)) allocate(self%grad_y_rho, mold=grad_y)
      self%grad_x_rho(:, :) = grad_x
      self%grad_y_rho(:, :) = grad_y
    case('p')
      if(.not. allocated(self%grad_x_p)) allocate(self%grad_x_p, mold=grad_x)
      if(.not. allocated(self%grad_y_p)) allocate(self%grad_y_p, mold=grad_y)
      self%grad_x_p(:, :) = grad_x
      self%grad_y_p(:, :) = grad_y
    end select

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi) &
    !$omp private(i, j, x, y, x_ij, y_ij, dx, dy, recon_var) &
    !$omp shared(reconstructed_var, self, grad_x, grad_y, primitive_var)
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi

        x_ij = self%grid%cell_centroid_x(i, j)
        y_ij = self%grid%cell_centroid_y(i, j)

        do p = 1, 8
          x = self%grid%cell_node_x(p, i, j)
          y = self%grid%cell_node_y(p, i, j)
          dx = grad_x(i, j) * (x - x_ij)
          dy = grad_y(i, j) * (y - y_ij)

          ! write(*, '(4(es16.6))') dx, dy, grad_x(i, j), grad_y(i, j)
          if(abs(dx) < 1e-9_rk) dx = 0.0_rk
          if(abs(dy) < 1e-9_rk) dy = 0.0_rk

          recon_var = primitive_var(i, j) + (dx + dy)

          if(abs(recon_var - primitive_var(i, j)) < 1e-12_rk) then
            reconstructed_var(p, i, j) = primitive_var(i, j)
          else
            reconstructed_var(p, i, j) = recon_var
          end if

        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(grad_x)
    deallocate(grad_y)
    deallocate(edge_values)
  end subroutine reconstruct
end module mod_piecewise_linear_reconstruction
