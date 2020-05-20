module mod_piecewise_linear_reconstruction

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print, n_ghost_layers
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_grid, only: grid_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use mod_floating_point_utils, only: equal
  use mod_gradients, only: green_gauss_gradient
  use mod_edge_reconstruction, only: reconstruct_edge_values

  implicit none

  private
  public :: piecewise_linear_reconstruction_t

  type, extends(abstract_reconstruction_t) :: piecewise_linear_reconstruction_t
    !< Implementation of a 2nd order piecewise-linear reconstruction operator
  contains
    procedure, public :: initialize
    procedure, public :: reconstruct
    final :: finalize
  end type

contains

  subroutine initialize(self, input, grid_target)
    !< Construct the piecewise_linear_reconstruction_t type

    class(piecewise_linear_reconstruction_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid_target

    integer(ik) :: alloc_status

    call debug_print('Initializing piecewise_linear_reconstruction_t', __FILE__, __LINE__)

    self%order = 2
    self%name = 'piecewise_linear'

    self%grid => grid_target

    call self%set_slope_limiter(name=input%limiter)

  end subroutine initialize

  subroutine finalize(self)
    !< Finalize the piecewise_linear_reconstruction_t type
    type(piecewise_linear_reconstruction_t), intent(inout) :: self
    integer(ik) :: alloc_status

    call debug_print('Running piecewise_linear_reconstruction_t%finalize()', __FILE__, __LINE__)

    if(associated(self%grid)) nullify(self%grid)
  end subroutine finalize

  subroutine reconstruct(self, primitive_var, reconstructed_var, lbounds)
    !< Reconstruct each corner/midpoint. This converts the cell centered conserved
    !< quantities [rho, rho u, rho v, e] to reconstructed primitive variables [rho, u, v, p]
    !< based on the chosen reconstruction order, e.g. using a piecewise-linear function based on the
    !< selected cell and it's neighbors. Rather than do it a point at a time, this reuses some
    !< of the data necessary, like the cell average and gradient
    class(piecewise_linear_reconstruction_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
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
    call reconstruct_edge_values(q=primitive_var, lbounds=lbounds, limiter=self%limiter, edge_values=edge_values)

    ! Now find the cell gradient
    edge_lbounds = lbound(edge_values)
    call green_gauss_gradient(edge_vars=edge_values, lbounds=edge_lbounds, grid=self%grid, &
                              grad_x=grad_x, grad_y=grad_y)

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
          if(dx < 1e-12_rk) dx = 0.0_rk
          if(dy < 1e-12_rk) dy = 0.0_rk

          recon_var = primitive_var(i, j) + dx + dy

          if(abs(recon_var - primitive_var(i, j)) < 1e-12_rk) then
            reconstructed_var(p, i, j) = primitive_var(i, j)
          else
            reconstructed_var(p, i, j) = recon_var
          end if

        end do
      end do
    end do

    ! write(*, '(a, es16.6)') "reconstructed_var(1, 204, 3) - reconstructed_var(1, 204, 2) -> ", reconstructed_var(1, 204, 3) - reconstructed_var(1, 204, 2)
    ! write(*, '(a, es16.6)') "grad_x(204, 3) - grad_x(204, 2) -> ", grad_x(204, 3) - grad_x(204, 2)
    ! write(*, '(a, es16.6)') "grad_y(204, 3) - grad_y(204, 2) -> ", grad_y(204, 3) - grad_y(204, 2)
    ! write(*, '(a, es16.6)') "cell_centroid_x(204, 3) - cell_centroid_x(204, 2) -> ", self%grid%cell_centroid_x(204, 3) - self%grid%cell_centroid_x(204, 2)
    ! write(*, '(a, es16.6)') "cell_centroid_y(204, 3) - cell_centroid_y(204, 2) -> ", self%grid%cell_centroid_y(204, 3) - self%grid%cell_centroid_y(204, 2)
    !$omp end do
    !$omp end parallel
  end subroutine reconstruct
end module mod_piecewise_linear_reconstruction
