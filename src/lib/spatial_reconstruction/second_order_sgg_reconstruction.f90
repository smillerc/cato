module mod_second_order_sgg_reconstruction

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print, n_ghost_layers
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_grid, only: grid_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use mod_floating_point_utils, only: equal
  use mod_gradients, only: green_gauss_gradient, get_smoothness

  implicit none

  private
  public :: second_order_sgg_reconstruction_t

  type, extends(abstract_reconstruction_t) :: second_order_sgg_reconstruction_t
    !< Implementation of a 2nd order piecewise-linear reconstruction operator. The
    !< "sgg" just means that it uses the standard Green-Gauss gradient form.
  contains
    procedure, public :: initialize
    procedure, public :: reconstruct
    ! procedure, public :: reconstruct_point
    procedure, private :: estimate_gradient
    procedure, public :: copy
    final :: finalize
  end type

contains

  subroutine initialize(self, input, grid_target)
    !< Construct the second_order_sgg_reconstruction_t type

    class(second_order_sgg_reconstruction_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid_target

    integer(ik) :: alloc_status

    call debug_print('Initializing second_order_sgg_reconstruction_t', __FILE__, __LINE__)

    self%order = 2
    self%name = 'piecewise_linear_reconstruction'

    self%grid => grid_target

    call self%set_slope_limiter(name=input%slope_limiter)

  end subroutine initialize

  subroutine finalize(self)
    !< Finalize the second_order_sgg_reconstruction_t type
    type(second_order_sgg_reconstruction_t), intent(inout) :: self
    integer(ik) :: alloc_status

    call debug_print('Running second_order_sgg_reconstruction_t%finalize()', __FILE__, __LINE__)

    if(associated(self%grid)) nullify(self%grid)
    ! if(associated(self%primitive_vars)) nullify(self%primitive_vars)
    ! if(allocated(self%cell_gradient)) then
    !   deallocate(self%cell_gradient, stat=alloc_status)
    !   if(alloc_status /= 0) then
    !     error stop "Unable to deallocate second_order_sgg_reconstruction_t%cell_gradient"
    !   end if
    ! end if
  end subroutine finalize

  subroutine copy(out_recon, in_recon)
    class(abstract_reconstruction_t), intent(in) :: in_recon
    class(second_order_sgg_reconstruction_t), intent(inout) :: out_recon

    ! call debug_print('Running second_order_sgg_reconstruction_t%copy()', __FILE__, __LINE__)

    ! if(associated(out_recon%grid)) nullify(out_recon%grid)
    ! out_recon%grid => in_recon%grid

    ! if(associated(out_recon%primitive_vars)) nullify(out_recon%primitive_vars)
    ! out_recon%primitive_vars => in_recon%primitive_vars

    ! if(allocated(out_recon%name)) deallocate(out_recon%name)
    ! allocate(out_recon%name, source=in_recon%name)

    ! if(allocated(out_recon%cell_gradient)) deallocate(out_recon%cell_gradient)
    ! allocate(out_recon%cell_gradient, source=in_recon%cell_gradient)

    ! out_recon%limiter = in_recon%limiter
    ! out_recon%domain_has_been_reconstructed = .false.
  end subroutine

  subroutine reconstruct(self, primitive_var, reconstructed_var, lbounds)
    !< Reconstruct each corner/midpoint. This converts the cell centered conserved
    !< quantities [rho, rho u, rho v, e] to reconstructed primitive variables [rho, u, v, p]
    !< based on the chosen reconstruction order, e.g. using a piecewise-linear function based on the
    !< selected cell and it's neighbors. Rather than do it a point at a time, this reuses some
    !< of the data necessary, like the cell average and gradient
    class(second_order_sgg_reconstruction_t), intent(inout) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: primitive_var !< (i,j); cell primitive variable to reconstruct
    real(rk), dimension(:, lbounds(1):, lbounds(2):), contiguous, intent(out) :: reconstructed_var
    !< ((corner1:midpoint4), i, j); reconstructed variable, the first index is 1:8, or (c1,m1,c2,m2,c3,m3,c4,m4), c:corner, m:midpoint

    real(rk), dimension(:, :), allocatable :: grad_x
    real(rk), dimension(:, :), allocatable :: grad_y

    integer(ik) :: i, j, p  !< cell i,j index
    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: ilo_bc, ihi_bc, jlo_bc, jhi_bc
    real(rk) :: U_cell_ave_max
    real(rk) :: U_cell_ave_min
    real(rk) :: U_recon_max
    real(rk) :: U_recon_min
    real(rk), dimension(2) :: grad_u_limited
    real(rk) :: beta_min
    real(rk) :: beta_max
    real(rk) :: phi_lim
    real(rk) :: x, y, x_ij, y_ij

    real(rk), dimension(8) :: reconstructed_cell !< reconstructed corner/midpoints for the current cell

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

    allocate(grad_x(ilo_bc:ihi_bc, jlo_bc:jhi_bc)) ! smaller than the primitive_var b/c of ghost regions
    allocate(grad_y(ilo_bc:ihi_bc, jlo_bc:jhi_bc)) ! smaller than the primitive_var b/c of ghost regions

    call self%estimate_gradient(primitive_var=primitive_var, grad_x=grad_x, grad_y=grad_y, lbounds=lbounds)

    !  Reconstruction points for each cell (corners and mid-points)
    !  C4---M3---C3
    !  |         |
    !  M4   x    M2
    !  |         |
    !  C1---M1---C2

    !!$omp parallel default(none), &
    !!$omp firstprivate(ilo, ihi, jlo, jhi) &
    !!$omp private(i, j, x, y, x_ij, y_ij) &
    !!$omp private(U_cell_ave_max, U_cell_ave_min, U_recon_max, U_recon_min) &
    !!$omp private(beta_min, beta_max, phi_lim) &
    !!$omp shared(reconstructed_cell, reconstructed_var, self, grad_x, grad_y, primitive_var)
    !!$omp do
    do j = jlo, jhi
      do i = ilo, ihi

        x_ij = self%grid%cell_centroid_x(i, j)
        y_ij = self%grid%cell_centroid_y(i, j)

        U_cell_ave_max = maxval(primitive_var(i - 1:i + 1, j - 1:j + 1))
        U_cell_ave_min = minval(primitive_var(i - 1:i + 1, j - 1:j + 1))

        ! First, find the unlimited interpolated values for each corner and midpoint
        do p = 1, 8
          x = self%grid%cell_node_x(p, i, j)
          y = self%grid%cell_node_y(p, i, j)
          reconstructed_cell(p) = primitive_var(i, j) + grad_x(i, j) * (x - x_ij) + &
                                  grad_y(i, j) * (y - y_ij)
        end do

        U_recon_max = maxval(reconstructed_cell)
        U_recon_min = minval(reconstructed_cell)

        associate(U_ave=>primitive_var(i, j))
          beta_min = 0.0_rk
          beta_max = 0.0_rk
          if(abs(U_recon_min - U_ave) > 0.0_rk) then
            beta_min = max(0.0_rk,(U_cell_ave_min - U_ave) / (U_recon_min - U_ave))
          end if

          if(abs(U_recon_max - U_ave) > 0.0_rk) then
            beta_max = max(0.0_rk,(U_cell_ave_max - U_ave) / (U_recon_max - U_ave))
          end if
          phi_lim = min(1.0_rk, beta_min, beta_max)
        end associate

        do p = 1, 8
          x = self%grid%cell_node_x(p, i, j)
          y = self%grid%cell_node_y(p, i, j)
          reconstructed_cell(p) = primitive_var(i, j) + phi_lim * grad_x(i, j) * (x - x_ij) + &
                                  phi_lim * grad_y(i, j) * (y - y_ij)
        end do

        reconstructed_var(:, i, j) = reconstructed_cell
      end do
    end do
    !!$omp end do
    !!$omp end parallel

    self%domain_has_been_reconstructed = .true.
  end subroutine reconstruct

  subroutine estimate_gradient(self, primitive_var, grad_x, grad_y, lbounds)
    !< Estimate the slope-limited gradient of the primitive variables in the cell (i,j). This assumes
    !< a quadrilateral structured grid
    class(second_order_sgg_reconstruction_t), intent(in) :: self
    integer(ik), dimension(2), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(in) :: primitive_var !< (i,j); data to estimate the gradient of
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(out) :: grad_x !< (i,j); data to estimate the gradient of
    real(rk), dimension(lbounds(1):, lbounds(2):), contiguous, intent(out) :: grad_y !< (i,j); data to estimate the gradient of

    integer(ik) :: i, j, k
    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: ilo_bc, ihi_bc, jlo_bc, jhi_bc

    real(rk), dimension(4) :: edge_lengths  !< length of each face
    real(rk), dimension(4) :: v_edge  !< length of each face
    real(rk), dimension(4) :: n_x  !< normal vectors of each face
    real(rk), dimension(4) :: n_y  !< normal vectors of each face
    real(rk) :: d_dx, d_dy

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
    !$omp private(n_x, n_y, v_edge, edge_lengths) &
    !$omp shared(grad_x, grad_y, primitive_var, self) &
    !$omp reduction(+:d_dx) &
    !$omp reduction(+:d_dy)
    !$omp do
    do j = jlo, jhi
      do i = ilo, ihi

        associate(center=>primitive_var(i, j), &      ! current cell
                  bottom=>primitive_var(i, j - 1), &  ! bottom cell
                  right=>primitive_var(i + 1, j), &   ! right cell
                  top=>primitive_var(i, j + 1), &     ! top cell
                  left=>primitive_var(i - 1, j), &      ! left cell
                  vol_center=>self%grid%cell_volume(i, j), &      ! current cell
                  vol_bottom=>self%grid%cell_volume(i, j - 1), &  ! bottom cell
                  vol_right=>self%grid%cell_volume(i + 1, j), &   ! right cell
                  vol_top=>self%grid%cell_volume(i, j + 1), &     ! top cell
                  vol_left=>self%grid%cell_volume(i - 1, j))

          ! Edge (face) interface data
          edge_lengths(1) = self%grid%cell_edge_lengths(1, i, j - 1)  ! bottom
          edge_lengths(2) = self%grid%cell_edge_lengths(2, i + 1, j)  ! right
          edge_lengths(3) = self%grid%cell_edge_lengths(3, i, j + 1)  ! top
          edge_lengths(4) = self%grid%cell_edge_lengths(4, i - 1, j)  ! left

          n_x(1) = self%grid%cell_edge_norm_vectors(1, 1, i, j - 1)  ! bottom
          n_x(2) = self%grid%cell_edge_norm_vectors(1, 2, i + 1, j)  ! right
          n_x(3) = self%grid%cell_edge_norm_vectors(1, 3, i, j + 1)  ! top
          n_x(4) = self%grid%cell_edge_norm_vectors(1, 4, i - 1, j)  ! left
          n_y(1) = self%grid%cell_edge_norm_vectors(2, 1, i, j - 1)  ! bottom
          n_y(2) = self%grid%cell_edge_norm_vectors(2, 2, i + 1, j)  ! right
          n_y(3) = self%grid%cell_edge_norm_vectors(2, 3, i, j + 1)  ! top
          n_y(4) = self%grid%cell_edge_norm_vectors(2, 4, i - 1, j)  ! left

          v_edge(1) = (center * vol_center + bottom * vol_bottom) / (vol_center + vol_bottom)  ! bottom
          v_edge(2) = (center * vol_center + right * vol_right) / (vol_center + vol_right)  ! right
          v_edge(3) = (center * vol_center + top * vol_top) / (vol_center + vol_top)  ! top
          v_edge(4) = (center * vol_center + left * vol_left) / (vol_center + vol_left)  ! left

          d_dx = sum(v_edge * n_x * edge_lengths)
          d_dy = sum(v_edge * n_y * edge_lengths)

        end associate
        grad_x(i, j) = d_dx / self%grid%cell_volume(i, j)
        grad_y(i, j) = d_dy / self%grid%cell_volume(i, j)
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine estimate_gradient

end module mod_second_order_sgg_reconstruction
