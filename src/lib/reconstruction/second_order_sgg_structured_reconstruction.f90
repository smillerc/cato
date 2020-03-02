module mod_second_order_sgg_structured_reconstruction

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_grid, only: grid_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use mod_floating_point_utils, only: equal
  use mod_gradients, only: green_gauss_gradient_limited, get_smoothness

  implicit none

  private
  public :: second_order_sgg_structured_reconstruction_t

  type, extends(abstract_reconstruction_t) :: second_order_sgg_structured_reconstruction_t
    !< Implementation of a 2nd order piecewise-linear reconstruction operator described
    !< in the original FVLEG paper
  contains
    procedure, public :: initialize
    procedure, public :: reconstruct_domain
    procedure, public :: reconstruct_point
    procedure, private :: estimate_gradients
    procedure, public :: copy
    final :: finalize
  end type

contains

  subroutine initialize(self, input, grid_target)
    !< Construct the second_order_sgg_structured_reconstruction_t type

    class(second_order_sgg_structured_reconstruction_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid_target

    integer(ik) :: alloc_status

    call debug_print('Initializing second_order_sgg_structured_reconstruction_t', __FILE__, __LINE__)

    self%order = 2
    self%name = 'piecewise_linear_reconstruction'

    self%grid => grid_target

    call self%set_slope_limiter(name=input%slope_limiter)

    associate(imin=>grid_target%ilo_bc_cell, imax=>grid_target%ihi_bc_cell, &
              jmin=>grid_target%jlo_bc_cell, jmax=>grid_target%jhi_bc_cell)

      allocate(self%cell_gradient(4, 2, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) then
        error stop "Unable to allocate second_order_sgg_structured_reconstruction_t%cell_gradient"
      end if
      self%cell_gradient = 0.0_rk
    end associate

  end subroutine initialize

  subroutine finalize(self)
    !< Finalize the second_order_sgg_structured_reconstruction_t type
    type(second_order_sgg_structured_reconstruction_t), intent(inout) :: self
    integer(ik) :: alloc_status

    call debug_print('Running second_order_sgg_structured_reconstruction_t%finalize()', __FILE__, __LINE__)

    if(associated(self%grid)) nullify(self%grid)
    if(associated(self%primitive_vars)) nullify(self%primitive_vars)
    if(allocated(self%cell_gradient)) then
      deallocate(self%cell_gradient, stat=alloc_status)
      if(alloc_status /= 0) then
        error stop "Unable to deallocate second_order_sgg_structured_reconstruction_t%cell_gradient"
      end if
    end if
  end subroutine finalize

  subroutine copy(out_recon, in_recon)
    class(abstract_reconstruction_t), intent(in) :: in_recon
    class(second_order_sgg_structured_reconstruction_t), intent(inout) :: out_recon

    call debug_print('Running second_order_sgg_structured_reconstruction_t%copy()', __FILE__, __LINE__)

    if(associated(out_recon%grid)) nullify(out_recon%grid)
    out_recon%grid => in_recon%grid

    if(associated(out_recon%primitive_vars)) nullify(out_recon%primitive_vars)
    out_recon%primitive_vars => in_recon%primitive_vars

    if(allocated(out_recon%name)) deallocate(out_recon%name)
    allocate(out_recon%name, source=in_recon%name)

    if(allocated(out_recon%cell_gradient)) deallocate(out_recon%cell_gradient)
    allocate(out_recon%cell_gradient, source=in_recon%cell_gradient)

    out_recon%limiter = in_recon%limiter
    out_recon%domain_has_been_reconstructed = .false.
  end subroutine

  function reconstruct_point(self, xy, cell_ij) result(V_bar)
    !< Reconstruct the value of the primitive variables (U) at location (x,y)
    !< withing a cell (i,j)

    class(second_order_sgg_structured_reconstruction_t), intent(in) :: self
    real(rk), dimension(2), intent(in) :: xy !< where should V_bar be reconstructed at?
    integer(ik), dimension(2), intent(in) :: cell_ij !< cell (i,j) indices to reconstruct within
    real(rk), dimension(4) :: V_bar  !< V_bar = reconstructed [rho, u, v, p]
    ! real(rk), dimension(2) :: centroid_xy !< (x,y) location of the cell centroid
    integer(ik) :: i, j

    i = cell_ij(1); j = cell_ij(2)
    ! centroid_xy = self%grid%get_cell_centroid_xy(i, j)
    V_bar = self%interpolate(i=i, j=j, x=xy(1), y=xy(2))
    ! if (i == 601 .or. i == 602) then
    ! print*, 'reconstructing point, ij:', cell_ij
    ! print*, 'reconstructing point, xy:', xy
    ! print*, 'grad_u(:,1)', self%cell_gradient(:,1,i,j)
    ! print*, 'grad_u(:,2)', self%cell_gradient(:,2,i,j)
    ! print*, 'cell_ave:', self%primitive_vars(:,i,j)
    ! print*, 'V_bar:', V_bar
    ! end if
  end function reconstruct_point

  subroutine reconstruct_domain(self, reconstructed_domain, lbounds)
    !< Reconstruct each corner/midpoint. This converts the cell centered conserved
    !< quantities [rho, rho u, rho v, e] to reconstructed primitive variables [rho, u, v, p]
    !< based on the chosen reconstruction order, e.g. using a piecewise-linear function based on the
    !< selected cell and it's neighbors. Rather than do it a point at a time, this reuses some
    !< of the data necessary, like the cell average and gradient
    class(second_order_sgg_structured_reconstruction_t), intent(inout) :: self
    integer(ik), dimension(5), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):, &
                        lbounds(4):, lbounds(5):), intent(out) :: reconstructed_domain

    !< ((rho, u ,v, p), point, node/midpoint, i, j);
    !< The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints

    integer(ik) :: i, j, l  !< cell i,j index
    integer(ik) :: m  !< midpoint index
    integer(ik) :: c  !< corner index
    integer(ik) :: n, p
    integer(ik) :: phi, nhi, ilo, ihi, jlo, jhi
    real(rk), dimension(4) :: U_cell_ave_max
    real(rk), dimension(4) :: U_cell_ave_min
    real(rk), dimension(4) :: U_recon_max
    real(rk), dimension(4) :: U_recon_min
    real(rk), dimension(4, 2) :: grad_u_limited
    real(rk), dimension(4) :: beta_min
    real(rk), dimension(4) :: beta_max
    real(rk), dimension(4) :: phi_lim

    real(rk), dimension(4, 4, 2) :: reconstructed_cell !< reconstructed corner/midpoints for the current cell

    ! Bounds do not include ghost cells. Ghost cells get their
    ! reconstructed values and gradients from the boundary conditions
    ilo = lbound(reconstructed_domain, dim=4) + 1
    ihi = ubound(reconstructed_domain, dim=4) - 1
    jlo = lbound(reconstructed_domain, dim=5) + 1
    jhi = ubound(reconstructed_domain, dim=5) - 1

    ! Find the unlimited gradients sets the self%cell_gradients array
    call self%estimate_gradients()
    do j = jlo, jhi
      do i = ilo, ihi
        ! Now reinterpolate with the limited gradient
        do n = 1, 2 ! First do corners, then to midpoints
          do p = 1, 4 ! Loop through each point (N1-N4, and M1-M4)
            associate(x=>self%grid%cell_node_xy(1, p, n, i, j), &
                      y=>self%grid%cell_node_xy(2, p, n, i, j))
              reconstructed_cell(:, p, n) = self%interpolate(i, j, x, y)
            end associate
          end do
        end do
        reconstructed_domain(:, :, :, i, j) = reconstructed_cell
      end do
    end do

    self%domain_has_been_reconstructed = .true.
  end subroutine reconstruct_domain

  subroutine estimate_gradients(self)
    !< Estimate the slope-limited gradient of the primitive variables in the cell (i,j). This assumes
    !< a quadrilateral structured grid
    class(second_order_sgg_structured_reconstruction_t), intent(inout) :: self
    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi

    real(rk), dimension(4, 5) :: prim_vars  !< primitive variables of current and neighbor cells
    real(rk), dimension(5) :: volumes  !< volume of each cell
    real(rk), dimension(4) :: edge_lengths  !< length of each face
    real(rk), dimension(2, 4) :: edge_normals  !< normal vectors of each face

    real(rk), dimension(4, 2) :: phi_lim !< slope limiter
    real(rk), dimension(4, 2) :: gradient
    real(rk), dimension(4) :: smooth_updown !< slope limiter scalar function
    real(rk), dimension(4) :: smooth_leftright !< slope limiter scalar function

    phi_lim = 0.0_rk

    ilo = lbound(self%primitive_vars, dim=2) + 1
    ihi = ubound(self%primitive_vars, dim=2) - 1
    jlo = lbound(self%primitive_vars, dim=3) + 1
    jhi = ubound(self%primitive_vars, dim=3) - 1

    do j = jlo, jhi
      do i = ilo, ihi
        ! current cell and neighbor cell [bottom, right, top, left] information for gradient estimation
        prim_vars(:, 1) = self%primitive_vars(:, i, j)      ! current
        prim_vars(:, 2) = self%primitive_vars(:, i, j - 1)  ! bottom
        prim_vars(:, 3) = self%primitive_vars(:, i + 1, j)  ! right
        prim_vars(:, 4) = self%primitive_vars(:, i, j + 1)  ! top
        prim_vars(:, 5) = self%primitive_vars(:, i - 1, j)  ! left

        volumes(1) = self%grid%cell_volume(i, j)      ! current
        volumes(2) = self%grid%cell_volume(i, j - 1)  ! bottom
        volumes(3) = self%grid%cell_volume(i + 1, j)  ! right
        volumes(4) = self%grid%cell_volume(i, j + 1)  ! top
        volumes(5) = self%grid%cell_volume(i - 1, j)  ! left

        ! Edge (face) interface data
        edge_lengths(1) = self%grid%cell_edge_lengths(1, i, j - 1)  ! bottom
        edge_lengths(2) = self%grid%cell_edge_lengths(2, i + 1, j)  ! right
        edge_lengths(3) = self%grid%cell_edge_lengths(3, i, j + 1)  ! top
        edge_lengths(4) = self%grid%cell_edge_lengths(4, i - 1, j)  ! left

        edge_normals(:, 1) = self%grid%cell_edge_norm_vectors(:, 1, i, j - 1)  ! bottom
        edge_normals(:, 2) = self%grid%cell_edge_norm_vectors(:, 2, i + 1, j)  ! right
        edge_normals(:, 3) = self%grid%cell_edge_norm_vectors(:, 3, i, j + 1)  ! top
        edge_normals(:, 4) = self%grid%cell_edge_norm_vectors(:, 4, i - 1, j)  ! left

        self%cell_gradient(:, :, i, j) = green_gauss_gradient_limited(prim_vars, volumes, edge_lengths, edge_normals)
      end do
    end do
  end subroutine estimate_gradients

end module mod_second_order_sgg_structured_reconstruction
