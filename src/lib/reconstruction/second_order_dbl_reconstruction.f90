
module mod_second_order_dbl_reconstruction
  !< Summary: Reconstruct using the discontinuous bilinear recovery (DBL)
  !< References:
  !< [1] Lukacova-Medvid'ová, M., Saibertová, J., and Warnecke, G, "Finite volume evolution galerkin methods for nonlinear hyperbolic systems",
  !<     https://doi.org/10.1006/jcph.2002.7207

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_grid, only: grid_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t
  use mod_eos, only: eos
  use mod_floating_point_utils, only: equal
  use mod_gradients, only: green_gauss_gradient, get_smoothness

  implicit none

  private
  public :: second_order_dbl_reconstruction_t

  type, extends(abstract_reconstruction_t) :: second_order_dbl_reconstruction_t
    !< Implementation of a 2nd order piecewise-linear reconstruction operator. The
    !< "sgg" just means that it uses the standard Green-Gauss gradient form.
  contains
    procedure, public :: initialize
    procedure, public :: reconstruct_domain
    procedure, public :: reconstruct_point
    procedure, private :: find_extrema
    procedure, private, nopass :: limit
    procedure, public :: copy
    final :: finalize
  end type

contains

  subroutine initialize(self, input, grid_target)
    !< Construct the second_order_dbl_reconstruction_t type

    class(second_order_dbl_reconstruction_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid_target

    integer(ik) :: alloc_status

    call debug_print('Initializing second_order_dbl_reconstruction_t', __FILE__, __LINE__)

    self%order = 2
    self%name = 'piecewise_linear_reconstruction'

    self%grid => grid_target

    call self%set_slope_limiter(name=input%slope_limiter)

    associate(imin=>grid_target%ilo_bc_cell, imax=>grid_target%ihi_bc_cell, &
              jmin=>grid_target%jlo_bc_cell, jmax=>grid_target%jhi_bc_cell)

      allocate(self%cell_gradient(4, 2, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) then
        error stop "Unable to allocate second_order_dbl_reconstruction_t%cell_gradient"
      end if
      self%cell_gradient = 0.0_rk
    end associate

  end subroutine initialize

  subroutine finalize(self)
    !< Finalize the second_order_dbl_reconstruction_t type
    type(second_order_dbl_reconstruction_t), intent(inout) :: self
    integer(ik) :: alloc_status

    call debug_print('Running second_order_dbl_reconstruction_t%finalize()', __FILE__, __LINE__)

    if(associated(self%grid)) nullify(self%grid)
    if(associated(self%primitive_vars)) nullify(self%primitive_vars)
    if(allocated(self%cell_gradient)) then
      deallocate(self%cell_gradient, stat=alloc_status)
      if(alloc_status /= 0) then
        error stop "Unable to deallocate second_order_dbl_reconstruction_t%cell_gradient"
      end if
    end if
  end subroutine finalize

  subroutine copy(out_recon, in_recon)
    class(abstract_reconstruction_t), intent(in) :: in_recon
    class(second_order_dbl_reconstruction_t), intent(inout) :: out_recon

    call debug_print('Running second_order_dbl_reconstruction_t%copy()', __FILE__, __LINE__)

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

  pure function reconstruct_point(self, xy, cell_ij) result(V_bar)
    !< Reconstruct the value of the primitive variables (U) at location (x,y)
    !< withing a cell (i,j)

    class(second_order_dbl_reconstruction_t), intent(in) :: self
    real(rk), dimension(2), intent(in) :: xy !< where should V_bar be reconstructed at?
    real(rk), dimension(4) :: V_bar  !< V_bar = reconstructed [rho, u, v, p]
    integer(ik), dimension(2), intent(in) :: cell_ij !< cell (i,j) indices to reconstruct within
    real(rk), dimension(2) :: centroid_xy !< (x,y) location of the cell centroid
    integer(ik) :: i, j

    i = cell_ij(1); j = cell_ij(2)
    centroid_xy = self%grid%get_cell_centroid_xy(i, j)

    if(.not. self%domain_has_been_reconstructed) then
      error stop "Error in second_order_dbl_reconstruction_t%reconstruct_point(), "// &
        "domain_has_been_reconstructed is false, but should be true"
    end if

    associate(dU_dx=>self%cell_gradient(:, 1, i, j), &
              dU_dy=>self%cell_gradient(:, 2, i, j), &
              cell_ave=>self%primitive_vars(:, i, j), &
              x=>xy(1), y=>xy(2), &
              x_ij=>centroid_xy(1), y_ij=>centroid_xy(2))

      V_bar = cell_ave + dU_dx * (x - x_ij) + dU_dy * (y - y_ij)
    end associate

  end function reconstruct_point

  subroutine reconstruct_domain(self, reconstructed_domain, lbounds)
    !< Reconstruct each corner/midpoint. This converts the cell centered conserved
    !< quantities [rho, rho u, rho v, e] to reconstructed primitive variables [rho, u, v, p]
    !< based on the chosen reconstruction order, e.g. using a piecewise-linear function based on the
    !< selected cell and it's neighbors. Rather than do it a point at a time, this reuses some
    !< of the data necessary, like the cell average and gradient
    class(second_order_dbl_reconstruction_t), intent(inout) :: self
    integer(ik), dimension(5), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):, &
                        lbounds(4):, lbounds(5):), intent(out) :: reconstructed_domain

    !< ((rho, u ,v, p), point, node/midpoint, i, j);
    !< The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints

    real(rk), dimension(2) :: centroid_xy !< (x,y) location of the cell centroid
    integer(ik) :: i, j, n, p  !< cell i,j index
    integer(ik) :: ihi, ilo, jhi, jlo
    integer(ik), parameter :: c = 1  !< corner index
    integer(ik), parameter :: m = 2  !< midpoint index
    real(rk), dimension(4, 4, 2) :: U_max
    real(rk), dimension(4, 4, 2) :: U_min
    real(rk), dimension(4) :: U_tilde
    real(rk), dimension(4) :: delta_x
    real(rk), dimension(4) :: delta_y
    real(rk), dimension(4) :: mu_x
    real(rk), dimension(4) :: mu_y
    real(rk), dimension(4, 4) :: face_prim_vars

    real(rk), dimension(2, 4, 2) :: node_pos
    real(rk), dimension(2, 4, 2) :: h !< cell size parameter
    real(rk), dimension(4) :: phi_lim !< slope limiter for (rho, u, v, p)

    ! Bounds do not include ghost cells. Ghost cells get their
    ! reconstructed values and gradients from the boundary conditions
    ilo = lbound(reconstructed_domain, dim=4) + 1
    ihi = ubound(reconstructed_domain, dim=4) - 1
    jlo = lbound(reconstructed_domain, dim=5) + 1
    jhi = ubound(reconstructed_domain, dim=5) - 1

    !  Reconstruction points for each cell (corners and mid-points)
    !          E3  (edge 3)
    !     C4---M3---C3
    !     |         |
    ! E4  M4   x    M2  E2
    !     |         |
    !     C1---M1---C2
    !          E1
    do j = jlo, jhi
      do i = ilo, ihi
        centroid_xy = self%grid%get_cell_centroid_xy(i=i, j=j)

        ! Volume average to get the interface values
        associate(v=>self%grid%cell_volume, U=>self%primitive_vars, U_f=>face_prim_vars)
          U_f(:, 1) = (U(:, i, j) * v(i, j) + U(:, i, j - 1) * v(i, j - 1)) / &
                      (v(i, j) + v(i, j - 1))  ! down

          U_f(:, 2) = (U(:, i, j) * v(i, j) + U(:, i + 1, j) * v(i + 1, j)) / &
                      (v(i, j) + v(i + 1, j))  ! right

          U_f(:, 3) = (U(:, i, j) * v(i, j) + U(:, i, j + 1) * v(i, j + 1)) / &
                      (v(i, j) + v(i, j + 1))  ! up

          U_f(:, 4) = (U(:, i, j) * v(i, j) + U(:, i - 1, j) * v(i - 1, j)) / &
                      (v(i, j) + v(i - 1, j))  ! left

          mu_x = 0.5_rk * (U_f(:, 2) + U_f(:, 4)) ! Average of face values
          mu_y = 0.5_rk * (U_f(:, 3) + U_f(:, 1)) ! Average of face values

          delta_x = U_f(:, 2) - U_f(:, 4) ! Difference in face values
          delta_y = U_f(:, 3) - U_f(:, 1) ! Difference in face values
        end associate

        node_pos = self%grid%cell_node_xy(:, :, :, i, j) !< (x,y), (p1,p2,p3,p4), (node/midpoint),

        ! Calculate the edge lengths (needed for non-uniform grid)
        associate(x_c1=>node_pos(1, 1, 1), y_c1=>node_pos(2, 1, 1), &
                  x_m1=>node_pos(1, 1, 2), y_m1=>node_pos(2, 1, 2), &
                  x_c2=>node_pos(1, 2, 1), y_c2=>node_pos(2, 2, 1), &
                  x_m2=>node_pos(1, 2, 2), y_m2=>node_pos(2, 2, 2), &
                  x_c3=>node_pos(1, 3, 1), y_c3=>node_pos(2, 3, 1), &
                  x_m3=>node_pos(1, 3, 2), y_m3=>node_pos(2, 3, 2), &
                  x_c4=>node_pos(1, 4, 1), y_c4=>node_pos(2, 4, 1), &
                  x_m4=>node_pos(1, 4, 2), y_m4=>node_pos(2, 4, 2))

          h(1, 1, c) = abs(x_c2 - x_c1) ! C1 x-component
          h(2, 1, c) = abs(y_c4 - y_c1) ! C1 y-component

          h(1, 1, m) = abs(x_c2 - x_c1) ! M1 x-component
          h(2, 1, m) = abs(y_m3 - y_m1) ! M1 y-component

          h(1, 2, c) = abs(x_c2 - x_c1) ! C2 x-component
          h(2, 2, c) = abs(y_c3 - y_c2) ! C2 y-component

          h(1, 2, m) = abs(x_m2 - x_m4) ! M2 x-component
          h(2, 2, m) = abs(y_c3 - y_c2) ! M2 y-component

          h(1, 3, c) = abs(x_c3 - x_c4) ! C3 x-component
          h(2, 3, c) = abs(y_c3 - y_c2) ! C3 y-component

          h(1, 3, m) = abs(x_c3 - x_c4) ! M3 x-component
          h(2, 3, m) = abs(y_m3 - y_m1) ! M3 y-component

          h(1, 4, c) = abs(x_c3 - x_c4) ! C4 x-component
          h(2, 4, c) = abs(y_c4 - y_c1) ! C4 y-component

          h(1, 4, m) = abs(x_m2 - x_m4) ! M4 x-component
          h(2, 4, m) = abs(y_c4 - y_c1) ! M4 y-component
        end associate

        call self%find_extrema(i, j, U_max, U_min)

        do n = 1, 2 ! corner / midpoint
          do p = 1, 4 ! point

            associate(x=>node_pos(1, p, n), y=>node_pos(2, p, n), h_x=>h(1, p, n), h_y=>h(2, p, n), &
                      U=>self%primitive_vars, x_ij=>centroid_xy(1), y_ij=>centroid_xy(2))

              U_tilde = ((x - x_ij) / h_x) * mu_x * mu_y * mu_y * delta_x + &
                        ((y - y_ij) / h_y) * mu_x * mu_x * mu_y * delta_y + &
                        (((x - x_ij) * (y - y_ij)) / (h_x * h_y)) * mu_x * mu_y * delta_x * delta_y

              phi_lim = self%limit(u_ij=U(:, i, j), u_tilde=U_tilde, u_max=U_max(:, p, n), u_min=U_min(:, p, n))
              reconstructed_domain(:, p, n, i, j) = U(:, i, j) + phi_lim * U_tilde
            end associate
          end do
        end do
      end do
    end do

    self%domain_has_been_reconstructed = .true.
  end subroutine reconstruct_domain

  real(rk) elemental function limit(U_ij, U_tilde, U_max, U_min) result(phi_lim)
    !< Apply a limiter. Based on Eq 3.9 in Ref [1]
    real(rk), intent(in) :: U_ij
    real(rk), intent(in) :: U_tilde
    real(rk), intent(in) :: U_max
    real(rk), intent(in) :: U_min

    if(equal(U_tilde, U_ij)) then
      phi_lim = 1.0_rk
    else if(U_tilde > U_ij) then
      phi_lim = min(1.0_rk,(U_max - U_ij) / (U_tilde - U_ij))
    else if(U_tilde < U_ij) then
      phi_lim = min(1.0_rk,(U_min - U_ij) / (U_tilde - U_ij))
    end if

  end function limit

  subroutine find_extrema(self, i, j, U_max, U_min)
    !< At each corner and minpoint, find the min/max
    class(second_order_dbl_reconstruction_t), intent(in) :: self
    real(rk), dimension(4, 4, 2), intent(out) :: U_max !< ((rho, u, v, p), (point 1 - 4), (corner=1/midpoint=2))
    real(rk), dimension(4, 4, 2), intent(out) :: U_min !< ((rho, u, v, p), (point 1 - 4), (corner=1/midpoint=2))
    integer(ik), intent(in) :: i, j
    integer(ik) :: l
    integer(ik), parameter :: c = 1 !< corner index
    integer(ik), parameter :: m = 2 !< midpoint index

    ! Find the extrema at each node point
    associate(U=>self%primitive_vars)

      ! C1
      do l = 1, 4
        U_max(1, l, c) = max(U(l, i, j), U(l, i - 1, j), U(l, i - 1, j - 1), U(l, i, j - 1))
        U_min(1, l, c) = min(U(l, i, j), U(l, i - 1, j), U(l, i - 1, j - 1), U(l, i, j - 1))
      end do

      ! C2
      do l = 1, 4
        U_max(2, l, c) = max(U(l, i, j), U(l, i + 1, j), U(l, i + 1, j - 1), U(l, i, j - 1))
        U_min(2, l, c) = min(U(l, i, j), U(l, i + 1, j), U(l, i + 1, j - 1), U(l, i, j - 1))
      end do

      ! C3
      do l = 1, 4
        U_max(3, l, c) = max(U(l, i, j), U(l, i + 1, j), U(l, i + 1, j + 1), U(l, i, j + 1))
        U_min(3, l, c) = min(U(l, i, j), U(l, i + 1, j), U(l, i + 1, j + 1), U(l, i, j + 1))
      end do

      ! C4
      do l = 1, 4
        U_max(4, l, c) = max(U(l, i, j), U(l, i, j + 1), U(l, i - 1, j + 1), U(l, i - 1, j))
        U_min(4, l, c) = min(U(l, i, j), U(l, i, j + 1), U(l, i - 1, j + 1), U(l, i - 1, j))
      end do

      ! M1
      do l = 1, 4
        U_max(1, l, m) = max(U(l, i, j), U(l, i, j - 1))
        U_min(1, l, m) = min(U(l, i, j), U(l, i, j - 1))
      end do

      ! M2
      do l = 1, 4
        U_max(2, l, m) = max(U(l, i, j), U(l, i + 1, j))
        U_min(2, l, m) = min(U(l, i, j), U(l, i + 1, j))
      end do

      ! M3
      do l = 1, 4
        U_max(3, l, m) = max(U(l, i, j), U(l, i, j + 1))
        U_min(3, l, m) = min(U(l, i, j), U(l, i, j + 1))
      end do

      ! M4
      do l = 1, 4
        U_max(4, l, m) = max(U(l, i, j), U(l, i - 1, j))
        U_min(4, l, m) = min(U(l, i, j), U(l, i - 1, j))
      end do

    end associate

  end subroutine find_extrema

end module mod_second_order_dbl_reconstruction
