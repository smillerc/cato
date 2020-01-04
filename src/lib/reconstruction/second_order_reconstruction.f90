#ifdef __DEBUG__
#define debug_write write
#else
#define debug_write ! write
#endif

module mod_second_order_reconstruction

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_grid, only: grid_t
  use mod_slope_limiter, only: slope_limiter_t
  use mod_input, only: input_t

  implicit none

  private
  public :: second_order_reconstruction_t

  type, extends(abstract_reconstruction_t) :: second_order_reconstruction_t
    !< Implementation of a 2nd order piecewise-linear reconstruction operator
  contains
    procedure, public :: initialize => init_second_order
    procedure, public :: reconstruct_domain
    procedure, public :: reconstruct_point
    procedure, private :: estimate_gradients
    procedure, private :: estimate_single_gradient
    procedure, public :: copy
    final :: finalize
  end type

contains

  subroutine init_second_order(self, input, grid_target)
    !< Construct the second_order_reconstruction_t type

    class(second_order_reconstruction_t), intent(inout) :: self
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid_target
    ! integer(ik), dimension(3), intent(in) :: lbounds
    ! real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), &
    !   intent(in), target :: conserved_vars_target

    integer(ik) :: alloc_status

    call debug_print('Initializing second_order_reconstruction_t', __FILE__, __LINE__)

    self%order = 2
    self%name = 'piecewise_linear_reconstruction'

    ! self%conserved_vars => conserved_vars_target
    self%grid => grid_target

    call self%set_slope_limiter(name=input%slope_limiter)

    associate(imin=>grid_target%ilo_bc_cell, imax=>grid_target%ihi_bc_cell, &
              jmin=>grid_target%jlo_bc_cell, jmax=>grid_target%jhi_bc_cell)

      allocate(self%cell_gradient(2, 4, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) then
        error stop "Unable to allocate second_order_reconstruction_t%cell_gradient"
      end if
      self%cell_gradient = 0.0_rk
    end associate

  end subroutine init_second_order

  subroutine finalize(self)
    !< Finalize the second_order_reconstruction_t type
    type(second_order_reconstruction_t), intent(inout) :: self
    integer(ik) :: alloc_status

    call debug_print('Calling second_order_reconstruction_t%finalize()', __FILE__, __LINE__)

    if(associated(self%grid)) nullify(self%grid)
    if(associated(self%conserved_vars)) nullify(self%conserved_vars)
    if(allocated(self%cell_gradient)) then
      deallocate(self%cell_gradient, stat=alloc_status)
      if(alloc_status /= 0) then
        error stop "Unable to deallocate second_order_reconstruction_t%cell_gradient"
      end if
    end if
  end subroutine finalize

  subroutine copy(out_recon, in_recon)
    class(abstract_reconstruction_t), intent(in) :: in_recon
    class(second_order_reconstruction_t), intent(inout) :: out_recon

    call debug_print('Calling second_order_reconstruction_t%copy()', __FILE__, __LINE__)

    if(associated(out_recon%grid)) nullify(out_recon%grid)
    ! allocate(out_recon%grid, source=in_recon%grid)
    out_recon%grid => in_recon%grid

    if(associated(out_recon%conserved_vars)) nullify(out_recon%conserved_vars)
    ! allocate(out_recon%conserved_vars, source=in_recon%conserved_vars)
    out_recon%conserved_vars => in_recon%conserved_vars

    if(allocated(out_recon%name)) deallocate(out_recon%name)
    allocate(out_recon%name, source=in_recon%name)

    if(allocated(out_recon%cell_gradient)) deallocate(out_recon%cell_gradient)
    allocate(out_recon%cell_gradient, source=in_recon%cell_gradient)

    out_recon%limiter = in_recon%limiter
  end subroutine

  pure function reconstruct_point(self, xy, cell_ij) result(U_bar)
    !< Reconstruct the value of the conserved variables (U) at location (x,y)
    !< withing a cell (i,j)

    class(second_order_reconstruction_t), intent(in) :: self
    ! real(rk), dimension(:, 0:, 0:), intent(in) :: conserved_vars
    real(rk), dimension(2), intent(in) :: xy !< where should U_bar be reconstructed at?
    real(rk), dimension(4) :: U_bar  !< U_bar = reconstructed [rho, u, v, p]
    integer(ik), dimension(2), intent(in) :: cell_ij !< cell (i,j) indices to reconstruct within
    real(rk), dimension(2) :: centroid_xy !< (x,y) location of the cell centroid
    integer(ik) :: i, j

    i = cell_ij(1); j = cell_ij(2)
    centroid_xy = self%grid%get_cell_centroid_xy(i, j)

    associate(dU_dx=>self%cell_gradient(1, :, i, j), &
              dU_dy=>self%cell_gradient(2, :, i, j), &
              cell_ave=>self%conserved_vars(:, i, j), &
              x=>xy(1), y=>xy(2), &
              x_ij=>centroid_xy(1), y_ij=>centroid_xy(2))

      U_bar = cell_ave + dU_dx * (x - x_ij) + dU_dy * (y - y_ij)
    end associate

  end function reconstruct_point

  pure subroutine reconstruct_domain(self, reconstructed_domain, lbounds)
    !< Reconstruct the entire domain. Rather than do it a point at a time, this reuses some
    !< of the data necessary, like the cell average and gradient
    class(second_order_reconstruction_t), intent(inout) :: self
    ! real(rk), dimension(:, 0:, 0:), intent(in) :: conserved_vars
    integer(ik), dimension(5), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):, &
                        lbounds(4):, lbounds(5):), intent(out) :: reconstructed_domain

    !< ((rho, u ,v, p), point, node/midpoint, i, j);
    !< The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints

    real(rk), dimension(2) :: centroid_xy !< (x,y) location of the cell centroid
    integer(ik) :: i, j  !< cell i,j index
    integer(ik) :: n  !< node index -> i.e. is it a corner (1), or midpoint (2)
    integer(ik) :: p  !< point index -> which point in the grid element (corner 1-4, or midpoint 1-4)
    integer(ik) :: phi, nhi, ilo, ihi, jlo, jhi

    ! Bounds do not include ghost cells. Ghost cells get their
    ! reconstructed values and gradients from the boundary conditions
    phi = ubound(reconstructed_domain, dim=2)
    nhi = ubound(reconstructed_domain, dim=3)
    ilo = lbound(reconstructed_domain, dim=4) + 1
    ihi = ubound(reconstructed_domain, dim=4) - 1
    jlo = lbound(reconstructed_domain, dim=5) + 1
    jhi = ubound(reconstructed_domain, dim=5) - 1

#ifdef __DEBUG__
    do j = jlo, jhi
      do i = ilo, ihi
#else
        do concurrent(j=jlo:jhi)
          do concurrent(i=ilo:ihi)
#endif

            self%cell_gradient(:, :, i, j) = self%estimate_gradients(i, j)
            centroid_xy = self%grid%get_cell_centroid_xy(i=i, j=j)

            ! First do corners, then to midpoints
            do concurrent(n=1:nhi)

              ! Loop through each point (N1-N4, and M1-M4)
              do concurrent(p=1:phi)

                associate(U_bar=>reconstructed_domain, &
                          cell_ave=>self%conserved_vars(:, i, j), &
                          x=>self%grid%cell_node_xy(1, p, n, i, j), &
                          y=>self%grid%cell_node_xy(2, p, n, i, j), &
                          dU_dx=>self%cell_gradient(1, :, i, j), &
                          dU_dy=>self%cell_gradient(2, :, i, j), &
                          x_ij=>centroid_xy(1), y_ij=>centroid_xy(2))

                  ! reconstructed_state(rho:p, point, node/midpoint, i, j)
                  U_bar(:, p, n, i, j) = cell_ave + dU_dx * (x - x_ij) + dU_dy * (y - y_ij)

                end associate

              end do
            end do
          end do
        end do

        end subroutine reconstruct_domain

        pure function estimate_gradients(self, i, j) result(gradients)
          !< Estimate the gradient of the conserved variables in the cell (i,j)
          class(second_order_reconstruction_t), intent(in) :: self
          real(rk), dimension(2, 4) :: gradients !< ([x,y], [rho,u,v,p])
          integer(ik), intent(in) :: i, j

          ! density
          gradients(:, 1) = self%estimate_single_gradient(i, j, var_idx=1)

          ! x velocity
          gradients(:, 2) = self%estimate_single_gradient(i, j, var_idx=2)

          ! y velocity
          gradients(:, 3) = self%estimate_single_gradient(i, j, var_idx=3)

          ! pressure
          gradients(:, 4) = self%estimate_single_gradient(i, j, var_idx=4)

        end function estimate_gradients

        pure function estimate_single_gradient(self, i, j, var_idx) result(grad_v)
          !< Find the gradient of a variable (v) within a cell at indices (i,j) based on the neighbor information.
          !< See Eq. 9 in https://doi.org/10.1016/j.jcp.2006.03.018. The slope limiter is set via the constructor
          !< of this derived type.

          class(second_order_reconstruction_t), intent(in) :: self
          integer(ik), intent(in) :: var_idx !< index of the variable to estimate the gradient
          integer(ik), intent(in) :: i, j !< cell index
          real(rk), dimension(2) :: grad_v !< (dV/dx, dV/dy) gradient of the variable

          associate(L=>self%limiter, &
                    U=>self%conserved_vars, v=>var_idx, &
                    volume=>self%grid%get_cell_volumes(i, j), &
                    n1=>self%grid%cell_edge_norm_vectors(:, 1, i, j), &
                    n2=>self%grid%cell_edge_norm_vectors(:, 2, i, j), &
                    n3=>self%grid%cell_edge_norm_vectors(:, 3, i, j), &
                    n4=>self%grid%cell_edge_norm_vectors(:, 4, i, j), &
                    delta_l1=>self%grid%cell_edge_lengths(1, i, j), &
                    delta_l2=>self%grid%cell_edge_lengths(2, i, j), &
                    delta_l3=>self%grid%cell_edge_lengths(3, i, j), &
                    delta_l4=>self%grid%cell_edge_lengths(4, i, j))

            grad_v = (1._rk / (2.0_rk * volume)) * &
                     (L%limit(U(v, i + 1, j) - U(v, i, j), U(v, i, j) - U(v, i - 1, j)) * (n2 * delta_l2 - n4 * delta_l4) + &
                      L%limit(U(v, i, j + 1) - U(v, i, j), U(v, i, j) - U(v, i, j - 1)) * (n3 * delta_l3 - n1 * delta_l1))
          end associate

        end function estimate_single_gradient

        end module mod_second_order_reconstruction
