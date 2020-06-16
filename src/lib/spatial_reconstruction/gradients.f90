module mod_gradients
  !< Summary: Provide procedures to calculate the gradient of a particular variable
  !< Date: 05/08/2020
  !< Author: Sam Miller
  !< Notes:
  !< References:
  !      [1]

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_floating_point_utils, only: neumaier_sum
  use mod_grid, only: grid_t
  use mod_globals, only: n_ghost_layers, plot_gradients
  use hdf5_interface, only: hdf5_file

  implicit none
  real(rk), parameter :: EPS = 2.0_rk * epsilon(1.0_rk)
  logical, parameter :: filter_small_values = .true.
  private
  public :: green_gauss_gradient

contains

  subroutine green_gauss_gradient(edge_vars, lbounds, grid, grad_x, grad_y, name, stage_name)
    !< Estimate the slope-limited gradient of the primitive variables in the cell (i,j). This assumes
    !< a quadrilateral structured grid
    class(grid_t), intent(in) :: grid
    integer(ik), dimension(3), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), contiguous, intent(in) :: edge_vars
    !< ((edge 1:n), i, j); primitive variable interpolated on the cell interface/edge

    real(rk), dimension(lbounds(2):, lbounds(3):), intent(inout) :: grad_x !< (i, j); gradient in the x-direction
    real(rk), dimension(lbounds(2):, lbounds(3):), intent(inout) :: grad_y !< (i, j); gradient in the y-direction

    character(len=*), intent(in) :: name !< what variable are we finding the gradient of?
    character(len=*), intent(in) :: stage_name !< what stage in the time integration are we
    logical :: plot !< plot to file?

    integer(ik) :: i, j
    integer(ik) :: ilo, ihi, jlo, jhi

    character(len=50) :: filename = ''

    real(rk), dimension(4) :: edge_lengths  !< length of each face
    real(rk), dimension(4) :: n_x  !< normal vectors of each face
    real(rk), dimension(4) :: n_y  !< normal vectors of each face
    real(rk) :: d_dx, d_dy
    real(rk) :: diff

    ilo = lbound(edge_vars, dim=2)
    ihi = ubound(edge_vars, dim=2)
    jlo = lbound(edge_vars, dim=3)
    jhi = ubound(edge_vars, dim=3)

    !$omp parallel default(none), &
    !$omp firstprivate(ilo, ihi, jlo, jhi, n_ghost_layers) &
    !$omp private(i, j) &
    !$omp private(n_x, n_y, edge_lengths, d_dx, d_dy, diff) &
    !$omp shared(grad_x, grad_y, edge_vars, grid)
    !$omp do
    do j = jlo + n_ghost_layers, jhi - n_ghost_layers
      do i = ilo + n_ghost_layers, ihi - n_ghost_layers

        ! Edge (face) interface data
        edge_lengths(1) = grid%cell_edge_lengths(1, i, j)  ! bottom
        edge_lengths(2) = grid%cell_edge_lengths(2, i, j)  ! right
        edge_lengths(3) = grid%cell_edge_lengths(3, i, j)  ! top
        edge_lengths(4) = grid%cell_edge_lengths(4, i, j)  ! left

        if(abs(edge_lengths(3) - edge_lengths(1)) < EPS) edge_lengths(3) = edge_lengths(1)
        if(abs(edge_lengths(4) - edge_lengths(2)) < EPS) edge_lengths(4) = edge_lengths(2)

        diff = maxval(edge_lengths) - minval(edge_lengths)
        if(diff < EPS) edge_lengths = maxval(edge_lengths)

        n_x(1) = grid%cell_edge_norm_vectors(1, 1, i, j)  ! bottom
        n_y(1) = grid%cell_edge_norm_vectors(2, 1, i, j)  ! bottom
        n_x(2) = grid%cell_edge_norm_vectors(1, 2, i, j)  ! right
        n_y(2) = grid%cell_edge_norm_vectors(2, 2, i, j)  ! right
        n_x(3) = grid%cell_edge_norm_vectors(1, 3, i, j)  ! top
        n_y(3) = grid%cell_edge_norm_vectors(2, 3, i, j)  ! top
        n_x(4) = grid%cell_edge_norm_vectors(1, 4, i, j)  ! left
        n_y(4) = grid%cell_edge_norm_vectors(2, 4, i, j)  ! left

        d_dx = sum(edge_vars(:, i, j) * n_x * edge_lengths)
        d_dy = sum(edge_vars(:, i, j) * n_y * edge_lengths)
        ! print*, 'edge n_x', n_x
        ! print*, 'edge n_y', n_y
        ! print*, 'edge_lengths', edge_lengths
        ! print*, 'edge vars', edge_vars(:,i,j)
        ! print*, 'd_dx', d_dx
        ! print*

        if(filter_small_values) then
          if(abs(d_dx) < EPS) then
            grad_x(i, j) = 0.0_rk
          else
            grad_x(i, j) = d_dx / grid%cell_volume(i, j)
          end if

          if(abs(d_dy) < EPS) then
            grad_y(i, j) = 0.0_rk
          else
            grad_y(i, j) = d_dy / grid%cell_volume(i, j)
          end if
        else
          grad_x(i, j) = d_dx / grid%cell_volume(i, j)
          grad_y(i, j) = d_dy / grid%cell_volume(i, j)
        end if

      end do
    end do
    !$omp end do
    !$omp end parallel

    if(plot_gradients) then
      io: block
        type(hdf5_file) :: h5
        logical :: file_exists

        filename = 'gradients_'//trim(stage_name)//'.h5'

        inquire(file=filename, exist=file_exists)

        if(file_exists) then
          call h5%initialize(filename=filename, status='old', action='rw', comp_lvl=6)
        else
          call h5%initialize(filename=filename, status='new', action='rw', comp_lvl=6)
        end if

        call h5%add("d_dx_"//name, grad_x)
        call h5%add("d_dy_"//name, grad_y)
        call h5%finalize()
      end block io
    end if
  end subroutine green_gauss_gradient

end module mod_gradients
