module mod_gradients
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_floating_point_utils, only: equal

  implicit none

  private
  public :: get_smoothness, green_gauss_gradient, modified_green_gauss_gradient
contains

  elemental function get_smoothness(current, left, right) result(R)
    !< Determine the smoothness of the current cell based on the neighbor cells.
    !< for x it's left/right, for y it's top/bottom

    real(rk) :: R  !< Smoothness -> R = (U_{i+1} - u_{i}) /( u_{i} - u_{i-1})
    real(rk), intent(in) :: current
    real(rk), intent(in) :: left
    real(rk), intent(in) :: right
    real(rk) :: denomenator, numerator

    numerator = right - current
    denomenator = current - left

    if(equal(denomenator, 0.0_rk)) then
      R = 0.0_rk
    else
      R = numerator / denomenator
    end if
  end function get_smoothness

  function green_gauss_gradient(prim_vars, volumes, edge_lengths, edge_normals) result(gradient)
    !< Summary: Estimate the gradient of the cell using the standar Green-Gauss reconstruction. This only
    !<          works well for orthogonal rectangular cells. Use the modified Green-Gauss for any other grids.
    !< Note: In the inputs below, the current cell in question is always the first in line

    real(rk), dimension(4, 5), intent(in) :: prim_vars
    !< ((rho, u, v, p), (current cell, neighbor_cell (1-n))) primitive
    real(rk), dimension(4), intent(in) :: edge_lengths  !< (n_edges)
    real(rk), dimension(4), intent(in) :: volumes  !< (n_cells)
    real(rk), dimension(2, 4), intent(in) :: edge_normals !< ((x,y), edge)

    real(rk), dimension(4, 2) :: gradient ! ((rho, u, v, p), (grad x, grad y))
    real(rk), dimension(4, 4) :: face_prim_vars !< value of the primitive variables at the edge

    integer(ik), parameter :: n_faces = 4
    integer(ik) :: i

    gradient = 0.0_rk

    ! associate(v => volumes, U => prim_vars, U_f => face_prim_vars)
    !   ! Volume weighted average on the face
    !   U_f(:,1) = (U(:,1) * v(1) + U(:,2) * v(2)) / (v(1) + v(2))  ! down
    !   U_f(:,2) = (U(:,1) * v(1) + U(:,3) * v(3)) / (v(1) + v(3))  ! right
    !   U_f(:,3) = (U(:,1) * v(1) + U(:,4) * v(4)) / (v(1) + v(4))  ! up
    !   U_f(:,4) = (U(:,1) * v(1) + U(:,5) * v(5)) / (v(1) + v(5))  ! left

    !   ! Smoothness indicators
    !   R_left_right = get_smoothness(U(:,1),U(:,5), U(:,3))
    !   R_up_down = get_smoothness(U(:,1),U(:,5), U(:,3))

    ! end associate

    do i = 1, n_faces
      ! Volume weighted average on the face
      face_prim_vars(:, i) = (prim_vars(:, 1) * volumes(1) + prim_vars(:, i + 1) * volumes(i + 1)) / (volumes(1) + volumes(i + 1))
      ! x-component
      gradient(:, 1) = gradient(:, 1) + face_prim_vars(:, i) * edge_normals(1, i) * edge_lengths(i)

      ! y-component
      gradient(:, 2) = gradient(:, 2) + face_prim_vars(:, i) * edge_normals(2, i) * edge_lengths(i)
    end do

    gradient = gradient / volumes(1)
  end function green_gauss_gradient

  function modified_green_gauss_gradient(prim_vars, centroids, volumes, &
                                         edge_lengths, edge_midpoints, edge_normals) result(gradient)

    !< Estimate the gradient of the cell using the Green-Gauss reconstruction. This only
    !< works well for orthogonal rectangular cells. Use the modified Green-Gauss for any other grids
    !< Note in these arrays, the current cell in question is always the first in line

    real(rk), dimension(:, :), intent(in) :: prim_vars
    !< ((rho, u, v, p), (current cell, neighbor_cell (1-n))) primitive
    real(rk), dimension(:, :), intent(in) :: centroids !< ((x,y), (current cell, neighbor_cell (1-n)))
    real(rk), dimension(:), intent(in) :: edge_lengths  !< (n_edges)
    real(rk), dimension(:), intent(in) :: volumes  !< (n_cells)
    real(rk), dimension(:, :), intent(in) :: edge_midpoints !< ((x,y), edge)
    real(rk), dimension(:, :), intent(in) :: edge_normals !< ((x,y), edge)
    real(rk), dimension(4, 2) :: gradient

    real(rk), dimension(2) :: r_f
    integer(ik) :: n_faces, f

    real(rk), parameter :: eps = 1.0e-8_rk ! convergence criteria

    gradient = 0.0_rk
    ! n_faces = size(prim_vars, dim=2) - 1
    ! allocate(face_prim_vars(4, n_faces))

    ! grad_c = 0.0_rk
    ! grad_nb = 0.0_rk

    ! k = 1

    ! iter_loop: do  ! loop until converged
    !   do f = 2, n_faces
    !     associate(x_c=>centroids(1), x_nb=>centroids(f), &
    !               x_f=>edge_midpoints(f), n_f=>edge_normals(f), &
    !               phi_c=>prim_vars(:, 1), phi_nb=>prim_vars(:, f))

    !       ! vector from edge midpoint to neighbor centroid
    !       r_f = (x_nb - x_f) / norm2(x_nb - x_f)

    !       ! orthogonality scaling factor
    !       alpha = dot_product(n_f, r_f)

    !       ! distance from the current centroid to neighbor centroid
    !       delta_r = norm2(x_nb - x_c)

    !       orth_contrib = alpha * (phi_nb - phi_c) / delta_r
    !       non_orth_contrib = 0.5_rk * dot_product((grad_c + grad_nb),(n_f - alpha * r_f))
    !       dphi_dn = orthogonal_contrib + non_orthogonal_contrib

    !     end associate
    !   end do

    !   if(norm2() < eps) exit
    !   k = k + 1
    ! end do

    ! gradient = grad_c
    ! deallocate(face_prim_vars)

  end function modified_green_gauss_gradient
end module mod_gradients
