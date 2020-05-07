module mod_gradients
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_floating_point_utils, only: equal

  implicit none

  private
  public :: get_smoothness, green_gauss_gradient, green_gauss_gradient_limited, modified_green_gauss_gradient
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
    real(rk), dimension(5), intent(in) :: volumes  !< (n_cells)
    real(rk), dimension(2, 4), intent(in) :: edge_normals !< ((x,y), edge)

    real(rk), dimension(4, 2) :: gradient ! ((rho, u, v, p), (grad x, grad y))
    real(rk), dimension(4, 4) :: face_prim_vars !< value of the primitive variables at the edge

    integer(ik), parameter :: n_faces = 4
    integer(ik) :: i

    gradient = 0.0_rk

    associate(v=>volumes, U=>prim_vars, U_f=>face_prim_vars)
      ! Volume weighted average on the face
      U_f(:, 1) = (U(:, 1) * v(1) + U(:, 2) * v(2)) / (v(1) + v(2))  ! bottom
      U_f(:, 2) = (U(:, 1) * v(1) + U(:, 3) * v(3)) / (v(1) + v(3))  ! right
      U_f(:, 3) = (U(:, 1) * v(1) + U(:, 4) * v(4)) / (v(1) + v(4))  ! top
      U_f(:, 4) = (U(:, 1) * v(1) + U(:, 5) * v(5)) / (v(1) + v(5))  ! left
    end associate

    do i = 1, n_faces
      ! Volume weighted average on the face
      ! x-component
      gradient(:, 1) = gradient(:, 1) + (face_prim_vars(:, i) * edge_normals(1, i) * edge_lengths(i))

      ! y-component
      gradient(:, 2) = gradient(:, 2) + (face_prim_vars(:, i) * edge_normals(2, i) * edge_lengths(i))
    end do

    gradient = gradient / volumes(1)
  end function green_gauss_gradient

  function green_gauss_gradient_limited(prim_vars, volume, edge_lengths, edge_normals) result(gradient)
    !< Summary: Estimate the gradient of the cell using the standard Green-Gauss reconstruction. This only
    !<          works well for orthogonal rectangular cells. Use the modified Green-Gauss for any other grids.
    !< Note: In the inputs below, the current cell in question is always the first in line

    real(rk), dimension(4, 5), intent(in) :: prim_vars !< ((rho, u, v, p), (center, bottom, right, top, left))
    !< ((rho, u, v, p), (current cell, neighbor_cell (1-n))) primitive
    real(rk), dimension(4), intent(in) :: edge_lengths  !< (n_edges)
    real(rk), dimension(2, 4), intent(in) :: edge_normals !< ((x,y), edge)
    real(rk), intent(in) :: volume

    real(rk), dimension(4, 2) :: gradient ! ((rho, u, v, p), (grad x, grad y))
    real(rk), dimension(4, 4) :: face_prim_vars !< value of the primitive variables at the edge

    integer(ik), parameter :: n_faces = 4
    integer(ik) :: i, l
    real(rk) :: u_center, u_other
    logical :: underflow_mode
    real(rk), dimension(4) :: phi_left_right, phi_up_down
    real(rk), parameter :: tiny_diff = 1e-10_rk

    gradient = 0.0_rk

    ! face_prim_vars indexing: 1 = bottom, 2 = right, 3 = top, 4 = left
    associate(V_edge=>face_prim_vars, &
              V_1=>prim_vars(:, 1), &
              V_2=>prim_vars(:, 2), & ! bot
              V_3=>prim_vars(:, 3), & ! right
              V_4=>prim_vars(:, 4), & ! top
              V_5=>prim_vars(:, 5))   ! left

      ! right - center, center - left
      phi_left_right = limit(V_3 - V_1, V_1 - V_5)
      ! bot - center, center - top
      phi_up_down = limit(V_2 - V_1, V_1 - V_4)

      V_edge(:, 2) = V_1 + 0.5_rk * phi_left_right  ! right
      V_edge(:, 4) = V_1 - 0.5_rk * phi_left_right  ! left

      V_edge(:, 1) = V_1 - 0.5_rk * phi_up_down  ! bottom
      V_edge(:, 3) = V_1 + 0.5_rk * phi_up_down  ! top

    end associate

    do i = 1, n_faces
      ! write(*,'(a, 4(es16.6, 2x))') 'U_f: ', face_prim_vars(:, i)
      ! x-component
      gradient(:, 1) = gradient(:, 1) + (face_prim_vars(:, i) * edge_normals(1, i) * edge_lengths(i))
      ! write(*,'(a, 4(es16.6, 2x))') 'i, x: ', (face_prim_vars(:, i) * edge_normals(1, i) * edge_lengths(i))

      ! y-component
      gradient(:, 2) = gradient(:, 2) + (face_prim_vars(:, i) * edge_normals(2, i) * edge_lengths(i))
    end do

    gradient = gradient / volume
    where(abs(gradient) < tiny_diff) gradient = 0.0_rk

    ! if (any(abs(gradient) > 1e-4_rk)) then
    !   print*, 'tiny_diff', tiny_diff
    ! write(*,'(a, 4(es16.6))') 'd/dx: ', gradient(:,1)
    ! write(*,'(a, 4(es16.6))') 'd/dy: ', gradient(:,2)
    ! print*
    ! end if

  end function green_gauss_gradient_limited

  impure elemental function limit(a_0, b_0) result(phi)
    real(rk), intent(in) :: a_0, b_0
    real(rk) :: a, b
    real(rk) :: phi
    real(rk) :: denom
    real(rk), parameter :: tiny_diff = epsilon(1.0_rk) * 5.0_rk

    a = a_0
    b = b_0
    if(abs(a_0) < tiny_diff) a = 0.0_rk
    if(abs(b_0) < tiny_diff) b = 0.0_rk

    if(abs(a - b) < tiny_diff) then
      phi = 0.0_rk
    else
      denom = a**2 + b**2
      if(denom > 0.0_rk) then
        phi = max(a * b, 0.0_rk) * (a + b) / denom
      else
        phi = 0.0_rk
      end if
    end if

    ! write(*,'(a, 6(es16.6))') 'limit: ', a, b, a**2 + b**2, abs(a-b), phi
  end function

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
