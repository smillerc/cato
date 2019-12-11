module mod_mach_cone_geometry
  !< Summary: This module contains the procedures to calculate the geometry
  !< of the Mach cone for the solution of the Euler equations using the
  !< characteristic lines. Refer to the Appendix of DOI: 10.1016/j.jcp.2009.04.001. The
  !< entire purpose of this module is to calculate the quantities \(\theta_{i,b}\)
  !< and \(\theta_{i,e}\)

  use iso_fortran_env, only: ik => int32, rk => real64
  use math_constants, only: pi
  use mod_vector, only: vector_t, operator(.unitnorm.), operator(.dot.), operator(.cross.), angle_between
  use mod_eos, only: eos

  implicit none

  private
  public :: mach_cone_geometry_t, new_cone

  type :: mach_cone_geometry_t
    !< Type to encapsulate and calculate the geometry of the mach cone for the evolution operators

    real(rk), dimension(4, 2) :: theta_ie = 0.0_rk
    !< ((cell_1:cell_4), (arc_1, arc_2)); arc end angle

    real(rk), dimension(4, 2) :: theta_ib = 0.0_rk
    !< ((cell_1:cell_4), (arc_1, arc_2)); arc begin angle

    real(rk), dimension(2) :: p_prime_xy = 0.0_rk
    !< (x,y); Location of P'

    real(rk), dimension(2) :: p_xy = 0.0_rk
    !< (x,y); Location of P

    integer(ik), dimension(2) :: p_prime_ij = [0, 0]
    !< x,y location of P' or the apex of the Mach cone (global, not relative to P0)

    integer(ik) :: n_neighbor_cells = 0
    !< How many cells influence this cone (up to 4 for now)

    real(rk) :: radius = 0.0_rk
    !< Radius of the cone

    integer(ik), dimension(4) :: n_arcs = 0
    !< (cell_1:cell_4); Number of arcs that each cell contributes (0, 1, or 2)

    logical, dimension(4) :: p_prime_in_cell = .false.
    !< (cell_1:cell_4); is P' inside the control volume?

    real(rk), dimension(4, 2, 4) :: cell_conserved_vars = 0.0_rk
    !< ((rho,u,v,p), (arc_1, arc_2), cell_1:cell_n)

    real(rk), dimension(4) :: reference_state = 0.0_rk
    !< (rho,u,v,a); Reference state (local cell average of U)

    real(rk) :: tau = 0.0_rk
    !< Time evolution increment

    real(rk), dimension(4, 4) :: reconstructed_state = 0.0_rk
    !< ((rho,u,v,p), (edge_vector_1:edge_vector_4))

  contains
    procedure, nopass, private :: calculate_l_parameter
    procedure, nopass, private :: calculate_theta_kl
    procedure, nopass, private :: determine_n_intersections
    procedure, nopass, private :: calculate_arc_thetas
    procedure, nopass, private :: determine_if_p_prime_is_in_cell
  end type

  interface mach_cone_geometry_t
    module procedure new_cone
  end interface

  ! Corner/midpoint index convention         Cell Indexing convention
  ! --------------------------------         ------------------------
  !
  !   C----M----C3----M----C
  !   |         |         |                             E3
  !   O    x    O    x    O                      N4-----M3----N3
  !   |         |         |                      |            |
  !   C4----M---C0----M----C2                E4  M4     C     M2  E2
  !   |         |         |                      |            |
  !   O    x    O    x    O                      N1----M1----N2
  !   |         |         |                            E1
  !   C----M----C1----M----C
  !
  ! For left/right midpoints, the edge vectors go left then right.
  ! The neighboring cells are above (i,j) and below (i,j-1)
  ! For quad cells, N - corner, M - midpoint, E - edge

  ! Corner vector order should be (vec1 (C0->C1), vec2 (C0->C2), vec3 (C0->C3), vec4 (C0->C4))

contains

  pure function new_cone(tau, edge_vectors, reconstructed_state, reference_state, cell_indices)
    !< Constructor for the Mach cone type

    type(mach_cone_geometry_t) :: new_cone

    real(rk), intent(in) :: tau
    !< time increment, tau -> 0 (very small number)

    real(rk), dimension(:, :, :), intent(in) :: edge_vectors
    !< ((x,y), (tail,head), (vector_1:vector_n)); set of vectors that define the cell edges

    integer(ik), dimension(:, :), intent(in) :: cell_indices
    !< ((i,j), cell_1:cell_n); set of indices for the neighboring cells -> needed to find P' i,j index

    real(rk), dimension(:, :), intent(in) :: reconstructed_state
    !< ((rho,u,v,p), cell); reconstructed state for point P.
    !< P has a different reconstruction for each cell

    real(rk), dimension(4), intent(in) :: reference_state  !< (rho, u, v, p)
    !< (rho, u, v, p); reference state of the point P

    ! real(rk), dimension(4) :: cone_reference_state
    !< (rho, u, v, a); reference state of the point P

    integer(ik) :: n_total_vectors  !< number of edge vectors (should be 2 or 4)
    integer(ik) :: neighbor_cell  !< index for looping through neighbor cells
    integer(ik), dimension(2, 4) :: edge_vector_ordering  !< index order for each set of vectors
    real(rk), dimension(2, 2) :: edge_vector_1, edge_vector_2

    type(vector_t) :: p_prime_vector
    integer(ik), dimension(2) :: n_intersections
    real(rk), dimension(2) :: theta_ib
    real(rk), dimension(2) :: theta_ie
    logical :: p_prime_in_cell
    integer(ik) :: n_arcs
    integer(ik), dimension(2) :: cell_ij

    edge_vector_ordering = 0
    n_total_vectors = size(edge_vectors, dim=3)
    new_cone%tau = tau

    ! In the cone reference state, index 4 is sound speed rather than pressure
    !< (rho, u, v, a) vs  !< (rho, u, v, p)
    new_cone%reference_state = reference_state
    ! //TODO: Move speed of sound calculation to the EOS module
    associate(a=>new_cone%reference_state(4), rho=>reference_state(1), &
              p=>reference_state(4), gamma=>eos%get_gamma())
      a = sqrt(gamma * p / rho)
      new_cone%radius = a * tau
    end associate

    ! the order of the vectors matters, mainly for the cross product to determine
    ! if P' is in the neighboring cell or not
    select case(n_total_vectors)
    case(2)
      edge_vector_ordering(:, 1:2) = reshape([[2, 1],[1, 2]], shape=[2, 2])
      new_cone%n_neighbor_cells = 2
    case(4)
      edge_vector_ordering = reshape([[4, 1],[1, 2],[2, 3],[3, 4]], shape=[2, 4])
      new_cone%n_neighbor_cells = 4
    case default
      error stop 'Code not set up yet to handle mach cones with other than 2 or 4 edge vectors'
    end select

    ! All the edge vectors take the point P as the origin (or tail) of their vectors
    associate(x=>edge_vectors(1, 1, 1), y=>edge_vectors(2, 1, 1), &
              u_tilde=>new_cone%reference_state(2), v_tilde=>new_cone%reference_state(3))

      ! This defines P' (x,y) globally, not with respect to P
      new_cone%p_xy = [x, y]
      new_cone%p_prime_xy = [x - u_tilde * tau, y - v_tilde * tau]
    end associate

    ! The P' vector points from P (tail), to P' (head)
    p_prime_vector = vector_t(x=[edge_vectors(1, 1, 1), new_cone%p_prime_xy(1)], &
                              y=[edge_vectors(2, 1, 1), new_cone%p_prime_xy(2)])

    ! Loop through each neighbor cell and determine intersections and angles
    do neighbor_cell = 1, new_cone%n_neighbor_cells
      p_prime_in_cell = .false.

      cell_ij = cell_indices(:, neighbor_cell)  ! cell_indices is indexed via ((i,j), cell_1:cell_n)

      ! indexing for edge_vectors is ((x,y), (tail,head), (vector_1:vector_n))
      ! edge_vector_ordering ((), ())
      edge_vector_1 = edge_vectors(:, :, edge_vector_ordering(1, neighbor_cell))
      edge_vector_2 = edge_vectors(:, :, edge_vector_ordering(2, neighbor_cell))

      ! print*, 'Neighbor cell: ', neighbor_cell
      call calculate_arc_segment(p_prime_vector=p_prime_vector, &
                                 edge_vector_1=edge_vector_1, edge_vector_2=edge_vector_2, &
                                 reference_state=new_cone%reference_state, tau=tau, &
                                 n_arcs=n_arcs, &
                                 theta_ib=theta_ib, theta_ie=theta_ie, &
                                 p_prime_in_cell=p_prime_in_cell)

      new_cone%n_arcs(neighbor_cell) = n_arcs
      new_cone%theta_ib(neighbor_cell, :) = theta_ib
      new_cone%theta_ie(neighbor_cell, :) = theta_ie

      ! arc 1 & 2 (if more than one arc)
      new_cone%cell_conserved_vars(:, 1, neighbor_cell) = reconstructed_state(:, neighbor_cell)
      new_cone%cell_conserved_vars(:, 2, neighbor_cell) = reconstructed_state(:, neighbor_cell)

      new_cone%p_prime_in_cell(neighbor_cell) = p_prime_in_cell
      if(p_prime_in_cell) new_cone%p_prime_ij = cell_ij

    end do

  end function

  pure subroutine calculate_arc_segment(p_prime_vector, edge_vector_1, edge_vector_2, reference_state, tau, &
                                       n_arcs, theta_ib, theta_ie, p_prime_in_cell)
    !< For each neighbor cell, this procedure calculates the angles, number of intersections,
    !< and whether P' is in this cell
    type(vector_t), intent(in) :: p_prime_vector  !< vector pointing to P' from P (this should be shifted to the origin by now)
    real(rk), dimension(2, 2), intent(in) :: edge_vector_1 !< 1st edge vector ((x,y), (head,tail))
    real(rk), dimension(2, 2), intent(in) :: edge_vector_2  !< 2nd edge vector ((x,y), (head,tail))
    real(rk), dimension(4), intent(in) :: reference_state !< (rho,u,v,a)
    real(rk), intent(in) :: tau  !< time increment
    integer(ik), dimension(2) :: n_intersections !< # of intersections
    integer(ik), intent(out) :: n_arcs !< # of arcs for each cell (0, 1, or 2)
    real(rk), dimension(2), intent(out) :: theta_ib
    real(rk), dimension(2), intent(out) :: theta_ie
    logical, intent(out) :: p_prime_in_cell

    type(vector_t), dimension(2) :: edge_vectors
    real(rk), dimension(2) :: b_k
    real(rk), dimension(2, 2) :: l_k
    real(rk), dimension(2, 2) :: theta_kl
    real(rk), dimension(2) :: sin_alpha
    real(rk), dimension(2) :: cos_alpha
    integer(ik) :: k

    ! Although these are specified as a pair of x,y points, they will get shifted to the origin
    edge_vectors(1) = vector_t(x=[edge_vector_1(1, 1), edge_vector_1(1, 2)], &
                               y=[edge_vector_1(2, 1), edge_vector_1(2, 2)])

    edge_vectors(2) = vector_t(x=[edge_vector_2(1, 1), edge_vector_2(1, 2)], &
                               y=[edge_vector_2(2, 1), edge_vector_2(2, 2)])

    do k = 1, 2
      associate(x=>edge_vectors(k)%x, y=>edge_vectors(k)%y, len=>edge_vectors(k)%length)
        sin_alpha(k) = (y / len)
        cos_alpha(k) = (x / len)
      end associate
    end do

    p_prime_in_cell = determine_if_p_prime_is_in_cell(edge_vectors(1), edge_vectors(2), p_prime_vector)
    b_k = calculate_b(sin_alpha, cos_alpha, reference_state)
    ! print*, 'b_k', b_k
    l_k = calculate_l_parameter(tau, sin_alpha, cos_alpha, reference_state, b_k)
    ! print*, 'l_k', l_k
    theta_kl = calculate_theta_kl(sin_alpha, cos_alpha, reference_state, l_k)
    ! print*, 'theta_kl', theta_kl
    n_intersections = determine_n_intersections(reference_state, b_k, l_k)
    ! print*, 'n_intersections',n_intersections
    call calculate_arc_thetas(theta_kl, n_intersections, p_prime_in_cell, n_arcs, theta_ib, theta_ie)

  end subroutine

  pure subroutine calculate_arc_thetas(theta, n_intersections, p_prime_in_cell, n_arcs, theta_ib, theta_ie)
    !< Implementation of the logic used to find the arc angles theta ib and ie (begin and end)

    real(rk), dimension(2, 2), intent(in) :: theta !< this is the theta k,l term in the appendix,
    integer(ik), dimension(2), intent(in) :: n_intersections
    logical, intent(in) :: p_prime_in_cell

    real(rk), dimension(2), intent(out) :: theta_ib
    real(rk), dimension(2), intent(out) :: theta_ie
    integer(ik), intent(out) :: n_arcs

    associate(n1=>n_intersections(1), n2=>n_intersections(2))

      if(n1 == 0 .and. n2 == 0 .and. p_prime_in_cell) then
        n_arcs = 1
        theta_ib = 0.0_rk
        theta_ie = 2 * pi

      else if(n1 == 1 .and. n2 == 1) then
        n_arcs = 1
        theta_ib = theta(1, 2)

        if(theta(2, 2) >= theta(1, 2)) then
          theta_ie = theta(2, 2)
        else
          theta_ie = theta(2, 2) + 2 * pi
        endif

      else if(n1 == 2 .and. n2 == 0) then
        n_arcs = 1
        theta_ib = theta(1, 2)

        if(theta(1, 1) >= theta(1, 2)) then
          theta_ie = theta(1, 1)
        else
          theta_ie = theta(1, 1) + 2 * pi
        endif

      else if(n1 == 0 .and. n2 == 2) then
        n_arcs = 1
        theta_ib = theta(2, 1)

        if(theta(2, 2) >= theta(2, 1)) then
          theta_ie = theta(2, 2)
        else
          theta_ie = theta(2, 2) + 2 * pi
        endif

      else if(n1 == 2 .and. n2 == 2) then
        ! There are 2 arcs for each edge vector
        n_arcs = 2

        theta_ib(1) = theta(1, 2)
        theta_ib(2) = theta(2, 1)

        ! 1st arc
        if(theta(2, 2) >= theta(1, 2)) then
          theta_ie(1) = theta(2, 2)
        else
          theta_ie(1) = theta(2, 2) + 2 * pi
        endif

        ! 2nd arc
        if(theta(1, 1) >= theta(2, 1)) then
          theta_ie(2) = theta(1, 1)
        else
          theta_ie(2) = theta(1, 1) + 2 * pi
        endif
      else
        ! This cell doesn't contribute at all (setting to 0 for the sake of clean-ness)...
        theta_ib = 0.0_rk
        theta_ie = 0.0_rk
        n_arcs = 0
      end if
    end associate

  end subroutine calculate_arc_thetas

  pure function calculate_b(sin_alpha, cos_alpha, ref_state) result(b)
    !< Determine if there are real-value solutions to l(k)
    real(rk), dimension(2) :: b
    real(rk), dimension(2), intent(in) :: sin_alpha
    real(rk), dimension(2), intent(in) :: cos_alpha
    real(rk), dimension(4), intent(in) :: ref_state  !< (rho, u, v, a)

    associate(u_tilde=>ref_state(2), v_tilde=>ref_state(3))
      b(1) = abs(u_tilde * sin_alpha(1) + v_tilde * cos_alpha(1))
      b(2) = abs(u_tilde * sin_alpha(2) + v_tilde * cos_alpha(2))
    end associate

  end function calculate_b

  pure function calculate_l_parameter(tau, sin_alpha, cos_alpha, ref_state, b) result(l)
    !< Calculate the \( \ell_{k} \)] parameter

    real(rk), dimension(2, 2) :: l
    real(rk), intent(in) :: tau
    real(rk), dimension(2), intent(in) :: b
    real(rk), dimension(2), intent(in) :: sin_alpha
    real(rk), dimension(2), intent(in) :: cos_alpha
    real(rk), dimension(4), intent(in) :: ref_state !< (rho, u, v, a)
    integer(ik) :: k

    l = 0.0_rk

    do concurrent(k=1:2)
      associate(a_tilde=>ref_state(4), &
                u_tilde=>ref_state(2), v_tilde=>ref_state(3))

        if(b(k) <= a_tilde) then
          l(k, 1) = -tau * (u_tilde * cos_alpha(k) + v_tilde * sin_alpha(k)) - &
                    tau * sqrt(a_tilde**2 - ((u_tilde * sin_alpha(k) - v_tilde * cos_alpha(k))**2))

          l(k, 2) = -tau * (u_tilde * cos_alpha(k) + v_tilde * sin_alpha(k)) + &
                    tau * sqrt(a_tilde**2 - ((u_tilde * sin_alpha(k) - v_tilde * cos_alpha(k))**2))
        end if

      end associate
    end do
  end function calculate_l_parameter

  pure function calculate_theta_kl(sin_alpha, cos_alpha, ref_state, l_k) result(theta)
    !< Implementation of the \(\theta_{k,\ell}\) term in the Appendix of the text

    real(rk), dimension(2, 2) :: theta
    real(rk), dimension(2), intent(in) :: sin_alpha
    real(rk), dimension(2), intent(in) :: cos_alpha
    real(rk), dimension(4), intent(in) :: ref_state  !< (rho, u, v, a)
    real(rk), dimension(2, 2), intent(in) :: l_k

    real(rk) :: sine_alpha_kl
    integer(ik) :: k, l

    do concurrent(k=1:2)
      do concurrent(l=1:2)
        associate(a_tilde=>ref_state(4), u_tilde=>ref_state(2), v_tilde=>ref_state(3))

          sine_alpha_kl = v_tilde + l_k(k, l) * sin_alpha(k) / a_tilde

          theta(k, l) = pi + sign(1.0_rk, sine_alpha_kl) * &
                        (acos((u_tilde + l_k(k, l) * cos_alpha(k)) / a_tilde) - pi)
        end associate
      end do
    end do
  end function calculate_theta_kl

  pure function determine_n_intersections(ref_state, b, l) result(n)
    !< Summary: Find the number of intersections based on the input parameters
    !< Refer to the appendix for this particular logic
    integer(ik), dimension(2) :: n
    real(rk), dimension(4), intent(in) :: ref_state  !< (rho, u, v, a)
    real(rk), dimension(2), intent(in) :: b
    real(rk), dimension(2, 2), intent(in) :: l
    integer(ik) :: k

    do concurrent(k=1:2)
      associate(a_tilde=>ref_state(4))
        if((b(k) > a_tilde) .or. &
           ((l(k, 1) <= l(k, 2)) .and. (l(k, 2) < 0.0_rk))) then
          n(k) = 0
        else if((l(k, 2) >= 0.0_rk) .and. (l(k, 1) < 0.0_rk)) then
          n(k) = 1
        else
          n(k) = 2
        end if
      end associate
    end do

  end function determine_n_intersections

  logical pure function determine_if_p_prime_is_in_cell(edge_vector_1, edge_vector_2, p_prime_vector) result(in_cell)
    !< Implementation of whether the P' point is inside the current cell/control volume. This
    !< uses the cross product of 2 vectors in 2d, which gives a scalar
    type(vector_t), intent(in) :: p_prime_vector
    type(vector_t), intent(in) :: edge_vector_1, edge_vector_2

    ! In the text, these inequalities are swapped, but I define the P' position vector
    ! as starting at P0 and going to P'
    ! write(*,*) 'p_prime_vector: ', p_prime_vector
    ! write(*,*) 'edge_vector_1:  ', edge_vector_1
    ! write(*,*) 'edge_vector_2:  ', edge_vector_2

    ! In the text is has >= instead of <= for some reason, but the following works
    ! like it's supposed to.
    if((p_prime_vector.cross.edge_vector_1) <= 0.0_rk .and. &
       (edge_vector_2.cross.p_prime_vector) <= 0.0_rk) then
      in_cell = .true.
    else
      in_cell = .false.
    end if

  end function determine_if_p_prime_is_in_cell

end module mod_mach_cone_geometry
