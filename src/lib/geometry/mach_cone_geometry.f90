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
  public :: mach_cone_geometry_t

  type :: mach_cone_geometry_t
    !< Type to encapsulate and calculate the geometry of the mach cone for the evolution operators
    ! private
    ! Public attributes
    real(rk), dimension(4, 2), public :: theta_ie
    real(rk), dimension(4, 2), public :: theta_ib

    real(rk), dimension(2), public :: p_prime_xy = 0.0_rk

    integer(ik), dimension(2), public :: p_prime_ij
    !< x,y location of P' or the apex of the Mach cone (global, not relative to P0)

    integer(ik), public :: neighbor_cells = 0
    integer(ik), dimension(2), public :: n_intersections = 0
    logical, dimension(4), public :: p_prime_in_cell = .false.  !< is P' inside the control volume?

    real(rk), dimension(4, 2, 4), public :: cell_conserved_vars
    !< ((rho,u,v,p), (intersection_1, intersection_2), cell_1:cell_n)

    real(rk), dimension(4), public :: reference_state !< Reference state (rho,u,v,a)
    ! real(rk), public :: tau
    real(rk), dimension(4, 4) :: reconstructed_state

    ! Private attributes
    ! real(rk), dimension(2, 2) :: theta
    ! real(rk), dimension(2) :: sin_alpha_k
    ! real(rk), dimension(2) :: cos_alpha_k
    ! real(rk) :: reference_speed_of_sound !< a_tilde
    ! real(rk), dimension(2) :: b !< b parameter (see appendix in text)
    ! real(rk), dimension(2, 2) :: l
    ! real(rk), dimension(2) :: reference_velocity !< [u_tilde, v_tilde]
    ! type(vector_t), dimension(2) :: edge_vectors !< edge vectors relative to P0
    ! type(vector_t) :: p_prime_pos_vector !< vector defining P' location relative to P0
  contains
    ! procedure, public :: initialize
    ! procedure, public :: get_p_state
    procedure, nopass, private :: calculate_l_parameter
    procedure, nopass, private :: calculate_theta_kl
    procedure, nopass, private :: determine_n_intersections
    procedure, nopass, private :: calculate_arc_thetas
    procedure, nopass, private :: determine_if_p_prime_is_in_cell
  end type

  interface mach_cone_geometry_t
    module procedure constructor
  end interface

contains

  ! subroutine initialize(self, tau)

  !   ! p_position_xy, cell_group_indices, reconstructed_state, reference_state, tau)
  !   ! // TODO: fix me!
  !   !< Implementation of the mach cone constructor
  !   class(mach_cone_geometry_t), intent(inout) :: self

  !   ! real(rk), dimension(2), intent(in) :: p_position_xy
  !   !   !< (x,y) position of P (aka Mach cone apex)

  !   ! integer(ik), dimension(:,:), intent(in) :: cell_group_indices
  !   !   !< List of the cells sharing an interface with this position P
  !   !   !< ((i,j), cell_1:cell_n)

  !   ! real(rk), dimension(:,:), intent(in) :: edge_points
  !   !   !< List of points that describe the edge vectors (x,y)
  !   !   !< The edge vectors are defined as (head:tail) by e1=([x1,y1]:[x2,y2]) and e2=([x3,y3]:[x1,y1])

  !   ! real(rk), dimension(:,:), intent(in) :: reconstructed_state
  !   !   !< Reconstructed U state for point P for each of the neighboring cells
  !   !   !< ((rho, u, v, p), cell_1:cell_n)

  !   ! real(rk), dimension(4), intent(in) :: reference_state
  !   !   !< Reference state (U_tilde) at point P (rho, u, v, p)
  !   !   !< This is typically just the average from the neighboring cells

  !   real(rk), intent(in) :: tau !< time increment

  !   integer(ik) :: k

  !   ! ! Although these are specified as a pair of x,y points, they will get shifted to the origin
  !   ! self%edge_vectors(1) = vector_t(x=edge_points, y=e1(:, 2))

  !   ! self%edge_vectors(2) = vector_t(x=e2(:, 1), y=e2(:, 2))

  !   ! self%reference_velocity = reference_state_uva(1:2)
  !   ! self%reference_speed_of_sound = reference_state_uva(3)

  !   ! ! Find the (x,y) location of the center of the Mach circle (P')
  !   ! associate(x=>e1(1, 1), y=>e1(1, 2), &
  !   !           u_tilde=>self%reference_velocity(1), v_tilde=>self%reference_velocity(2))

  !   !   ! This defines P' (x,y) globally, not with respect to P0
  !   !   self%p_prime_xy = [x - u_tilde * tau, y - v_tilde * tau]
  !   ! end associate

  !   ! ! Find the angle with respect to the x-axis for each edge vector
  !   ! do k = 1, 2
  !   !   associate(x=>self%edge_vectors(k)%x, y=>self%edge_vectors(k)%y, len=>self%edge_vectors(k)%length)

  !   !     self%sin_alpha_k(k) = y / len
  !   !     self%cos_alpha_k(k) = x / len
  !   !   end associate
  !   ! end do

  !   ! call self%determine_if_p_prime_is_in_cell()
  !   ! call self%calculate_l_parameter()
  !   ! call self%calculate_theta_kl()
  !   ! call self%determine_n_intersections()
  !   ! call self%calculate_arc_thetas()
  ! end subroutine initialize

  type(mach_cone_geometry_t) pure function constructor(tau, edge_vectors, &
                                                       reconstructed_state, reference_state, cell_indices) result(cone)

    real(rk), intent(in) :: tau  !< time increment, tau -> 0 (very small number)
    real(rk), dimension(:, :, :), intent(in) :: edge_vectors
    !< ((x,y), (tail,head), (vector_1:vector_n)); set of vectors that define the cell edges

    integer(ik), dimension(:, :), intent(in) :: cell_indices
    !< ((i,j), cell_1:cell_n); set of indices for the neighboring cells -> needed to find P' i,j index

    real(rk), dimension(:, :), intent(in) :: reconstructed_state
    !< ((rho,u,v,p), cell); reconstructed state for point P.
    !< P has a different reconstruction for each cell

    real(rk), dimension(4), intent(in) :: reference_state
    !< (rho, u, v, p); reference state of the point P

    real(rk), dimension(4) :: cone_reference_state
    !< (rho, u, v, a); reference state of the point P

    integer(ik) :: n_total_vectors  !< number of edge vectors (should be 2 or 4)
    integer(ik) :: n_neighbor_cells  !< number of cells contained in this mach cone
    integer(ik) :: neighbor_cell  !< index for looping through neighbor cells
    integer(ik), dimension(2, 4) :: edge_vector_ordering  !< index order for each set of vectors
    real(rk), dimension(2, 2) :: edge_vector_1, edge_vector_2

    type(vector_t) :: p_prime_vector
    integer(ik), dimension(2) :: n_intersections
    real(rk), dimension(2) :: theta_ib
    real(rk), dimension(2) :: theta_ie
    logical :: p_prime_in_cell
    integer(ik) :: n_arcs
    integer(ik) :: cell_ij

    cone%tau = tau
    n_total_vectors = size(edge_vectors, dim=3)
    n_neighbor_cells = n_total_vectors / 2

    ! In the cone reference state, index 4 is sound speed rather than pressure
    cone_reference_state = reference_state
    associate(a=>cone_reference_state(4), rho=>reference_state(1), &
              p=>reference_state(4), gamma=>eos%get_gamma())
      a = sqrt(gamma * p / rho)
    end associate

    ! the order of the vectors matters, mainly for the cross product to determine
    ! if P' is in the neighboring cell or not
    select case(n_total_vectors)
    case(2)
      edge_vector_ordering(:, 1:2) = reshape([[2, 1],[1, 2]], shape=[2, 2])
    case(4)
      edge_vector_ordering = reshape([[4, 1],[1, 2],[2, 3],[3, 4]], shape=[2, 4])
    case default
      error stop 'Code not set up yet to handle mach cones with other than 2 or 4 edge vectors'
    end select

    ! All the edge vectors take the point P as the origin (or tail) of their vectors
    associate(x=>edge_vectors(1, 1, 1), y=>edge_vectors(2, 1, 1), &
              u_tilde=>cone_reference_state(1), v_tilde=>cone_reference_state(2))

      ! This defines P' (x,y) globally, not with respect to P
      cone%p_prime_xy = [x - u_tilde * tau, y - v_tilde * tau]
    end associate

    ! The P' vector points from P (tail), to P' (head)
    p_prime_vector = vector_t(x=[edge_vectors(1, 1, 1), cone%p_prime_xy(1)], &
                              y=[edge_vectors(2, 1, 1), cone%p_prime_xy(2)])

    ! Loop through each neighbor cell and determine intersections and angles
    do neighbor_cell = 1, n_neighbor_cells

      cell_ij = cell_indices(:, neighbor_cell)
      edge_vector_1 = edge_vectors(:, :, edge_vector_ordering(neighbor_cell, 1))
      edge_vector_2 = edge_vectors(:, :, edge_vector_ordering(neighbor_cell, 2))

      call calculate_arc_segment(p_prime_vector=p_prime_vector, edge_vector_1=edge_vector_1, edge_vector_2=edge_vector_2, &
                                 reference_state=cone_reference_state, tau=tau, &
                                 n_intersections=n_intersections, n_arcs=n_arcs, &
                                 theta_ib=theta_ib, theta_ie=theta_ie, &
                                 p_prime_in_cell=p_prime_in_cell, cell_index=cell_ij)

      cone%n_intersections = n_intersections
      cone%theta_ib(neighbor_cell, :) = theta_ib
      cone%theta_ie(neighbor_cell, :) = theta_ie

      ! arc 1 & 2 (if more than one arc)
      cone%cell_conserved_vars(:, 1, neighbor_cell) = reconstructed_state(:, neighbor_cell)
      cone%cell_conserved_vars(:, 2, neighbor_cell) = reconstructed_state(:, neighbor_cell)

      cone%p_prime_in_cell(neighbor_cell) = p_prime_in_cell
      if(p_prime_in_cell) cone%p_prime_xy = cell_ij

    end do

  end function

  pure subroutine calculate_arc_segment(p_prime_vector, edge_vector_1, edge_vector_2, reference_state, tau, &
                                        n_intersections, n_arcs, theta_ib, theta_ie, p_prime_in_cell)
    !< For each neighbor cell, this procedure calculates the angles, number of intersections,
    !< and whether P' is in this cell
    type(vector_t), intent(in) :: p_prime_vector
    real(rk), dimension(2, 2), intent(in) :: edge_vector_1
    real(rk), dimension(2, 2), intent(in) :: edge_vector_2
    real(rk), dimension(4), intent(in) :: reference_state
    real(rk), intent(in) :: tau
    integer(ik), dimension(2), intent(out) :: n_intersections
    integer(ik), intent(out) :: n_arcs
    real(rk), dimension(2), intent(out) :: theta_ib
    real(rk), dimension(2), intent(out) :: theta_ie
    logical, intent(out) :: p_prime_in_cell

    type(vector_t), dimension(2) :: edge_vectors

    real(rk), dimension(2) :: b_k
    real(rk), dimension(2, 2) :: l_k
    real(rk), dimension(2, 2) :: theta_kl
    integer(ik) :: k
    real(rk), dimension(2) :: sin_alpha
    real(rk), dimension(2) :: cos_alpha

    ! Although these are specified as a pair of x,y points, they will get shifted to the origin
    edge_vectors(1) = vector_t(x=edge_vector_1(:, 2), y=edge_vector_1(:, 2))
    edge_vectors(2) = vector_t(x=edge_vector_2(:, 2), y=edge_vector_2(:, 2))

    do k = 1, 2
      associate(x=>edge_vectors(k)%x, y=>edge_vectors(k)%y, len=>edge_vectors(k)%length)
        sin_alpha(k) = y / len
        cos_alpha(k) = x / len
      end associate
    end do

    p_prime_in_cell = determine_if_p_prime_is_in_cell(edge_vectors(1), edge_vectors(2), p_prime_vector)
    b_k = calculate_b(sin_alpha, cos_alpha, reference_state)
    l_k = calculate_l_parameter(tau, sin_alpha, cos_alpha, reference_state, b_k)
    theta_kl = calculate_theta_kl(sin_alpha, cos_alpha, reference_state, l_k)
    n_intersections = determine_n_intersections(reference_state, b_k, l_k)
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
        n_arcs = 0
      end if
    end associate

  end subroutine calculate_arc_thetas

  pure function calculate_b(sin_alpha, cos_alpha, ref_state) result(b)
    !< Determine if there are real-value solutions to l(k)
    real(rk), dimension(2) :: b
    real(rk), dimension(2), intent(in) :: sin_alpha
    real(rk), dimension(2), intent(in) :: cos_alpha
    real(rk), dimension(4), intent(in) :: ref_state

    associate(u_tilde=>ref_state(1), v_tilde=>ref_state(2))
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
    real(rk), dimension(4), intent(in) :: ref_state
    integer(ik) :: k

    l = 0.0_rk

    do concurrent(k=1:2)
      associate(a_tilde=>ref_state(4), &
                u_tilde=>ref_state(1), v_tilde=>ref_state(2))

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
    real(rk), dimension(4), intent(in) :: ref_state
    real(rk), dimension(2, 2), intent(in) :: l_k

    real(rk) :: sine_alpha_kl
    integer(ik) :: k, l

    do concurrent(k=1:2)
      do concurrent(l=1:2)
        associate(a_tilde=>ref_state(4), u_tilde=>ref_state(1), v_tilde=>ref_state(2))

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
    real(rk), dimension(4), intent(in) :: ref_state
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
    if(((p_prime_vector.cross.edge_vector_1) <= 0.0_rk) .and. &
       ((edge_vector_2.cross.p_prime_vector) <= 0.0_rk)) then
      in_cell = .true.
    else
      in_cell = .false.
    end if
  end function determine_if_p_prime_is_in_cell

end module mod_mach_cone_geometry
