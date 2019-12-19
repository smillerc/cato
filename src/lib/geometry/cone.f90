module mod_cone
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use math_constants, only: pi
  use mod_floating_point_utils, only: near_zero
  use mod_eos, only: eos
  use mod_vector, only: vector_t, operator(.unitnorm.), operator(.dot.), operator(.cross.), angle_between

  implicit none

  type :: cone_t
    !< Type to encapsulate and calculate the geometry of the mach cone for the evolution operators

    real(rk), dimension(4, 2) :: theta_ie = 0.0_rk
    !< ((cell_1:cell_4), (arc_1, arc_2)); arc end angle

    real(rk), dimension(4, 2) :: theta_ib = 0.0_rk
    !< ((cell_1:cell_4), (arc_1, arc_2)); arc begin angle

    real(rk), dimension(2) :: p_prime_xy = 0.0_rk
    !< (x,y); Location of P'

    real(rk), dimension(2) :: p_xy = 0.0_rk
    !< (x,y); Location of P

    integer(ik), dimension(2) :: p_prime_ij = 0
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
    !< ((rho,u,v,p), (cell_1:cell_4))

  contains
    procedure, nopass, private :: determine_if_p_prime_is_in_cell
    procedure, pass :: write => write_cone
    generic, public :: write(formatted) => write
  end type

contains

  pure function new_cone(tau, edge_vectors, reconstructed_state, reference_state, cell_indices)
    !< Constructor for the Mach cone type

    type(cone_t) :: new_cone

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

    integer(ik) :: n_total_vectors  !< number of edge vectors (should be 2 or 4)
    integer(ik) :: neighbor_cell  !< index for looping through neighbor cells
    integer(ik), dimension(2, 4) :: edge_vector_ordering  !< index order for each set of vectors
    real(rk), dimension(2, 2, 2) :: single_cell_edge_vectors

    type(vector_t) :: p_prime_vector
    logical :: p_prime_in_cell
    integer(ik) :: n_arcs
    integer(ik), dimension(2) :: cell_ij
    real(rk), dimension(2, 2) :: theta_ib_ie

    p_prime_in_cell = .false.
    n_arcs = 0
    cell_ij = 0
    theta_ib_ie = 0.0_rk
    n_total_vectors = 0
    neighbor_cell = 0
    edge_vector_ordering = 0
    single_cell_edge_vectors = 0.0_rk
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

    new_cone%p_prime_xy = 0.0_rk
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
      single_cell_edge_vectors(:, :, 1) = edge_vectors(:, :, edge_vector_ordering(1, neighbor_cell))
      single_cell_edge_vectors(:, :, 2) = edge_vectors(:, :, edge_vector_ordering(2, neighbor_cell))

      p_prime_in_cell = determine_if_p_prime_is_in_cell(single_cell_edge_vectors, p_prime_vector)
      call get_arc_segments(lines=single_cell_edge_vectors, origin_in_cell=p_prime_in_cell, &
                            circle_xy=new_cone%p_prime_xy, circle_radius=new_cone%radius, &
                            arc_segments=theta_ib_ie, n_arcs=n_arcs)

      new_cone%n_arcs(neighbor_cell) = n_arcs
      new_cone%theta_ib(neighbor_cell, :) = theta_ib_ie(1, :)
      new_cone%theta_ie(neighbor_cell, :) = theta_ib_ie(2, :)

      ! arc 1 & 2 (if more than one arc)
      new_cone%cell_conserved_vars(:, 1, neighbor_cell) = reconstructed_state(:, neighbor_cell)
      new_cone%cell_conserved_vars(:, 2, neighbor_cell) = reconstructed_state(:, neighbor_cell)

      new_cone%p_prime_in_cell(neighbor_cell) = p_prime_in_cell
      if(p_prime_in_cell) new_cone%p_prime_ij = cell_ij
      ! print *
    end do

  end function

  subroutine write_cone(self, unit, iotype, v_list, iostat, iomsg)
    !< Implementation of `write(*,*) vector_t`

    class(cone_t), intent(in) :: self     !< cone class
    integer, intent(in) :: unit           !< input/output unit
    character(*), intent(in) :: iotype    !< LISTDIRECTED or DTxxx
    integer, intent(in) :: v_list(:)      !< parameters from fmt spec.
    integer, intent(out) :: iostat        !< non zero on error, etc.
    character(*), intent(inout) :: iomsg  !< define if iostat non zero.

    integer(ik) :: i

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) 'Mach Cone Details'//new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) '================='//new_line('a')

    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "Tau: ", self%tau, new_line('a')
    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "Radius: ", self%radius, new_line('a')
    write(unit, '(a, 2(f7.3,1x), a)', iostat=iostat, iomsg=iomsg) "P (x,y): (", self%p_xy, ")"//new_line('a')
    write(unit, '(a, 2(f7.3,1x), a)', iostat=iostat, iomsg=iomsg) "P'(x,y): (", self%p_prime_xy, ")"//new_line('a')

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a, f0.3, a)', iostat=iostat, iomsg=iomsg) "Reference density:     ", self%reference_state(1), new_line('a')
    write(unit, '(a, f0.3, a)', iostat=iostat, iomsg=iomsg) "Reference x velocity:  ", self%reference_state(2), new_line('a')
    write(unit, '(a, f0.3, a)', iostat=iostat, iomsg=iomsg) "Reference y velocity:  ", self%reference_state(3), new_line('a')
    write(unit, '(a, f0.3, a)', iostat=iostat, iomsg=iomsg) "Reference sound speed: ", self%reference_state(4), new_line('a')

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) 'Angles:  Theta_ib [deg]     Theta_ie [deg]'//new_line('a')
    do i = 1, 4
write(unit, '(a, i0, 2(a,2(f6.2, 1x)), a)', iostat=iostat, iomsg=iomsg) 'Cell: ', i, ' [ ', self%theta_ib(i, :) * (180.0_rk / pi), &
        '], [ ', self%theta_ie(i, :) * (180.0_rk / pi), ']'//new_line('a')
    end do

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) 'Cell Conserved Vars state [rho,u,v,p]'//new_line('a')
    do i = 1, 4
      write(unit, '(a, i0, a, 4(f6.2, 1x), a)', iostat=iostat, iomsg=iomsg) 'Cell: ', i, ' [ ', self%cell_conserved_vars(:,1,i) ,'] (arc 1)' // new_line('a')
      write(unit, '(a, 4(f6.2, 1x), a)', iostat=iostat, iomsg=iomsg) '        [ ', self%cell_conserved_vars(:, 2, i), '] (arc 2)'//new_line('a')
    end do

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) 'Reconstructed state [rho,u,v,p]'//new_line('a')
    do i = 1, 4
      write(unit, '(a, i0, a, 4(f6.2, 1x), a)', iostat=iostat, iomsg=iomsg) &
        'Cell: ', i, ' [ ', self%reconstructed_state(:, i), ']'//new_line('a')
    end do

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) "P' in cell?"//new_line('a')
    do i = 1, 4
      write(unit, '(a, i0, a, l2, a)', iostat=iostat, iomsg=iomsg) 'Cell: ', i, ' [', self%p_prime_in_cell(i), ' ]'//new_line('a')
    end do

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    ! write(unit, '(a, (f0.4))', iostat=iostat) 'Theta_ib', self%theta_ib
    ! write(unit, '((a,f0.4))', iostat=iostat) 'Theta_ie', self%theta_ie
    ! write(unit, '((a,f0.4))', iostat=iostat) 'Theta_ie'
  end subroutine write_cone

  logical pure function determine_if_p_prime_is_in_cell(edge_vectors, p_prime_vector) result(in_cell)
    !< Implementation of whether the P' point is inside the current cell/control volume. This
    !< uses the cross product of 2 vectors in 2d, which gives a scalar
    type(vector_t), intent(in) :: p_prime_vector
    real(rk), dimension(2, 2, 2), intent(in) :: edge_vectors

    type(vector_t) :: edge_vector_1, edge_vector_2

    in_cell = .false.
    edge_vector_1 = vector_t(x=[edge_vectors(1, 1, 1), edge_vectors(1, 2, 1)], &
                             y=[edge_vectors(2, 1, 1), edge_vectors(2, 2, 1)])

    edge_vector_2 = vector_t(x=[edge_vectors(1, 1, 2), edge_vectors(1, 2, 2)], &
                             y=[edge_vectors(2, 1, 2), edge_vectors(2, 2, 2)])

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

  pure subroutine get_arc_segments(lines, origin_in_cell, circle_xy, circle_radius, arc_segments, n_arcs)
    !< Given 2 lines and a circle, find their intersections and starting/ending angles for each arc

    real(rk), dimension(2, 2, 2), intent(in) :: lines
    !< ((x,y), (tail,head), (line_1:line_2)); set of vectors that define the lines to intersect with the circle

    logical, intent(in) :: origin_in_cell  !< is the circle origin in the cell defined by the two lines?
    real(rk), dimension(2), intent(in) :: circle_xy       !< (x,y) origin of the circle
    real(rk), intent(in) :: circle_radius                 !< radius of the circle
    real(rk), dimension(2, 2), intent(out):: arc_segments !< ((start,end), (arc1, arc2))
    integer(ik), intent(out) :: n_arcs !< Number of arcs with a valid start and end angle

    integer(ik) :: n_intersections_per_line !< # intersections for a single line
    integer(ik), dimension(2) :: n_intersections !< # intersections for all lines
    real(rk), dimension(2, 2) :: line
    real(rk), dimension(2) :: intersection_angles_per_line !< intersection angles
    real(rk), dimension(2, 2) :: intersection_angles !< intersection angles for all lines
    integer(ik) :: i

    n_intersections_per_line = 0
    n_intersections = 0
    line = 0.0_rk
    intersection_angles_per_line = 0.0_rk
    intersection_angles = 0.0_rk

    do i = 1, 2
      line = lines(:, :, i)
      call get_intersection_angles(line_xy=line, circle_xy=circle_xy, circle_radius=circle_radius, &
                                   arc_angles=intersection_angles_per_line, n_intersections=n_intersections_per_line)
      intersection_angles(:, i) = intersection_angles_per_line
      n_intersections(i) = n_intersections_per_line
    end do
    call get_theta_start_end(thetas=intersection_angles, origin_in_cell=origin_in_cell, &
                             n_intersections=n_intersections, &
                             theta_start_end=arc_segments, n_arcs=n_arcs)
  end subroutine

  pure subroutine get_theta_start_end(thetas, origin_in_cell, n_intersections, theta_start_end, n_arcs)
    real(rk), dimension(2, 2), intent(in) :: thetas  !< ((intersection 1, intersection 2), (line 1, line 2))
    logical, intent(in) :: origin_in_cell  !< is the circle center in the cell?
    integer(ik), dimension(2), intent(in) :: n_intersections
    real(rk), dimension(2, 2), intent(out) :: theta_start_end  !< ((start, end), (arc1, arc2))
    integer(ik), intent(out) :: n_arcs

    n_arcs = 0

    associate(theta_ib=>theta_start_end(1, :), theta_ie=>theta_start_end(2, :), &
              n1=>n_intersections(1), n2=>n_intersections(2))

      ! Default to zero contribution
      n_arcs = 0
      theta_ib = 0.0_rk
      theta_ie = 0.0_rk

      if(n1 == 0 .and. n1 == 0 .and. origin_in_cell) then
        n_arcs = 0
        theta_ib = 0.0_rk
        theta_ie = 2.0_rk * pi

      else if(n1 == 0 .and. n2 == 2) then
        n_arcs = 1
        theta_ib = thetas(2, 1)
        theta_ie = thetas(2, 2)

      else if(n1 == 1 .and. n2 == 1) then
        n_arcs = 1
        theta_ib = thetas(1, 1)
        theta_ie = thetas(2, 1)

      else if(n1 == 2 .and. n1 == 0) then
        n_arcs = 1
        theta_ib = thetas(1, 2)
        theta_ie = thetas(1, 1)

      else if(n1 == 2 .and. n1 == 2) then
        n_arcs = 2
        theta_ib(1) = thetas(1, 2)
        theta_ie(1) = thetas(2, 2)

        theta_ib(2) = thetas(2, 1)
        theta_ie(2) = thetas(1, 1)

      end if
    end associate
  end subroutine

  pure subroutine get_intersection_angles(line_xy, circle_xy, circle_radius, arc_angles, n_intersections)
    !< Given a arbitrary line from (x1,y1) to (x2,y2) and a circle at (x,y) with a given radius, find
    !< the angle that a the vector from the circle's center to the intersection point(s) has with respect
    !< to the x-axis

    real(rk), dimension(2, 2), intent(in) :: line_xy !< ((x,y), (point_1, point_2))
    real(rk), dimension(2), intent(in) :: circle_xy !< (x,y)
    real(rk), intent(in) :: circle_radius
    integer(ik), intent(out) :: n_intersections
    real(rk), dimension(2), intent(out):: arc_angles

    real(rk), dimension(2, 2) :: intersection_xy !< ((x,y), (point_1, point_2))
    logical, dimension(2) :: valid_intersection !< (point_1, point_2)

    intersection_xy = 0.0_rk
    valid_intersection = .false.

    call vector_circle_intersect(line_xy, circle_xy, circle_radius, intersection_xy, valid_intersection)
    n_intersections = count(valid_intersection)
    arc_angles = intersection_angle_from_x_axis(circle_xy, intersection_xy)
  end subroutine

  pure subroutine vector_circle_intersect(line_xy, circle_xy, circle_radius, intersection_xy, valid_intersection)
    !< Find the intersections between an arbitrary line and circle. There can be 0, 1, or 2 intersections.

    real(rk), dimension(2, 2), intent(in) :: line_xy !< ((x,y), (point_1, point_2))
    real(rk), dimension(2), intent(in) :: circle_xy !< (x,y)
    real(rk), intent(in) :: circle_radius
    real(rk), dimension(2, 2), intent(out) :: intersection_xy !< ((x,y), (point_1, point_2))
    logical, dimension(2), intent(out) :: valid_intersection !< (point_1, point_2)

    real(rk) :: discriminiant  !< term under the square root in the quadratic formula
    integer(ik) :: i
    real(rk), dimension(2) :: t !< scale factor (should be between 0 and 1) of where the intersection point is along the line
    real(rk) :: a, b, c  !< quadratic formula variables

    t = 0.0_rk
    discriminiant = 0.0_rk
    a = 0.0_rk
    b = 0.0_rk
    c = 0.0_rk

    associate(x0=>line_xy(1, 1), x1=>line_xy(1, 2), y0=>line_xy(2, 1), y1=>line_xy(2, 2), &
              r=>circle_radius, h=>circle_xy(1), k=>circle_xy(2))
      a = (x1 - x0)**2 + (y1 - y0)**2
      b = 2 * (x1 - x0) * (x0 - h) + 2 * (y1 - y0) * (y0 - k)
      c = (x0 - h)**2 + (y0 - k)**2 - r**2
    end associate

    discriminiant = b**2 - 4 * a * c

    if(discriminiant > 0.0_rk) then
      ! This used the alternative quadratic formula better suited for floating point operations
      if(near_zero(-b + sqrt(discriminiant))) then
        t(1) = 0.0_rk
      else
        t(1) = (2 * c) / (-b + sqrt(discriminiant))
      end if

      if(near_zero(-b - sqrt(discriminiant))) then
        t(2) = 1.0_rk
      else
        t(2) = (2 * c) / (-b - sqrt(discriminiant))
      end if

    end if

    ! The scale factor t must be between 0 and 1, otherwise there is no intersection
    valid_intersection = .false.
    do i = 1, 2
      if(t(i) >= 0.0_rk .and. t(i) <= 1.0_rk) then
        valid_intersection(i) = .true.
      end if
    end do

    intersection_xy = 0.0_rk
    do i = 1, 2
      associate(x0=>line_xy(1, 1), x1=>line_xy(1, 2), y0=>line_xy(2, 1), y1=>line_xy(2, 2), &
                x_t=>intersection_xy(1, i), y_t=>intersection_xy(2, i))
        x_t = (x1 - x0) * t(i) + x0
        y_t = (y1 - y0) * t(i) + y0
      end associate
    end do
  end subroutine

  pure function intersection_angle_from_x_axis(circle_origin_xy, intersection_xy) result(angles)
    !< Given the the intersection point and origin of the circle it intersected, determine the
    !< angle with respect to the x axis from 0 to 2pi

    real(rk), dimension(2) :: angles
    real(rk), dimension(2), intent(in) :: circle_origin_xy !< (x,y)
    real(rk), dimension(2, 2), intent(in) :: intersection_xy !< ((x,y), (point_1, point_2))

    integer(ik) :: i
    real(rk) :: dx1, dy1
    real(rk) :: dx2, dy2
    real(rk) :: test

    test = 0.0_rk
    angles = 0.0_rk
    dx1 = 0.0_rk
    dy1 = 0.0_rk
    dx2 = 0.0_rk
    dy2 = 0.0_rk

    dx1 = intersection_xy(1, 1) - circle_origin_xy(1)
    dy1 = intersection_xy(2, 1) - circle_origin_xy(2)

    dx2 = intersection_xy(1, 2) - circle_origin_xy(1)
    dy2 = intersection_xy(2, 2) - circle_origin_xy(2)

    angles(1) = atan2(dy1, dx1)
    angles(2) = atan2(dy2, dx2)

    do i = 1, 2
      if(angles(i) < 0.0_rk) then
        angles(i) = angles(i) + 2.0_rk * pi
      end if
    end do
  end function

end module mod_cone
