module mod_cone
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use math_constants, only: pi, rad2deg
  use mod_floating_point_utils, only: near_zero, equal
  use mod_eos, only: eos
  use mod_vector, only: vector_t, operator(.unitnorm.), operator(.dot.), operator(.cross.), angle_between

  implicit none

  type :: cone_t
    !< Type to encapsulate and calculate the geometry of the mach cone for the evolution operators

    real(rk), dimension(2, 4) :: theta_ie = 0.0_rk
    !< ((arc_1, arc_2), (cell_1:cell_4)); arc end angle

    real(rk), dimension(2, 4) :: theta_ib = 0.0_rk
    !< ((arc_1, arc_2), (cell_1:cell_4)); arc begin angle

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

  type(cone_t) function new_cone(tau, edge_vectors, reconstructed_state, reference_state, cell_indices)
    !< Constructor for the Mach cone type
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
    integer(ik) :: n_arcs, arc
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
              p=>reference_state(4))
      if(rho < 0.0) then
        print *, reference_state
        print *, edge_vectors(:, :, 1)
        print *, edge_vectors(:, :, 2)
        ! print*, edge_vectors(:,:,3)
        ! print*, edge_vectors(:,:,4)
      end if
      a = eos%calc_sound_speed(pressure=p, density=rho)
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
      ! debug_write(*,*) x, u_tilde, tau, y, v_tilde
    end associate

    ! The P' vector points from P (tail), to P' (head)
    ! debug_write(*,*) "P'", new_cone%p_prime_xy
    p_prime_vector = vector_t(x=[edge_vectors(1, 1, 1), new_cone%p_prime_xy(1)], &
                              y=[edge_vectors(2, 1, 1), new_cone%p_prime_xy(2)])

    ! Loop through each neighbor cell and determine intersections and angles
    do neighbor_cell = 1, new_cone%n_neighbor_cells
      p_prime_in_cell = .false.

      cell_ij = cell_indices(:, neighbor_cell)  ! cell_indices is indexed via ((i,j), cell_1:cell_n)

      ! indexing for edge_vectors is ((x,y), (tail,head), (vector_1:vector_n))
      single_cell_edge_vectors(:, :, 1) = edge_vectors(:, :, edge_vector_ordering(1, neighbor_cell))
      single_cell_edge_vectors(:, :, 2) = edge_vectors(:, :, edge_vector_ordering(2, neighbor_cell))

      ! P' can only be in 1 cell
      p_prime_in_cell = determine_if_p_prime_is_in_cell(single_cell_edge_vectors, p_prime_vector)
      if(.not. any(new_cone%p_prime_in_cell)) then
        new_cone%p_prime_in_cell(neighbor_cell) = p_prime_in_cell
        if(p_prime_in_cell) new_cone%p_prime_ij = cell_ij
      end if

      call get_arc_segments(lines=single_cell_edge_vectors, origin_in_cell=p_prime_in_cell, &
                            circle_xy=new_cone%p_prime_xy, circle_radius=new_cone%radius, &
                            arc_segments=theta_ib_ie, n_arcs=n_arcs)

      if(n_arcs > 0) then
        do arc = 1, n_arcs
          new_cone%cell_conserved_vars(:, arc, neighbor_cell) = reconstructed_state(:, neighbor_cell)
        end do
      end if
      ! If all the arcs already add up to 2pi, then skip the remaining neighbor cells
      if(equal(sum(new_cone%theta_ie - new_cone%theta_ib), 2.0_rk * pi, epsilon=1e-10_rk)) exit
      ! //TODO: test the above statement in unit testing
      new_cone%n_arcs(neighbor_cell) = n_arcs
      new_cone%theta_ib(:, neighbor_cell) = theta_ib_ie(1, :)
      new_cone%theta_ie(:, neighbor_cell) = theta_ib_ie(2, :)

    end do

    if(count(new_cone%n_arcs >= 2) > 1) then
      error stop "Too many arcs in the mach cone (count(cone%n_arcs >= 2) > 1)"
    end if

    if(.not. equal(sum(new_cone%theta_ie - new_cone%theta_ib), 2.0_rk * pi, epsilon=1e-10_rk)) then
      ! print *, 'sum(new_cone%theta_ie - new_cone%theta_ib)', sum(new_cone%theta_ie - new_cone%theta_ib) - 2.0_rk * pi
      write(*, *) new_cone
      ! debug_write(*, *) new_cone
      error stop "Cone arcs do not add up to 2pi"
    end if
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
    write(unit, '(a, 2(es10.3,1x), a)', iostat=iostat, iomsg=iomsg) "P (x,y): (", self%p_xy, ")"//new_line('a')
    write(unit, '(a, 2(es10.3,1x), a)', iostat=iostat, iomsg=iomsg) "P'(x,y): (", self%p_prime_xy, ")"//new_line('a')
    write(unit, '(a, 2(es10.3,1x), a)', iostat=iostat, iomsg=iomsg) &
      "P'(x,y) - P(x,y): (", self%p_prime_xy - self%p_xy, ")"//new_line('a')

    write(unit, '(a, 2(i7,1x), a)', iostat=iostat, iomsg=iomsg) "P'(i,j): (", self%p_prime_ij, ")"//new_line('a')

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a, f0.3, a)', iostat=iostat, iomsg=iomsg) "Reference density:     ", self%reference_state(1), new_line('a')
    write(unit, '(a, f0.3, a)', iostat=iostat, iomsg=iomsg) "Reference x velocity:  ", self%reference_state(2), new_line('a')
    write(unit, '(a, f0.3, a)', iostat=iostat, iomsg=iomsg) "Reference y velocity:  ", self%reference_state(3), new_line('a')
    write(unit, '(a, f0.3, a)', iostat=iostat, iomsg=iomsg) "Reference sound speed: ", self%reference_state(4), new_line('a')

    write(unit, '(a)', iostat=iostat) new_line('a')
    write(unit, '(a, 4(i0, 1x),a)', iostat=iostat) '# of valid arcs in each cell: [', self%n_arcs, ']'
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) 'Angles:  Theta_ib [deg]     Theta_ie [deg]'//new_line('a')
    do i = 1, 4
      write(unit, '(a, i0, 2(a,2(f7.2, 1x)), 2(a,2(f7.2, 1x)), a)', iostat=iostat, iomsg=iomsg) &
        'Cell: ', i, ' [ ', rad2deg(self%theta_ib(:, i)), &
        '], [ ', rad2deg(self%theta_ie(:, i)), '], delta theta = ', &
        rad2deg(self%theta_ie(:, i) - self%theta_ib(:, i)), ' '//new_line('a')
    end do
    write(unit, '(a, f0.2, a)', iostat=iostat, iomsg=iomsg) 'sum(arc delta theta) = ', &
      rad2deg(sum(self%theta_ie - self%theta_ib)), ' '//new_line('a')

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) 'Angles:  Theta_ib [rad]     Theta_ie [rad]'//new_line('a')
    do i = 1, 4
      write(unit, '(a, i0, 2(a,2(f7.2, 1x)), 2(a,2(f7.2, 1x)), a)', iostat=iostat, iomsg=iomsg) &
        'Cell: ', i, ' [ ', self%theta_ib(:, i), &
        '], [ ', self%theta_ie(:, i), '], delta theta = ', self%theta_ie(:, i) - self%theta_ib(:, i), ' '//new_line('a')
    end do
    write(unit, '(a, f0.2, a)', iostat=iostat, iomsg=iomsg) &
      'sum(arc delta theta) = ', sum(self%theta_ie - self%theta_ib), ' '//new_line('a')

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) 'Cell Conserved Vars state [rho,u,v,p]'//new_line('a')
    do i = 1, 4
      write(unit, '(a, i0, a, 4(es10.3, 1x), a)', iostat=iostat, iomsg=iomsg) &
        'Cell: ', i, ' [ ', self%cell_conserved_vars(:, 1, i), '] (arc 1)'//new_line('a')
      write(unit, '(a, 4(es10.3, 1x), a)', iostat=iostat, iomsg=iomsg) &
        '        [ ', self%cell_conserved_vars(:, 2, i), '] (arc 2)'//new_line('a')
    end do

    ! write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    ! write(unit, '(a)', iostat=iostat, iomsg=iomsg) 'Reconstructed state [rho,u,v,p]'//new_line('a')
    ! do i = 1, 4
    !   write(unit, '(a, i0, a, 4(f6.2, 1x), a)', iostat=iostat, iomsg=iomsg) &
    !     'Cell: ', i, ' [ ', self%reconstructed_state(:, i), ']'//new_line('a')
    ! end do

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

  pure logical function determine_if_p_prime_is_in_cell(edge_vectors, p_prime_vector) result(in_cell)
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

    ! Input
    real(rk), dimension(2, 2, 2), intent(in) :: lines
    !< ((x,y), (tail,head), (line_1:line_2)); set of vectors that define the lines to intersect with the circle
    logical, intent(in) :: origin_in_cell  !< is the circle origin in the cell defined by the two lines?
    real(rk), dimension(2), intent(in) :: circle_xy       !< (x,y) origin of the circle
    real(rk), intent(in) :: circle_radius                 !< radius of the circle

    ! Output
    real(rk), dimension(2, 2), intent(out):: arc_segments !< ((start,end), (arc1, arc2))
    integer(ik), intent(out) :: n_arcs !< Number of arcs with a valid start and end angle

    ! Dummy
    integer(ik) :: n_intersections_per_line !< # intersections for a single line
    integer(ik), dimension(2) :: n_intersections !< # intersections for all lines
    real(rk), dimension(2, 2) :: line !< ((x,y), (tail,head)); Single line point locations
    real(rk), dimension(2) :: intersection_angles_per_line !< intersection angles
    real(rk), dimension(2, 2) :: intersection_angles !< intersection angles for all lines
    integer(ik) :: i
    logical, dimension(2) :: valid_intersections
    logical, dimension(2, 2) :: total_valid_intersections !< ((intersection 1, intersection 2), (line 1, line 2)); valid intersections for each line

    valid_intersections = .false.
    total_valid_intersections = .false.
    n_intersections_per_line = 0
    n_intersections = 0
    line = 0.0_rk
    intersection_angles_per_line = 0.0_rk
    intersection_angles = 0.0_rk

    do i = 1, 2
      line = lines(:, :, i)

      ! For a given line & circle intersection, find the angle with respect to the x-axis for each
      ! intersection point
      call get_intersection_angles(line_xy=line, circle_xy=circle_xy, circle_radius=circle_radius, &
                                   arc_angles=intersection_angles_per_line, valid_intersections=valid_intersections)
      intersection_angles(i, :) = intersection_angles_per_line
      n_intersections(i) = count(valid_intersections)
      total_valid_intersections(:, i) = valid_intersections
    end do
    ! print*, 'total_valid_intersections', total_valid_intersections
    ! print*, 'intersection_angles', rad2deg(intersection_angles)
    ! Find arc starting and ending angles
    call get_theta_start_end(thetas=intersection_angles, origin_in_cell=origin_in_cell, &
                             valid_intersections=total_valid_intersections, &
                             n_intersections=n_intersections, &
                             theta_start_end=arc_segments, n_arcs=n_arcs)
    ! print*, 'theta_start_end', rad2deg(arc_segments)
  end subroutine

  pure subroutine get_theta_start_end(thetas, origin_in_cell, valid_intersections, n_intersections, theta_start_end, n_arcs)
    !< Find the starting and ending angle of the arc

    ! Input
    real(rk), dimension(2, 2), intent(in) :: thetas  !< ((intersection 1, intersection 2), (line 1, line 2))
    logical, intent(in) :: origin_in_cell  !< is the circle center in the cell?
    logical, dimension(2, 2), intent(in) :: valid_intersections
    integer(ik), dimension(2), intent(in) :: n_intersections

    ! Output
    real(rk), dimension(2, 2), intent(out) :: theta_start_end  !< ((start, end), (arc1, arc2))
    integer(ik), intent(out) :: n_arcs

    n_arcs = 0
    theta_start_end = 0.0_rk

    associate(theta_ib=>theta_start_end(1, :), theta_ie=>theta_start_end(2, :), &
              n1=>n_intersections(1), n2=>n_intersections(2))

      if(origin_in_cell) then
        if(n1 == 0 .and. n2 == 0) then
          n_arcs = 1
          theta_ib(1) = 0.0_rk
          theta_ie(1) = 2.0_rk * pi
        end if
      else ! origin not in cell
        if(n1 == 0 .and. n2 == 0) then
          n_arcs = 0
          theta_ib = 0.0_rk
          theta_ie = 0.0_rk
        end if
      end if

      if(n1 == 0 .and. n2 == 2) then
        n_arcs = 1
        theta_ib(1) = thetas(2, 1)
        theta_ie(1) = thetas(2, 2)
        if(thetas(2, 2) < thetas(2, 1)) theta_ie(1) = theta_ie(1) + 2.0_rk * pi

      else if(n1 == 2 .and. n2 == 0) then
        n_arcs = 1
        theta_ib(1) = thetas(1, 2)
        theta_ie(1) = thetas(1, 1)
        if(thetas(1, 1) < thetas(1, 2)) theta_ie(1) = theta_ie(1) + 2.0_rk * pi

      else if(n1 == 1 .and. n2 == 1) then
        n_arcs = 1
        theta_ib(1) = thetas(1, 2)
        theta_ie(1) = thetas(2, 2)
        if(thetas(2, 2) < thetas(1, 2)) theta_ie(1) = theta_ie(1) + 2.0_rk * pi

      else if(n1 == 2 .and. n2 == 2) then
        ! print*, rad2deg(thetas(:,1))
        ! print*, rad2deg(thetas(:,2))
        n_arcs = 2

        theta_ib(1) = thetas(2, 1)
        theta_ie(1) = thetas(1, 1)

        if(theta_ie(1) < theta_ib(1)) then
          if(theta_ie(1) < 0.0_rk) then
            theta_ie(1) = theta_ie(1) + 2 * pi
          end if
        end if

        ! if(theta_ie(1) < theta_ib(1)) then
        !   print*, 'a'
        !   if(theta_ie(1) < 0.0_rk) then
        !     theta_ib(1) = theta_ie(1)
        !     theta_ie(1) = theta_ie(1) + 2*pi
        !   end if
        ! end if

        ! if(thetas(2, 1) < thetas(1, 1)) then
        !   ! if(thetas(2, 1) < 0.0_rk .or. thetas(1, 1) < 0.0_rk) then
        !     theta_ib(1) = thetas(2, 1)
        !     theta_ie(1) = thetas(1, 1)
        !   ! end if
        ! end if

        theta_ib(2) = thetas(1, 2)
        theta_ie(2) = thetas(2, 2)

        if(theta_ie(2) < theta_ib(2)) then
          if(theta_ie(2) < 0.0_rk) then
            theta_ie(2) = theta_ie(2) + 2 * pi
          end if
        end if

      end if ! if there are intersections

    end associate
  end subroutine

  pure subroutine get_intersection_angles(line_xy, circle_xy, circle_radius, arc_angles, valid_intersections)
    !< Given a arbitrary line from (x1,y1) to (x2,y2) and a circle at (x,y) with a given radius, find
    !< the angle that a the vector from the circle's center to the intersection point(s) has with respect
    !< to the x-axis

    ! Input
    real(rk), dimension(2, 2), intent(in) :: line_xy !< ((x,y), (point_1, point_2))
    real(rk), dimension(2), intent(in) :: circle_xy !< (x,y)
    real(rk), intent(in) :: circle_radius

    ! Output
    logical, dimension(2), intent(out) :: valid_intersections !< (point_1, point_2); .true. or .false.
    real(rk), dimension(2), intent(out):: arc_angles

    ! Dummy
    integer(ik) :: i
    real(rk), dimension(2, 2) :: intersection_xy !< ((x,y), (point_1, point_2))

    intersection_xy = 0.0_rk
    valid_intersections = .false.

    ! Find the intersection (x,y) locations (if any)
    call find_line_circle_intersections(line_xy, circle_xy, circle_radius, intersection_xy, valid_intersections)

    arc_angles = 0.0_rk
    do i = 1, 2
      if(valid_intersections(i)) then
        arc_angles(i) = intersection_angle_from_x_axis(circle_xy, intersection_xy(:, i))
      end if
    end do
  end subroutine

  pure subroutine find_line_circle_intersections(line_xy, circle_xy, circle_radius, intersection_xy, valid_intersection)
    !< Find the intersections between an arbitrary line and circle. There can be 0, 1, or 2 intersections.

    ! Input
    real(rk), dimension(2, 2), intent(in) :: line_xy !< ((x,y) (point_1, point_2)); Line location
    real(rk), dimension(2), intent(in) :: circle_xy !< (x,y); Circle origin
    real(rk), intent(in) :: circle_radius !< Radius of the circle (duh)

    ! Output
    real(rk), dimension(2, 2), intent(out) :: intersection_xy !< ((x,y), (point_1, point_2)); Intersection point
    logical, dimension(2), intent(out) :: valid_intersection !< (point_1, point_2); Is there an intersection or not?

    ! Dummy
    real(rk) :: discriminiant  !< term under the square root in the quadratic formula
    real(rk) :: sqrt_discriminiant  !< term under the square root in the quadratic formula
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
      sqrt_discriminiant = sqrt(discriminiant)

      ! This used the alternative quadratic formula better suited for floating point operations
      if(near_zero(-b + sqrt_discriminiant)) then
        t(1) = 0.0_rk ! t_1 -> 0 when intersection is at the vector start point
      else
        t(1) = (2 * c) / (-b + sqrt_discriminiant)
      end if

      if(near_zero(-b - sqrt_discriminiant)) then
        t(2) = 1.0_rk  ! t_2 -> 1 when intersection is at the vector end point
      else
        t(2) = (2 * c) / (-b - sqrt_discriminiant)
      end if
    end if

    ! The scale factor t must be between 0 and 1, otherwise there is no intersection
    valid_intersection = .false.
    if(discriminiant > 0.0_rk) then
      do i = 1, 2
        if(t(i) >= 0.0_rk .and. t(i) <= 1.0_rk) then
          valid_intersection(i) = .true.
        end if
      end do
    end if

    intersection_xy = 0.0_rk
    if(discriminiant > 0.0_rk) then
      do i = 1, 2
        associate(x0=>line_xy(1, 1), x1=>line_xy(1, 2), y0=>line_xy(2, 1), y1=>line_xy(2, 2), &
                  x_t=>intersection_xy(1, i), y_t=>intersection_xy(2, i))
          x_t = (x1 - x0) * t(i) + x0
          y_t = (y1 - y0) * t(i) + y0
        end associate
      end do
    end if
  end subroutine

  pure real(rk) function intersection_angle_from_x_axis(circle_origin_xy, intersection_xy) result(angle)
    !< Given the the intersection point and origin of the circle it intersected, determine the
    !< angle with respect to the x axis from 0 to 2pi

    real(rk), dimension(2), intent(in) :: circle_origin_xy !< (x,y)
    real(rk), dimension(2), intent(in) :: intersection_xy !< ((x,y)

    angle = atan2(y=intersection_xy(2) - circle_origin_xy(2), &
                  x=intersection_xy(1) - circle_origin_xy(1))

    ! if(angle < 0.0_rk) angle = angle + 2.0_rk * pi
  end function

  pure subroutine set_arc_segments(theta_ib, theta_ie)
    !< Sometimes, especially with cells that contain 2 arcs
    real(rk), dimension(4, 2), intent(inout) :: theta_ib
    real(rk), dimension(4, 2), intent(inout) :: theta_ie
  end subroutine
end module mod_cone
