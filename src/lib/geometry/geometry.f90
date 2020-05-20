module mod_geometry

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use math_constants, only: pi, rad2deg
  use mod_floating_point_utils, only: near_zero, equal, EPS
  use mod_intersections, only: intersect

  implicit none

  ! private
  ! public :: find_line_circle_intersections, intersection_angle_from_x_axis, &
  !           super_circle, circle_inside_circle, get_arc_segments

contains

  subroutine find_midpoint_arc_segments(edge_vector, origin_in_cell, &
                                        circle_xy, circle_radius, arc_segments, n_arcs_in_cell)
    !< Since the midpoint mach cone is always on a single straight edge, we can avoid
    !< some unnecessary duplicate work.

    real(rk), dimension(2, 2), intent(in) :: edge_vector
    !< ((x,y), (point 1, point 2)); position of the line points

    logical, dimension(2), intent(in) :: origin_in_cell
    !< (cell 1, cell 2); is the mach cone origin in this cell?

    real(rk), dimension(2), intent(in) :: circle_xy
    !< (x,y); position of the mach cone origin

    real(rk), intent(in) :: circle_radius
    !< radius of the mach cone

    real(rk), dimension(2, 2), intent(in) :: arc_segments
    !< ((theta begin, theta end), (cell 1, cell 2))

    integer(ik), dimension(2), intent(out) :: n_arcs_in_cell
    !< (cell 1, cell 2); number of arcs contained in each cell

  end subroutine find_midpoint_arc_segments

  subroutine find_arc_segments(origin, vector_1, vector_2, origin_in_cell, &
                               circle_xy, circle_radius, arc_segments, n_arcs)
    !< Given 2 lines and a circle, find their intersections and starting/ending angles for each arc

    real(rk), dimension(2), intent(in) :: origin   !< (x,y) origin of both vectors
    real(rk), dimension(2), intent(in) :: vector_1 !< (x,y) head of the 1st vector
    real(rk), dimension(2), intent(in) :: vector_2 !< (x,y) head of the 2nd vector

    logical, intent(in) :: origin_in_cell  !< is the circle origin in the cell defined by the two lines?
    real(rk), dimension(2), intent(in) :: circle_xy       !< (x,y) origin of the circle
    real(rk), intent(in) :: circle_radius                 !< radius of the circle

    ! Output
    real(rk), dimension(2, 2), intent(out):: arc_segments !< ((start,end), (arc1, arc2))
    integer(ik), intent(out) :: n_arcs !< Number of arcs with a valid start and end angle

    ! Dummy
    integer(ik) :: n_intersections_per_line !< # intersections for a single line
    integer(ik), dimension(2) :: n_intersections !< (vector 1, vector 2); # intersections for all lines
    real(rk), dimension(2, 2) :: line !< ((x,y), (tail,head)); Single line point locations
    real(rk), dimension(2) :: intersection_angles_per_line !< intersection angles
    real(rk), dimension(2, 2) :: intersection_angles !< intersection angles for all lines
    integer(ik) :: i, n_int
    logical, dimension(2) :: valid_intersections
    logical, dimension(2, 2) :: total_valid_intersections !< ((intersection 1, intersection 2), (line 1, line 2)); valid intersections for each line

    valid_intersections = .false.
    total_valid_intersections = .false.
    n_intersections_per_line = 0
    n_intersections = 0
    intersection_angles_per_line = 0.0_rk
    intersection_angles = 0.0_rk

    line(:, 1) = origin ! both lines share the same origin, so only do this once

    ! Line 1
    line(:, 2) = vector_1
    ! For a given line & circle intersection, find the angle with respect to the x-axis for each intersection point
    call get_intersection_angles(line_xy=line, circle_xy=circle_xy, circle_radius=circle_radius, &
                                 arc_angles=intersection_angles_per_line, &
                                 valid_intersections=valid_intersections, &
                                 n_intersections=n_int)
    intersection_angles(1, :) = intersection_angles_per_line
    n_intersections(1) = n_int
    total_valid_intersections(:, 1) = valid_intersections

    ! Line 2
    line(:, 2) = vector_2
    ! For a given line & circle intersection, find the angle with respect to the x-axis for each intersection point
    call get_intersection_angles(line_xy=line, circle_xy=circle_xy, circle_radius=circle_radius, &
                                 arc_angles=intersection_angles_per_line, &
                                 valid_intersections=valid_intersections, &
                                 n_intersections=n_int)
    intersection_angles(2, :) = intersection_angles_per_line
    n_intersections(2) = n_int
    total_valid_intersections(:, 2) = valid_intersections

    ! Find arc starting and ending angles
    call get_theta_start_end(thetas=intersection_angles, origin_in_cell=origin_in_cell, &
                             valid_intersections=total_valid_intersections, &
                             n_intersections=n_intersections, &
                             theta_start_end=arc_segments, n_arcs=n_arcs)
  end subroutine find_arc_segments

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
        n_arcs = 2

        theta_ib(1) = thetas(2, 1)
        theta_ie(1) = thetas(1, 1)

        if(theta_ie(1) < theta_ib(1)) then
          if(theta_ie(1) < 0.0_rk) then
            theta_ie(1) = theta_ie(1) + 2 * pi
          end if
        end if

        theta_ib(2) = thetas(1, 2)
        theta_ie(2) = thetas(2, 2)

        if(theta_ie(2) < theta_ib(2)) then
          if(theta_ie(2) < 0.0_rk) then
            theta_ie(2) = theta_ie(2) + 2 * pi
          end if
        end if

      end if ! if there are intersections

    end associate
  end subroutine get_theta_start_end

  subroutine get_intersection_angles(line_xy, circle_xy, circle_radius, arc_angles, valid_intersections, n_intersections)
    !< Given a arbitrary line from (x1,y1) to (x2,y2) and a circle at (x,y) with a given radius, find
    !< the angle that a the vector from the circle's center to the intersection point(s) has with respect
    !< to the x-axis

    ! Input
    real(rk), dimension(2, 2), intent(in) :: line_xy !< ((x,y), (point_1, point_2))
    real(rk), dimension(2), intent(in) :: circle_xy !< (x,y)
    real(rk), intent(in) :: circle_radius

    ! Output
    integer(ik), intent(out) :: n_intersections
    logical, dimension(2), intent(out) :: valid_intersections !< (point_1, point_2); .true. or .false.
    real(rk), dimension(2), intent(out):: arc_angles

    real(rk), dimension(2) :: arc_angles_new
    logical, dimension(2) :: valid_intersections_new

    ! Dummy
    integer(ik) :: i, j

    real(rk), dimension(2, 2) :: intersection_xy !< ((x,y), (point_1, point_2))
    real(rk), dimension(2, 2) :: intersection_xy_new !< ((x,y), (point_1, point_2))

    intersection_xy = 0.0_rk
    valid_intersections = .false.

    ! ! Find the intersection (x,y) locations (if any)
    ! call intersect(line_xy, circle_xy, circle_radius, intersection_xy_new, valid_intersections_new, arc_angles_new)
    call intersect(line_xy, circle_xy, circle_radius, intersection_xy, valid_intersections, arc_angles)
    ! call find_line_circle_intersections(line_xy, circle_xy, circle_radius, intersection_xy, valid_intersections)

    ! do i = 1, 2
    !   do j = 1, 2
    !     associate(new=>intersection_xy_new(i,j), old=>intersection_xy(i,j))
    !       if (abs(new-old) > 2.0_rk * epsilon(1.0_rk)) then
    !         print*, i, j, new, old, abs(new-old)
    !         print*, "new /= old"
    !       end if
    !     end associate
    !   end do
    ! end do

    ! arc_angles = 0.0_rk
    n_intersections = 0
    do i = 1, 2
      if(valid_intersections(i)) then
        ! arc_angles(i) = intersection_angle_from_x_axis(circle_xy, intersection_xy(:, i))
        n_intersections = n_intersections + 1
      end if
    end do

    ! do i = 1, 2
    !   associate(new=>arc_angles_new(i), old=>arc_angles(i))
    !     if (abs(new-old) > epsilon(1.0_rk)) then
    !       print*, "new /= old"
    !       print*, new, old, abs(new-old)
    !       print*, rad2deg(new), rad2deg(old)
    !       ! error stop
    !     end if
    !   end associate
    ! end do

    ! do i = 1, 2
    !   associate(new=>valid_intersections_new(i), old=>valid_intersections(i))
    !     if (new .neqv. old) then
    !       print*, 'Checking valid_intersections'
    !       print*, i, j
    !       print*, 'Orig: ', valid_intersections
    !       print*, 'New:  ', valid_intersections_new

    !       write(*,'(a, 4(es16.6))') 'Orig (x,y): ', intersection_xy
    !       write(*,'(a, 4(es16.6))') 'New  (x,y): ', intersection_xy_new
    !       write(*,'(a, 4(es16.6))') 'Orig (rad): ', arc_angles
    !       write(*,'(a, 4(es16.6))') 'New  (rad): ', arc_angles_new
    !       error stop "new /= old"
    !     end if
    !   end associate
    ! end do
  end subroutine get_intersection_angles

  subroutine find_line_circle_intersections(line_xy, circle_xy, circle_radius, intersection_xy, valid_intersection)
    !< Find the intersections between an arbitrary line and circle. There can be 0, 1, or 2 intersections.

    real(rk), dimension(2, 2), intent(in) :: line_xy !< ((x,y) (point_1, point_2)); Line location
    real(rk), dimension(2), intent(in) :: circle_xy !< (x,y); Circle origin
    real(rk), intent(in) :: circle_radius !< Radius of the circle (duh)
    real(rk), dimension(2, 2), intent(out) :: intersection_xy !< ((x,y), (point_1, point_2)); Intersection point
    logical, dimension(2), intent(out) :: valid_intersection !< (point_1, point_2); Is there an intersection or not?

    real(rk) :: discriminiant  !< term under the square root in the quadratic formula
    real(rk) :: sqrt_discriminiant  !< term under the square root in the quadratic formula
    integer(ik) :: i
    real(rk), dimension(2) :: t !< scale factor (should be between 0 and 1) of where the intersection point is along the line
    real(rk) :: a, b, c  !< quadratic formula variables
    real(rk) :: x0, x1, dx
    real(rk) :: y0, y1, dy
    real(rk) :: root_1, root_2

    t = 0.0_rk
    discriminiant = 0.0_rk
    a = 0.0_rk
    b = 0.0_rk
    c = 0.0_rk

    x0 = line_xy(1, 1)
    x1 = line_xy(1, 2)
    dx = x1 - x0

    y0 = line_xy(2, 1)
    y1 = line_xy(2, 2)
    dy = y1 - y0

    ! print *, "line_xy(:,1) ", line_xy(:, 1)
    ! print *, "line_xy(:,2) ", line_xy(:, 2)
    ! print *, "circle_xy    ", circle_xy
    ! print *, "circle_radius", circle_radius

    associate(r=>circle_radius, h=>circle_xy(1), k=>circle_xy(2))
      a = dx**2 + dy**2
      b = 2.0_rk * dx * (x0 - h) + 2.0_rk * dy * (y0 - k)
      c = (x0 - h)**2 + (y0 - k)**2 - r**2

      ! print *, 'a', a
      ! print *, 'b', b
      ! print *, 'c', c
      ! print *, 'r**2', r**2

      ! print *, 'dx', dx, 'dy', dy
      ! print *, "(x0 - h), (y0 - k) ",(x0 - h),(y0 - k)
    end associate

    discriminiant = b**2 - 4.0_rk * a * c
    ! write(*, '(3(es16.6))') discriminiant, b**2, 4.0_rk * a * c

    if(discriminiant > 0.0_rk) then
      sqrt_discriminiant = sqrt(discriminiant)
      root_1 = -b + sqrt_discriminiant
      root_2 = -b - sqrt_discriminiant

      ! print *, 'root_1', root_1
      ! print *, 'root_2', root_2
      ! This used the alternative quadratic formula better suited for floating point operations
      if(abs(root_1) < EPS) then
        t(1) = 0.0_rk ! t_1 -> 0 when intersection is at the vector start point
      else
        t(1) = (2.0_rk * c) / root_1
      end if

      if(abs(root_2) < EPS) then
        t(2) = 1.0_rk  ! t_2 -> 1 when intersection is at the vector end point
      else
        t(2) = (2.0_rk * c) / root_2
      end if
    end if

    ! The scale factor t must be between 0 and 1, otherwise there is no intersection
    valid_intersection = .false.
    if(discriminiant > 0.0_rk) then
      do i = 1, 2 ! TODO: Can we use a where block?
        if(t(i) < 0.0_rk .or. t(i) > 1.0_rk) then
          valid_intersection(i) = .false.
        else
          valid_intersection(i) = .true.
        end if
      end do
    end if

    intersection_xy = 0.0_rk
    if(discriminiant > 0.0_rk) then
      do i = 1, 2
        ! x_t
        intersection_xy(1, i) = dx * t(i) + x0
        ! y_t
        intersection_xy(2, i) = dy * t(i) + y0
      end do
    end if
  end subroutine find_line_circle_intersections

  pure real(rk) function intersection_angle_from_x_axis(circle_origin_xy, intersection_xy) result(angle)
    !< Given the the intersection point and origin of the circle it intersected, determine the
    !< angle with respect to the x axis from 0 to 2pi

    real(rk), dimension(2), intent(in) :: circle_origin_xy !< (x,y)
    real(rk), dimension(2), intent(in) :: intersection_xy !< ((x,y)

    angle = atan2(y=intersection_xy(2) - circle_origin_xy(2), &
                  x=intersection_xy(1) - circle_origin_xy(1))
  end function intersection_angle_from_x_axis

  pure subroutine super_circle(origins, radii, new_origin, new_radius)
    !< Given 2 circles with origins and radii, find a new circle that will
    !< superscribe both of them (e.g. contain both within it)
    !< See http://dx.doi.org/10.4208/cicp.220210.020710a Section 3.4 (Entropy Fix) for more details

    real(rk), dimension(2, 2), intent(in) :: origins   !< ((x,y), (circle1, circle2)); circle origin
    real(rk), dimension(2), intent(in) :: radii        !< (circle1, circle2); radius of each circle
    real(rk), dimension(2), intent(out) :: new_origin  !< (x,y); center for the new circle
    real(rk), intent(out) :: new_radius                 !< radius of the new circle (duh)

    logical :: is_inside
    integer(ik) :: largest_circle_idx
    real(rk) :: d

    is_inside = circle_inside_circle(origins=origins, radii=radii)
    largest_circle_idx = maxloc(radii, dim=1)

    if(is_inside) then
      new_origin = origins(:, largest_circle_idx)
      new_radius = radii(largest_circle_idx)
    else
      associate(x2=>origins(:, 2), x1=>origins(:, 1), &
                r2=>radii(2), r1=>radii(1))
        d = norm2(x2 - x1)
        new_radius = ((r1 + r2 + d) / 2.0_rk)
        new_origin = x1 + (new_radius - r1) * (x2 - x1) / d
      end associate
    end if

    ! if(abs(new_origin(1)) < TINY_DIST) new_origin(1) = 0.0_rk
    ! if(abs(new_origin(2)) < TINY_DIST) new_origin(2) = 0.0_rk
  end subroutine super_circle

  logical pure function circle_inside_circle(origins, radii) result(is_inside)
    !< Given 2 circles with origins and radii, determine if the 2nd circle is
    !< fully inside the 1st
    real(rk), dimension(2, 2), intent(in) :: origins  !< ((x,y), (circle1, circle2)); circle origin
    real(rk), dimension(2), intent(in) :: radii       !< (circle1, circle2); radius of each circle

    real(rk) :: distance_squared
    integer(ik) :: larger_circle, smaller_circle

    larger_circle = maxloc(radii, dim=1)
    if(larger_circle == 1) then
      smaller_circle = 2
    else
      smaller_circle = 1
    end if

    associate(x1=>origins(1, larger_circle), y1=>origins(2, larger_circle), &
              x2=>origins(1, smaller_circle), y2=>origins(2, smaller_circle), &
              r2=>radii(smaller_circle), r1=>radii(larger_circle))

      if(equal(r1, r2) .and. all(equal(origins(:, 1), origins(:, 2)))) then
        ! If the origin and radius are the same (to w/in a small tolerance)
        is_inside = .true.
      else
        distance_squared = sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)))

        if(distance_squared + r2 <= r1) then
          is_inside = .true.
        else
          is_inside = .false.
        end if
      end if
    end associate
  end function circle_inside_circle

end module mod_geometry
