module mod_geometry

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use math_constants, only: pi, rad2deg
  use mod_globals, only: TINY_DIST
  use mod_floating_point_utils, only: near_zero, equal

  implicit none

  ! private
  ! public :: find_line_circle_intersections, intersection_angle_from_x_axis, &
  !           super_circle, circle_inside_circle, get_arc_segments

contains

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

    ! Find arc starting and ending angles
    call get_theta_start_end(thetas=intersection_angles, origin_in_cell=origin_in_cell, &
                             valid_intersections=total_valid_intersections, &
                             n_intersections=n_intersections, &
                             theta_start_end=arc_segments, n_arcs=n_arcs)
    ! print*, 'theta_start_end', rad2deg(arc_segments)
  end subroutine get_arc_segments

pure subroutine get_arc_segments_new(origin, vec_1_head, vec_2_head, origin_in_cell, circle_xy, circle_radius, arc_segments, n_arcs)
    !< Given 2 lines and a circle, find their intersections and starting/ending angles for each arc

    ! Input
    ! real(rk), dimension(2, 2, 2), intent(in) :: lines
    !< ((x,y), (tail,head), (line_1:line_2)); set of vectors that define the lines to intersect with the circle

    real(rk), dimension(2), intent(in) :: origin     !< (x,y) origin of both vectors
    real(rk), dimension(2), intent(in) :: vec_1_head !< (x,y) head of the 1st vector
    real(rk), dimension(2), intent(in) :: vec_2_head !< (x,y) head of the 2nd vector

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
    intersection_angles_per_line = 0.0_rk
    intersection_angles = 0.0_rk

    line(:, 1) = origin ! both lines share the same origin, so only do this once

    ! Line 1
    line(:, 2) = vec_1_head
    ! For a given line & circle intersection, find the angle with respect to the x-axis for each intersection point
    call get_intersection_angles(line_xy=line, circle_xy=circle_xy, circle_radius=circle_radius, &
                                 arc_angles=intersection_angles_per_line, valid_intersections=valid_intersections)
    intersection_angles(1, :) = intersection_angles_per_line
    n_intersections(1) = count(valid_intersections)
    total_valid_intersections(:, 1) = valid_intersections

    ! Line 2
    line(:, 2) = vec_2_head
    ! For a given line & circle intersection, find the angle with respect to the x-axis for each intersection point
    call get_intersection_angles(line_xy=line, circle_xy=circle_xy, circle_radius=circle_radius, &
                                 arc_angles=intersection_angles_per_line, valid_intersections=valid_intersections)
    intersection_angles(2, :) = intersection_angles_per_line
    n_intersections(2) = count(valid_intersections)
    total_valid_intersections(:, 2) = valid_intersections

    ! Find arc starting and ending angles
    call get_theta_start_end(thetas=intersection_angles, origin_in_cell=origin_in_cell, &
                             valid_intersections=total_valid_intersections, &
                             n_intersections=n_intersections, &
                             theta_start_end=arc_segments, n_arcs=n_arcs)
  end subroutine get_arc_segments_new

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
  end subroutine get_intersection_angles

  pure subroutine find_line_circle_intersections(line_xy, circle_xy, circle_radius, intersection_xy, valid_intersection)
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
      associate(x0=>line_xy(1, 1), x1=>line_xy(1, 2), &
                y0=>line_xy(2, 1), y1=>line_xy(2, 2))
        ! x_t
        intersection_xy(1, :) = (x1 - x0) * t + x0

        ! y_t
        intersection_xy(2, :) = (y1 - y0) * t + y0
      end associate
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

    if(abs(new_origin(1)) < TINY_DIST) new_origin(1) = 0.0_rk
    if(abs(new_origin(2)) < TINY_DIST) new_origin(2) = 0.0_rk
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
