module mod_cone
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use math_constants, only: pi
  use mod_floating_point_utils, only: near_zero
  use mod_eos, only: eos

  implicit none

contains

  pure subroutine get_arc_segments(lines, origin_in_cell, circle_xy, circle_radius, arc_segments, n_arcs)
    !< Given 2 lines and a circle, find their intersections and starting/ending angles for each arc

    real(rk), dimension(2, 2, 2), intent(in) :: lines
    !< ((x,y), (tail,head), (line_1:line_2)); set of vectors that define the lines to intersect with the circle

    logical, intent(in) :: origin_in_cell  !< is the circle origin in the cell defined by the two lines?
    real(rk), dimension(2), intent(in) :: circle_xy       !< (x,y) origin of the circle
    real(rk), intent(in) :: circle_radius                 !< radius of the circle
    real(rk), dimension(2, 2), intent(out):: arc_segments !< ((start,end), (arc1, arc2))
    integer(ik), intent(out) :: n_arcs !< Number of arcs with a valid start and end angle

    integer(ik) :: n_intersections_per_line      !< # intersections for a single line
    integer(ik), dimension(2) :: n_intersections !< # intersections for all lines
    real(rk), dimension(2, 2) :: line

    real(rk), dimension(2) :: intersection_angles_per_line !< intersection angles
    real(rk), dimension(2, 2) :: intersection_angles          !< intersection angles for all lines
    integer(ik) :: i

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
