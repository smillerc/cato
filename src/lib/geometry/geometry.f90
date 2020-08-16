! MIT License
! Copyright (c) 2019 Sam Miller
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

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

  ! subroutine find_midpoint_arc_segments(edge_vector, origin_in_cell, &
  !                                       circle_xy, circle_radius, arc_segments, n_arcs_in_cell)
  !   !< Since the midpoint mach cone is always on a single straight edge, we can avoid
  !   !< some unnecessary duplicate work.

  !   real(rk), dimension(2, 2), intent(in) :: edge_vector
  !   !< ((x,y), (point 1, point 2)); position of the line points

  !   logical, dimension(2), intent(in) :: origin_in_cell
  !   !< (cell 1, cell 2); is the mach cone origin in this cell?

  !   real(rk), dimension(2), intent(in) :: circle_xy
  !   !< (x,y); position of the mach cone origin

  !   real(rk), intent(in) :: circle_radius
  !   !< radius of the mach cone

  !   real(rk), dimension(2, 2), intent(in) :: arc_segments
  !   !< ((theta begin, theta end), (cell 1, cell 2))

  !   integer(ik), dimension(2), intent(out) :: n_arcs_in_cell
  !   !< (cell 1, cell 2); number of arcs contained in each cell

  ! end subroutine find_midpoint_arc_segments

  pure subroutine find_arc_segments(origin, vector_1, vector_2, origin_in_cell, &
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
    integer(ik) :: i
    logical, dimension(2) :: valid_intersections
    logical, dimension(2, 2) :: total_valid_intersections !< ((intersection 1, intersection 2), (line 1, line 2)); valid intersections for each line

    real(rk), dimension(2, 2) :: intersection_xy
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
    call intersect(line_xy=line, circle_xy=circle_xy, circle_radius=circle_radius, &
                   intersection_xy=intersection_xy, valid_intersections=valid_intersections, &
                   angles=intersection_angles_per_line)
    do i = 1, 2
      if(valid_intersections(i)) then
        n_intersections(1) = n_intersections(1) + 1
      end if
    end do
    intersection_angles(1, :) = intersection_angles_per_line
    total_valid_intersections(:, 1) = valid_intersections

    ! Line 2
    line(:, 2) = vector_2
    ! For a given line & circle intersection, find the angle with respect to the x-axis for each intersection point
    call intersect(line_xy=line, circle_xy=circle_xy, circle_radius=circle_radius, &
                   intersection_xy=intersection_xy, valid_intersections=valid_intersections, &
                   angles=intersection_angles_per_line)
    do i = 1, 2
      if(valid_intersections(i)) then
        n_intersections(2) = n_intersections(2) + 1
      end if
    end do
    intersection_angles(2, :) = intersection_angles_per_line
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

    associate(theta_ib => theta_start_end(1, :), theta_ie => theta_start_end(2, :), &
              n1 => n_intersections(1), n2 => n_intersections(2))

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
      associate(x2 => origins(:, 2), x1 => origins(:, 1), &
                r2 => radii(2), r1 => radii(1))
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

    associate(x1 => origins(1, larger_circle), y1 => origins(2, larger_circle), &
              x2 => origins(1, smaller_circle), y2 => origins(2, smaller_circle), &
              r2 => radii(smaller_circle), r1 => radii(larger_circle))

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
