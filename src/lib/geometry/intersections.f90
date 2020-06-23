module mod_intersections
  !< Summary: Provide procedures to calculate the intersection of a line segment and a circle
  !< Date: 05/20/2020
  !< Author: Sam Miller
  !< References:
  !<      [1] https://cp-algorithms.com/geometry/circle-line-intersection.html
  !< Notes: Reference [1] had some errors in it, but otherwise this is better than using
  !<        the quadratic formula to find the intersects. Using the quadratic formula is by far
  !<        the most common algorithim online.

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  implicit none

  private
  public :: intersect, is_on

  real(rk), parameter :: EPS = epsilon(1.0_rk)
  real(rk), parameter :: EPS_2X = 1e-15_rk

contains

  pure subroutine intersect(line_xy, circle_xy, circle_radius, intersection_xy, valid_intersections, angles)
    !< Find the interesection(s) of a line and a circle. This subroutine shifts and scales everything so that the
    !< origin of the circle is at (0,0) and the radius is 1. It provides the interesection (x,y) points, a boolean
    !< if the intersections are valid or not (it will find an intersection based on an infinite line), and the angle
    !< between the intersection point and the x-axis. This was modified from [1] (which had some errors in it)

    real(rk), dimension(2, 2), intent(in) :: line_xy !< ((x,y) (point_1, point_2)); Line location
    real(rk), dimension(2), intent(in) :: circle_xy !< (x,y); Circle origin
    real(rk), intent(in) :: circle_radius !< Radius of the circle
    real(rk), dimension(2, 2), intent(out) :: intersection_xy !< ((x,y), (point_1, point_2)); Intersection point
    logical, dimension(2), intent(out) :: valid_intersections !< (point_1, point_2); Is there an intersection or not?
    real(rk), dimension(2), intent(out) :: angles !< (point_1, point_2); Angle w/respect to the x-axis

    ! Locals
    integer(ik) :: n_intersections
    !< number of intersections (not all may be valid yet, since it assumes the line is infinite at first)

    real(rk), parameter :: r_sq = 1.0_rk !< r**2. This makes the code easier to read
    real(rk) :: x1, y1, x2, y2
    real(rk) :: ax, ay !< 1st intersection point
    real(rk) :: bx, by !< 2nd intersection point
    real(rk) :: x0, y0, d, m
    real(rk) :: A, B, C !< Ax + By = C
    real(rk) :: threshold !< filter to check for roundoff

    ! Defaults
    ax = 0.0_rk
    bx = 0.0_rk
    ay = 0.0_rk
    by = 0.0_rk
    angles = 0.0_rk
    n_intersections = 0
    intersection_xy = 0.0_rk
    valid_intersections = .false.

    x1 = line_xy(1, 1) - circle_xy(1)
    threshold = abs(line_xy(1, 1) + circle_xy(1)) * EPS
    if(abs(x1) < threshold) then
      x1 = 0.0_rk
    else
      x1 = x1 / circle_radius
    end if

    y1 = line_xy(2, 1) - circle_xy(2)
    threshold = abs(line_xy(2, 1) + circle_xy(2)) * EPS
    if(abs(y1) < threshold) then
      y1 = 0.0_rk
    else
      y1 = y1 / circle_radius
    end if

    x2 = line_xy(1, 2) - circle_xy(1)
    threshold = abs(line_xy(1, 2) + circle_xy(1)) * EPS
    if(abs(x2) < threshold) then
      x2 = 0.0_rk
    else
      x2 = x2 / circle_radius
    end if

    y2 = line_xy(2, 2) - circle_xy(2)
    threshold = abs(line_xy(2, 2) + circle_xy(2)) * EPS
    if(abs(y2) < threshold) then
      y2 = 0.0_rk
    else
      y2 = y2 / circle_radius
    end if

    A = y2 - y1
    if(abs(A) < epsilon(1.0_rk)) A = 0.0_rk

    B = x1 - x2
    if(abs(B) < epsilon(1.0_rk)) B = 0.0_rk

    C = x1 * y2 - x2 * y1
    if(abs(C) < epsilon(1.0_rk)) C = 0.0_rk

    if(C**2 > (A**2 + B**2)) then
      n_intersections = 0
    else if(abs(C**2 - (A**2 + B**2)) < epsilon(1.0_rk)) then
      n_intersections = 1
      x0 = (A * C) / (A**2 + B**2)
      y0 = (B * C) / (A**2 + B**2)
      ax = x0
      ay = y0
    else
      n_intersections = 2
      x0 = (A * C) / (A**2 + B**2)
      y0 = (B * C) / (A**2 + B**2)
      d = sqrt(r_sq - ((C**2) / (A**2 + B**2)))
      m = sqrt(d**2 / (A**2 + B**2))
      ax = x0 + B * m
      bx = x0 - B * m
      ay = y0 - A * m
      by = y0 + A * m
    end if

    ! Get rid of very small numbers and signed 0's, i.e. -0
    if(abs(ax) < tiny(1.0_rk)) ax = 0.0_rk
    if(abs(ay) < tiny(1.0_rk)) ay = 0.0_rk
    if(abs(bx) < tiny(1.0_rk)) bx = 0.0_rk
    if(abs(by) < tiny(1.0_rk)) by = 0.0_rk

    intersection_xy(:, 1) = ([ax, ay] * circle_radius) + circle_xy
    intersection_xy(:, 2) = ([bx, by] * circle_radius) + circle_xy

    if(n_intersections > 0) then
      valid_intersections(1) = is_on(a=[x1, y1], b=[x2, y2], c=[ax, ay])
      valid_intersections(2) = is_on(a=[x1, y1], b=[x2, y2], c=[bx, by])

      if(valid_intersections(1)) angles(1) = atan2(y=ay, x=ax)
      if(valid_intersections(2)) angles(2) = atan2(y=by, x=bx)
    end if

  end subroutine intersect

  pure logical function between_zero_and_one(val)
    !< Stupidly simple function to check if a value is between 0 and 1 with
    !< machine epsilon checks included

    real(rk), intent(in) :: val

    if(val < 0.0_rk .or. (val - 1.0_rk) > EPS) then
      between_zero_and_one = .false.
    else
      between_zero_and_one = .true.
    end if
  end function between_zero_and_one

  pure logical function is_on(a, b, c)
    ! Check to see if the point c on the line segment ab

    real(rk), dimension(2), intent(in) :: a !< (x,y);
    real(rk), dimension(2), intent(in) :: b !< (x,y);
    real(rk), dimension(2), intent(in) :: c !< (x,y);

    real(rk) :: ab_dx !< delta x for the line segment AB
    real(rk) :: ab_dy !< delta y for the line segment AB
    real(rk) :: cb_dx !< delta x for the line segment CB
    real(rk) :: cb_dy !< delta y for the line segment CB

    real(rk) :: m1, m2 !< parameterization factor
    logical :: ab_is_horizontal
    logical :: ab_is_vertical

    ab_is_horizontal = .false.
    ab_is_vertical = .false.

    ab_dx = a(1) - b(1)
    ab_dy = a(2) - b(2)
    cb_dy = c(2) - b(2)
    cb_dx = c(1) - b(1)

    if(abs(ab_dy) < EPS_2X) ab_dy = 0.0_rk
    if(abs(cb_dx) < EPS_2X) cb_dx = 0.0_rk

    if(abs(ab_dy) < EPS_2X) then
      ab_is_horizontal = .true.
      ab_dy = 0.0_rk
    end if

    if(abs(ab_dx) < EPS_2X) then
      ab_is_vertical = .true.
      ab_dx = 0.0_rk
    end if

    if(abs(ab_dx) < EPS_2X .and. abs(ab_dy) < EPS_2X) then
      ! a and b are the same point,
      ! so check that c is the same as a and b
      is_on = (abs(cb_dx) < EPS_2X) .and. (abs(cb_dy) < EPS_2X)
    else
      if(ab_is_vertical) then
        m2 = cb_dy / ab_dy
        is_on = (abs(cb_dx) < EPS_2X) .and. between_zero_and_one(m2)
      else if(ab_is_horizontal) then
        m1 = cb_dx / ab_dx
        is_on = (abs(cb_dy) < EPS_2X) .and. between_zero_and_one(m1)
      else
        m1 = cb_dx / ab_dx
        if(.not. between_zero_and_one(m1)) then
          is_on = .false.
        else
          m2 = cb_dy / ab_dy
          is_on = abs(m2 - m1) < EPS_2X
        end if
      end if
    end if
  end function is_on

end module mod_intersections
