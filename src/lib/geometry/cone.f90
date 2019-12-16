module mod_cone
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use math_constants, only: pi
  use mod_floating_point_utils, only: near_zero
  use mod_eos, only: eos

  implicit none

contains

  subroutine vector_circle_intersect(line_xy, circle_xy, circle_radius, intersection_xy, valid_intersection)
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

end module mod_cone
