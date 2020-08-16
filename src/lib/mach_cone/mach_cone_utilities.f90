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

module mod_mach_cone_utilties

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use math_constants, only: pi, rad2deg
  use mod_vector, only: vector_t, operator(.cross.)

  implicit none

  real(rk), parameter :: TINY_DIST_RATIO = 1.0e-8_rk
contains

  subroutine get_cone_extents(tau, xy, vel, sound_speed, origin, radius)
    !< Given a velocity and sound speed, determine the extents of the mach cone

    real(rk), intent(in) :: tau                 !< time increment
    real(rk), dimension(2), intent(in) :: xy    !< P xy location
    real(rk), dimension(2), intent(in) :: vel   !< xy velocity
    real(rk), intent(in) :: sound_speed         !< sound speed

    real(rk), dimension(2), intent(out) :: origin      !< origin of the cone/circle
    real(rk), intent(out) :: radius      !< radius of the cone/circle
    real(rk) :: dx, dy

    associate(x => xy(1), y => xy(2), &
              u => vel(1), v => vel(2))
      origin = [x - tau * u, y - tau * v]
      radius = sound_speed * tau

      dx = abs(origin(1) - x)
      dy = abs(origin(2) - y)

      ! if(abs(dx / radius) < TINY_DIST_RATIO) origin(1) = x
      ! if(abs(dy / radius) < TINY_DIST_RATIO) origin(2) = y
      ! if(ieee_is_nan(radius)) then
      !   print *, 'Error: radius is NaN'
      !   print *, "origin: ", origin
      !   print *, "radius: ", radius
      !   print *, "dx    : ", dx
      !   print *, "dy    : ", dy
      !   print *, "cs    : ", sound_speed
      !   print *, "tau   : ", tau
      !   error stop "Radius is NaN!"
      ! end if
    end associate
  end subroutine get_cone_extents

  pure logical function determine_if_p_prime_is_in_cell(origin, vec_1_head, vec_2_head, p_prime_vector) result(in_cell)
    !< Implementation of whether the P' point is inside the current cell/control volume. This
    !< uses the cross product of 2 vectors in 2d, which gives a scalar
    type(vector_t), intent(in) :: p_prime_vector
    real(rk), dimension(2), intent(in) :: origin     !< (x,y) origin of both vectors
    real(rk), dimension(2), intent(in) :: vec_1_head !< (x,y) head of the 1st vector
    real(rk), dimension(2), intent(in) :: vec_2_head !< (x,y) head of the 2nd vector

    type(vector_t) :: edge_vector_1, edge_vector_2

    in_cell = .false.
    edge_vector_1 = vector_t(x=[origin(1), vec_1_head(1)], &
                             y=[origin(2), vec_1_head(2)])

    edge_vector_2 = vector_t(x=[origin(1), vec_2_head(1)], &
                             y=[origin(2), vec_2_head(2)])

    ! In the text is has >= instead of <= for some reason, but the following works
    ! like it's supposed to.
    if((p_prime_vector.cross.edge_vector_1) <= 0.0_rk .and. &
       (edge_vector_2.cross.p_prime_vector) <= 0.0_rk) then
      in_cell = .true.
    else
      in_cell = .false.
    end if
  end function determine_if_p_prime_is_in_cell

end module mod_mach_cone_utilties
