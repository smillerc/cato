module mod_mach_cone_utilties

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use math_constants, only: pi, rad2deg
  use mod_globals, only: TINY_DIST
  use mod_vector, only: vector_t, operator(.cross.)

  implicit none

  logical :: enable_tau_scaling = .false. !< automatically scale the tau factor to ensure nice numbers near 1 for floating point accuracy
  logical :: enable_edge_vector_scaling = .true. !< scale the edge vector lengths so that they are near 1 (based off of the average length)
  real(rk), parameter :: tau_scaling = 0.01 !< scaling factor for tau that ensures P' is within the neighboring cells

contains

  pure subroutine get_cone_extents(tau, xy, vel, sound_speed, origin, radius)
    !< Given a velocity and sound speed, determine the extents of the mach cone

    real(rk), intent(in) :: tau                 !< time increment
    real(rk), dimension(2), intent(in) :: xy    !< P xy location
    real(rk), dimension(2), intent(in) :: vel   !< xy velocity
    real(rk), intent(in) :: sound_speed         !< sound speed

    real(rk), dimension(2), intent(out) :: origin      !< origin of the cone/circle
    real(rk), intent(out) :: radius      !< radius of the cone/circle

    associate(x=>xy(1), y=>xy(2), &
              u=>vel(1), v=>vel(2))
      origin = [x - tau * u, y - tau * v]
      radius = sound_speed * tau

      ! For very small distances, just make it coincide with P
      if(abs(origin(1) - x) < TINY_DIST) origin(1) = x
      if(abs(origin(2) - y) < TINY_DIST) origin(2) = y
    end associate
  end subroutine get_cone_extents

end module mod_mach_cone_utilties
