module mod_mach_cone_geometry
  !< Summary: This module contains the procedures to calculate the geometry
  !< of the Mach cone for the solution of the Euler equations using the 
  !< characteristic lines. Refer to the Appendix of DOI: 10.1016/j.jcp.2009.04.001. The 
  !< entire purpose of this module is to calculate the quantities \(\theta_{i,b}\)
  !< and \(\theta_{i,e}\)
  
  use iso_fortran_env, only : int32, real64
  use math_constants, only : pi
  implicit none
  
  ! private
  ! public :: 
  

  real(real64), dimension(2) :: l


contains
  

  pure subroutine calculate_l_parameter(u_tilde,v_tilde,tau,alpha_k,speed_of_sound,l_parameter)
    !< Calculate the \( \ell_{k} \)] parameter

    real(real64), intent(in) :: u_tilde, v_tilde !<
    real(real64), intent(in) :: tau !< Time increment, e.g. (x,y,t+tau)
    real(real64), intent(in) :: alpha_k !< angle of e_hat with respect to the x-axis 
    real(real64), intent(in) :: speed_of_sound

    real(real64), dimension(2), intent(out) :: l_parameter

    real(real64) :: cos_alpha, sin_alpha

    ! Precompute for reuse
    cos_alpha = cos(alpha_k)
    sin_alpha = sin(alpha_k)

    associate(a_tilde => speed_of_sound, l=> l_parameter)
    
      l(1) = -tau * (u_tilde*cos_alpha + v_tilde*sin_alpha) - &
                tau * sqrt(a_tilde**2 - ((u_tilde*sin_alpha - v_tilde*cos_alpha)**2))

      l(2) = -tau * (u_tilde*cos_alpha + v_tilde*sin_alpha) + &
                tau * sqrt(a_tilde**2 - ((u_tilde*sin_alpha - v_tilde*cos_alpha)**2))
    
    end associate


  end subroutine

  pure function calculate_theta_kl(alpha_k, speed_of_sound, u, l_k) result(theta_kl)

    real(real64), intent(in) :: alpha_k
    real(real64), intent(in) :: speed_of_sound
    real(real64), intent(in) :: u
    real(real64), intent(in) :: l_k
    real(real64) :: theta_kl

    theta_kl = pi + sign(1.0_real64,sin(alpha_k)) * &
               (acos((u + l_k*cos(alpha_k)) / speed_of_sound) -  pi)

  end function

  pure function determine_n_intersections(b_k, speed_of_sound, l_ki) result(n_intersections)
    !< Summary: Find the number of intersections based on the input parameters
    !< Refer to the appendix for this particular logic. It makes a bit more sense in context

    real(real64), intent(in) :: speed_of_sound
    real(real64), intent(in) :: b_k
    real(real64), dimension(2), intent(in) :: l_ki

    integer(int32) :: n_intersections

    if ((b_k > speed_of_sound) .or. &
        ((l_ki(1) <= l_ki(2)) .and. (l_ki(2) < 0.0_real64))) then
      n_intersections = 0
    else if ((l_ki(2) >= 0.0_real64) .and. (l_ki(1) < 0.0_real64)) then
      n_intersections = 1
    else
      n_intersections = 2
    end if

  end function

end module mod_mach_cone_geometry