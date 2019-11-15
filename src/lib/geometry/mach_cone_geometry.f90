module mod_mach_cone_geometry
  !< Summary: This module contains the procedures to calculate the geometry
  !< of the Mach cone for the solution of the Euler equations using the
  !< characteristic lines. Refer to the Appendix of DOI: 10.1016/j.jcp.2009.04.001. The
  !< entire purpose of this module is to calculate the quantities \(\theta_{i,b}\)
  !< and \(\theta_{i,e}\)

  use iso_fortran_env, only: ik => int32, rk => real64
  use math_constants, only: pi
  use mod_vector, only: vector_t, operator(.unitnorm.), operator(.dot.), operator(.cross.), angle_between

  implicit none

  private
  public :: mach_cone_geometry_t

  type :: mach_cone_geometry_t
    !< Type to encapsulate and calculate the geometry of the mach cone for the evolution operators
    private
    real(rk), dimension(2), public :: theta_ie = 0.0_rk
    real(rk), dimension(2), public :: theta_ib = 0.0_rk
    real(rk), dimension(2), public :: p_prime_xy = 0.0_rk
    !< x,y location of P' or the apex of the Mach cone (global, not relative to P0)
    integer(ik), dimension(2), public :: n_intersections = 0
    logical :: p_prime_in_control_volume = .false.  !< is P' inside the control volume?

    real(rk) :: tau
    real(rk), dimension(2, 2) :: theta
    ! real(rk), dimension(2) :: alpha_k
    real(rk), dimension(2) :: sin_alpha_k
    real(rk), dimension(2) :: cos_alpha_k
    real(rk) :: reference_speed_of_sound !< a_tilde
    real(rk), dimension(2) :: b !< b parameter (see appendix in text)
    real(rk), dimension(2, 2) :: l
    real(rk), dimension(2) :: reference_velocity !< [u_tilde, v_tilde]
    type(vector_t), dimension(2) :: edge_vectors !< edge vectors relative to P0
    type(vector_t) :: p_prime_pos_vector !< vector defining P' location relative to P0
  contains
    procedure, public :: initialize
    procedure, private :: calculate_l_parameter
    procedure, private :: calculate_theta_kl
    procedure, private :: determine_n_intersections
    procedure, private :: calculate_thetas
    procedure, private :: determine_if_p_prime_is_in_cell
  end type

contains

  pure function initialize(self, e1, e2, reference_state_uva, tau)
    !< Implementation of the mach cone constructor
    class(mach_cone_geometry_t), intent(inout) :: self
    real(rk), dimension(2, 2), intent(in) :: e1 !< 1st edge vector [[x0,y0],[x1,y1]]
    real(rk), dimension(2, 2), intent(in) :: e2 !< 2nd edge vector [[x0,y0],[x1,y1]]
    real(rk), dimension(3), intent(in) :: reference_state_uva !< [u, v, speed_of_sound]
    real(rk), intent(in) :: tau !< time increment

    ! Although these are specified as a pair of x,y points, they will get shifted to the origin
    self%edge_vectors(1) = vector_t(x=e1[:, 1], y=e1[:, 2])
    self%edge_vectors(2) = vector_t(x=e2[:, 1], y=e2[:, 2])

    self%reference_velocity = reference_state_uva(1:2)
    self%reference_speed_of_sound = reference_state_uva(3)

    ! Find the (x,y) location of the center of the Mach circle (P')
    associate(x=>e1(1, 1), y=>e1(1, 2), &
              u_tilde=>self%reference_velocity(1), v_tilde=>self%reference_velocity(2))

      ! This defines P' (x,y) globally, not with respect to P0
      self%p_prime_xy = [x - u_tilde * tau, y - v_tilde * tau]
    end associate

    ! Find the angle with respect to the x-axis for each edge vector
    do k = 1, 2
      associate(x=>self%edge_vectors(k)%x, y=>self%edge_vectors(k)%y, len=>self%edge_vectors(k))

        self%sin_alpha_k(k) = y / len
        self%cos_alpha_k(k) = x / len
      end associate
    end do

    call self%determine_if_p_prime_is_in_cell()
    call self%calculate_l_parameter()
    call self%calculate_theta_kl()
    call self%determine_n_intersections()
    call self%calculate_arc_thetas()
  end function initialize

  subroutine calculate_arc_thetas(self)
    !< Implementation of the logic used to find the arc angles theta ib and ie (begin and end)
    class(mach_cone_geometry_t), intent(inout) :: self

    associate(b=>self%b, n1=>self%n_intersections(1), n2=>self%n_intersections(2), &
              theta=>self%theta, &  ! this is the theta k,l term in the appendix
              theta_ib=>self%theta_ib, theta_ie=>self%theta_ie)

      if(n1 == 0 .and. n2 == 0 .and. self%p_prime_in_cell) then
        theta_ib = 0.0_rk
        theta_ie = 2 * pi

      else if(n1 == 1 .and. n2 == 1) then
        theta_ib = theta(1, 2)

        if(theta(2, 2) >= theta(1, 2)) then
          theta_ie = theta(2, 2)
        else
          theta_ie = theta(2, 2) + 2 * pi
        endif

      else if(n1 == 2 .and. n2 == 0) then
        theta_ib = theta(1, 2)

        if(theta(1, 1) >= theta(1, 2)) then
          theta_ie = theta(1, 1)
        else
          theta_ie = theta(1, 1) + 2 * pi
        endif

      else if(n1 == 0 .and. n2 == 2) then
        theta_ib = theta(2, 1)

        if(theta(2, 2) >= theta(2, 1)) then
          theta_ie = theta(2, 2)
        else
          theta_ie = theta(2, 2) + 2 * pi
        endif

      else if(n1 == 2 .and. n2 == 2) then
        ! There are 2 arcs for each edge vector

        theta_ib(1) = theta(1, 2)
        theta_ib(2) = theta(2, 1)

        ! 1st arc
        if(theta(2, 2) >= theta(1, 2)) then
          theta_ie(1) = theta(2, 2)
        else
          theta_ie(1) = theta(2, 2) + 2 * pi
        endif

        ! 2nd arc
        if(theta(1, 1) >= theta(2, 1)) then
          theta_ie(2) = theta(1, 1)
        else
          theta_ie(2) = theta(1, 1) + 2 * pi
        endif

      end if
    end associate

  end subroutine calculate_arc_thetas

  subroutine calculate_l_parameter(self)
    !< Calculate the \( \ell_{k} \)] parameter
    class(mach_cone_geometry_t), intent(inout) :: self
    real(rk) :: cos_alpha, sin_alpha

    ! Determine if there are real-value solutions to l(k)
    do k = 1, 2
      associate(sin_alpha_k=>self%sin_alpha_k(k), cos_alpha_k=>self%cos_alpha_k(k), &
                u_tilde=>self%reference_velocity(1), v_tilde=>self%reference_velocity(2))
        b(k) = abs(u_tilde * sin_alpha(k) + v_tilde * cos_alpha(k))
      end associate
    end do

    self%l = 0.0_rk

    do k = 1, 2
      associate(a_tilde=>self%reference_speed_of_sound, l=>self%l_k, &
                tau=>self%tau, &
                cos_alpha=>self%cos_alpha_k, sin_alpha=>self%sin_alpha_k, &
                u_tilde=>self%reference_velocity(1), v_tilde=>self%reference_velocity(2))

        if(b(k) <= a_tilde) then
          l(k, 1) = -tau * (u_tilde * cos_alpha(k) + v_tilde * sin_alpha(k)) - &
                    tau * sqrt(a_tilde**2 - ((u_tilde * sin_alpha(k) - v_tilde * cos_alpha(k))**2))

          l(k, 2) = -tau * (u_tilde * cos_alpha(k) + v_tilde * sin_alpha(k)) + &
                    tau * sqrt(a_tilde**2 - ((u_tilde * sin_alpha(k) - v_tilde * cos_alpha(k))**2))
        end if

      end associate
    end do
  end subroutine calculate_l_parameter

  subroutine calculate_theta_kl(self)
    !< Implementation of the \(\theta_{k,\ell}\) term in the Appendix of the text
    class(mach_cone_geometry_t), intent(inout) :: self
    real(rk) :: sine_alpha_kl
    integer(ik) :: k, l

    do k = 1, 2
      do l = 1, 2
        associate(a_tilde=>self%reference_speed_of_sound, l=>self%l, &
                  alpha_k=>self%alpha_k, theta=>self%theta, &
                  u_tilde=>self%reference_velocity(1), v_tilde=>self%reference_velocity(2), &
                  sin_alpha=>self%sin_alpha_k, cos_alpha=>self%cos_alpha_k)

          sine_alpha_kl = v_tilde + l(k, l) * sin_alpha(k) / a_tilde

          theta(k, l) = pi + sign(1.0_rk, sine_alpha_kl) * &
                        (acos((u_tilde + l(k, l) * cos(alpha_k)) / a_tilde) - pi)
        end associate
      end do
    end do
  end subroutine calculate_theta_kl

  subroutine determine_n_intersections(self)
    !< Summary: Find the number of intersections based on the input parameters
    !< Refer to the appendix for this particular logic
    class(mach_cone_geometry_t), intent(inout) :: self
    integer(ik) :: k

    do k = 1, 2
      associate(b=>self%b, a_tilde=>self%reference_speed_of_sound, l=>self%l, &
                n=>self%n_intersections)

        if((b(k) > a_tilde) .or. &
           ((l(k, 1) <= l(k, 2)) .and. (l(k, 2) < 0.0_rk))) then
          n(k) = 0
        else if((l(k, 2) >= 0.0_rk) .and. (l(k, 1) < 0.0_rk)) then
          n(k) = 1
        else
          n(k) = 2
        end if
      end associate
    end do

  end subroutine determine_n_intersections

  subroutine determine_if_p_prime_is_in_cell(self)
    !< Implementation of whether the P' point is inside the current cell/control volume. This
    !< uses the cross product of 2 vectors in 2d, which gives a scalar
    class(mach_cone_geometry_t), intent(inout) :: self
    type(vector_t) :: ref_vel

    ref_vel = vector_t(x=self%p_prime_xy(1), y=self%p_prime_xy(2))

    ! In the text, these inequalities are swapped, but I define the P' position vector
    ! as starting at P0 and going to P'
    if((ref_vel.cross.self%edge_vectors(1) <= 0.0_rk) .and. &
       (self%edge_vectors(2) .cross.ref_vel <= 0.0_rk)) then
      self%p_prime_in_control_volume = .true.
    else
      self%p_prime_in_control_volume = .false.
    end if
  end subroutine determine_if_p_prime_is_in_cell

end module mod_mach_cone_geometry
