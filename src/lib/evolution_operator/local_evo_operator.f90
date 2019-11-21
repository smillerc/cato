module mod_local_evo_operator
  use iso_fortran_env, only: ik => int32, rk => real64
  use math_constants, only: pi
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_input, only: input_t
  use mod_mach_cone_geometry, only: mach_cone_geometry_t
  use mod_grid, only: grid_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_reconstruction_factory, only: reconstruction_factory_t

  implicit none

  private
  public :: local_evo_operator_t

  type(reconstruction_factory_t) :: recon_factory

  type, extends(abstract_evo_operator_t) :: local_evo_operator_t
    ! private
    ! ! character(:), allocatable :: reconstruction_type

    ! real(rk), dimension(4) :: reference_state !> Refrence state of [rho,u,v,a]
    ! real(rk), dimension(4) :: reconstructed_primative_state !> [rho,u,v,p] evaluated at corner or midpoint node

  contains
    procedure, public :: evolve
    procedure, private :: get_pressure
    procedure, private :: get_density
    procedure, private :: get_x_velocity
    procedure, private :: get_y_velocity
    ! procedure, private :: find_p_prime_location
    final :: finalize
  end type local_evo_operator_t

  interface local_evo_operator
    module procedure :: constructor
  end interface
contains

  type(local_evo_operator_t) function constructor(input, grid, conserved_vars) result(operator)
    !< Constructor for the FVLEG operator
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid
    real(rk), dimension(:, :, :), target :: conserved_vars

    ! ! allocate(operator)
    operator%name = "FVLEG"
    call operator%set_tau(input%tau)
    operator%grid => grid
    operator%conserved_vars => conserved_vars

    ! p_prime_reconstruction_operator = recon_factory

  end function constructor

  pure function evolve(self, reconstructed_state, i, j, edge, loc) result(U_bar)
    !< Evolve the interface (edge) quantities based on their local reconstruction and neighboring cells

    real(rk), dimension(4) :: U_bar !< Conserved variable values specified location

    class(local_evo_operator_t), intent(in) :: self
    real(rk), dimension(:, :, :, :, :), intent(in) :: reconstructed_state
    !< ((rho, u ,v, p), point, node/midpoint, i, j);
    !< The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints

    integer(ik), intent(in) :: i, j !< Cell indices
    integer(ik), intent(in) :: edge !< Cell edge index
    character(len=3), intent(in) :: loc !< where is this evolution taking place? Needs to be one of the following ['k,1','k,c','k,2']

    real(rk) :: rho, u, v, p
    real(rk) :: x, y
    integer(ik) :: alloc_stat

    ! select case(loc)
    ! case("k,1") ! 1st corner
    !   self%mach_cone = mach_cone_geometry_t(edge_points=self%grid%cell_edge_vectors(:,1,edge,i,j))
    ! case("k,c") ! midpoint
    !   self%mach_cone = mach_cone_geometry_t(edge_points=self%grid%cell_edge_vectors(:,2,edge,i,j))
    ! case("k,2") ! 2nd corner
    !   self%mach_cone = mach_cone_geometry_t(edge_points=self%grid%cell_edge_vectors(:,3,edge,i,j))
    ! case default
    !   error stop "Invalid location for the local evolution operator E0 "// &
    !     "(Needs to be one of the following ['k,1','k,c','k,2'])"
    ! end select

    ! p = self%get_pressure()
    ! rho = self%get_density(reconstruction_operator)
    ! u = self%get_x_velocity()
    ! v = self%get_y_velocity

    ! U_bar = [rho, u, v, p]

    ! if (allocated(self%p_prime_reconstruction_operator)) then
    !   deallocate(self%p_prime_reconstruction_operator, stat=alloc_stat)
    !   if (alloc_stat /= 0) error stop 'Unable to deallocate local_evo_operator_t%p_prime_reconstruction_operator'
    ! end if
  end function evolve

  subroutine finalize(self)
    !< Cleanup the operator type
    type(local_evo_operator_t), intent(inout) :: self
    ! deallocate(self%mach_cone)
    nullify(self%grid)
    nullify(self%conserved_vars)
    nullify(self%reconstruction_operator)
  end subroutine finalize

  pure function get_pressure(self) result(pressure)
    !< Implementation of the p(P) within the local evolution operator (Eq 45 in the text)

    class(local_evo_operator_t), intent(in) :: self
    real(rk) :: pressure
    integer(ik) :: i, arc

    pressure = 0.0_rk

    do arc = 1, self%mach_cone%n_arcs
      do i = 1, self%mach_cone%n_intersections(arc)

        associate(theta_ie=>self%mach_cone%theta_ie(arc, i), theta_ib=>self%mach_cone%theta_ib(arc, i), &
                  p_i=>self%mach_cone%arc_conserved_vars(4, arc, i), &
                  u_i=>self%mach_cone%arc_conserved_vars(2, arc, i), &
                  v_i=>self%mach_cone%arc_conserved_vars(1, arc, i), &
                  rho_tilde=>self%mach_cone%reference_state(1), &
                  a_tilde=>self%mach_cone%reference_state(4))

          pressure = pressure + p_i * (theta_ie - theta_ib) - &
                     rho_tilde * a_tilde * u_i * (sin(theta_ie) - sin(theta_ib)) + &
                     rho_tilde * a_tilde * v_i * (cos(theta_ie) - cos(theta_ib))
        end associate
      end do
    end do

    pressure = pressure / (2.0_rk * pi)

  end function get_pressure

  pure function get_density(self) result(density)
    !< Implementation of the rho(P) within the local evolution operator (Eq. 42 in the text)
    class(local_evo_operator_t), intent(in) :: self

    real(rk) :: density  !< rho(P)
    integer(ik) :: i, arc
    integer(ik), dimension(2) :: p_prime_ij
    real(rk), dimension(4) :: p_prime_u_bar !< [rho, u, v, p] at P'
    logical :: p_prime_found

    density = 0.0_rk

    do arc = 1, self%mach_cone%n_arcs
      do i = 1, self%mach_cone%n_intersections(arc)

        associate(theta_ie=>self%mach_cone%theta_ie(arc, i), theta_ib=>self%mach_cone%theta_ib(arc, i), &
                  p_i=>self%mach_cone%arc_conserved_vars(4, arc, i), &
                  u_i=>self%mach_cone%arc_conserved_vars(2, arc, i), &
                  v_i=>self%mach_cone%arc_conserved_vars(1, arc, i), &
                  rho_tilde=>self%mach_cone%reference_state(1), &
                  a_tilde=>self%mach_cone%reference_state(4))

          density = density + (p_i / a_tilde**2) * (theta_ie - theta_ib) - &
                    (rho_tilde / a_tilde) * u_i * (sin(theta_ie) - sin(theta_ib)) + &
                    (rho_tilde / a_tilde) * v_i * (cos(theta_ie) - cos(theta_ib))
        end associate
      end do
    end do

    density = density / (2.0_rk * pi)

    ! Reconstruct at P' to get rho(P') and p(P')
    p_prime_u_bar = self%reconstruction_operator%reconstruct_point(xy=self%mach_cone%p_prime_xy, &
                                                                   cell_ij=self%mach_cone%p_prime_ij)

    associate(rho_p_prime=>p_prime_u_bar(1), pressure_p_prime=>p_prime_u_bar(4), &
              a_tilde=>self%mach_cone%reference_state(4))
      density = density + rho_p_prime - (pressure_p_prime / a_tilde**2)
    end associate

  end function get_density

  pure function get_x_velocity(self) result(u)
    !< Implementation of the u(P) within the local evolution operator (Eq 43 in the text)

    class(local_evo_operator_t), intent(in) :: self
    real(rk) :: u !< x velocity at P
    integer(ik) :: i
    integer(ik) :: arc  !< arc index within the mach cone

    u = 0.0_rk

    do arc = 1, self%mach_cone%n_arcs
      do i = 1, self%mach_cone%n_intersections(arc)

        associate(theta_ie=>self%mach_cone%theta_ie(arc, i), theta_ib=>self%mach_cone%theta_ib(arc, i), &
                  p_i=>self%mach_cone%arc_conserved_vars(4, arc, i), &
                  u_i=>self%mach_cone%arc_conserved_vars(2, arc, i), &
                  v_i=>self%mach_cone%arc_conserved_vars(1, arc, i), &
                  rho_tilde=>self%mach_cone%reference_state(1), &
                  a_tilde=>self%mach_cone%reference_state(4))

          u = u - (p_i / (rho_tilde * a_tilde)) * (sin(theta_ie) - sin(theta_ib)) + &
              u_i * (0.5_rk * (theta_ie - theta_ib) + 0.25_rk * (sin(2 * theta_ie) - sin(2 * theta_ib))) - &
              v_i * 0.25_rk * (cos(2 * theta_ie) - cos(2 * theta_ib))

        end associate
      end do
    end do

    u = u / pi

  end function get_x_velocity

  pure function get_y_velocity(self) result(v)
    !< Implementation of the v(P) within the local evolution operator (Eq 44 in the text)

    class(local_evo_operator_t), intent(in) :: self
    real(rk) :: v !< y velocity at P
    integer(ik) :: i
    integer(ik) :: arc  !< arc index within the mach cone

    v = 0.0_rk

    do arc = 1, self%mach_cone%n_arcs
      do i = 1, self%mach_cone%n_intersections(arc)

        associate(theta_ie=>self%mach_cone%theta_ie(arc, i), theta_ib=>self%mach_cone%theta_ib(arc, i), &
                  p_i=>self%mach_cone%arc_conserved_vars(4, arc, i), &
                  u_i=>self%mach_cone%arc_conserved_vars(2, arc, i), &
                  v_i=>self%mach_cone%arc_conserved_vars(1, arc, i), &
                  rho_tilde=>self%mach_cone%reference_state(1), &
                  a_tilde=>self%mach_cone%reference_state(4))

          v = v + (p_i / (rho_tilde * a_tilde)) * (cos(theta_ie) - cos(theta_ib)) + &
              u_i * 0.25_rk * (cos(2 * theta_ie) - cos(2 * theta_ib)) + &
              v_i * (0.5_rk * (theta_ie - theta_ib) + 0.25_rk * (sin(2 * theta_ie) - sin(2 * theta_ib)))
        end associate
      end do
    end do

    v = v / pi
  end function get_y_velocity

end module mod_local_evo_operator
