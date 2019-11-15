module mod_local_evo_operator
  use iso_fortran_env, only: ik => int32, rk => real64
  use math_constants, only: pi
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_input, only: input_t
  use mod_mach_cone_geometry, only: mach_cone_geometry_t
  use mod_grid, only: grid_t
  use mod_conserved_vars, only: conserved_vars_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  implicit none

  private
  public :: local_evo_operator_t

  type, extends(abstract_evo_operator_t) :: local_evo_operator_t
    private
    type(mach_cone_geometry_t) :: mach_cone !< type used to determine mach cone geometry
    class(grid_t), pointer :: grid !< grid type to track neighbor info
    class(conserved_vars_t), pointer :: conserved_vars !< aka U, the vector of conserved variables
    real(rk), dimension(4) :: reference_state !> Refrence state of [rho,u,v,a]
    real(rk), dimension(4) :: reconstructed_primative_state !> [rho,u,v,p] evaluated at corner or midpoint node

    class(abstract_reconstruction_t), allocatable :: p_prime_reconstruction_operator
    !< Reconstruction used at P', when P' is outside of the current cell; The current cell reconstruction operator is passed
    !< to the evolve function already. This is just used when needed
  contains
    procedure, public :: evolve
    procedure, private :: get_pressure
    procedure, private :: get_density
    procedure, private :: get_x_velocity
    procedure, private :: get_y_velocity
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
    class(conserved_vars_t), intent(in), target :: conserved_vars

    ! allocate(operator)
    operator%name = "FVLEG"
    call operator%set_tau(input%tau)
    operator%grid => grid
    operator%conserved_vars => conserved_vars

  end function constructor

  function evolve(self, i, j, reconstruction_operator, edge, loc) result(U_bar)
    !< Evolve the interface (edge) quantities based on their local reconstruction and neighboring cells

    class(local_evo_operator_t), intent(in) :: self
    integer(ik), intent(in) :: i, j !< Cell indices
    class(abstract_reconstruction_t), pointer, intent(in) :: reconstruction_operator

    ! While these are optional based on the interface, they are needed for this type of reconstruction (maybe I need to rework the design)
    real(rk), dimension(4) :: U_bar !< Conserved variable values at the interface (edge)
    integer(ik), intent(in), optional :: edge !< Cell edge index
    character(len=3), intent(in), optional :: loc !< where is this evolution taking place? Needs to be one of the following ['k,1','k,c','k,2']

    real(rk) :: rho, u, v, p
    real(rk) :: x, y
    real(rk), dimension(2) :: e1, e2 !< edge
    integer(ik), dimension(3) :: edge_1_map = [1, 2, 3]  !< indices for N1, M1, and N2 in the cell_node_xy array
    integer(ik), dimension(3) :: edge_2_map = [3, 4, 5]  !< indices for N2, M2, and N3 in the cell_node_xy array
    integer(ik), dimension(3) :: edge_3_map = [5, 6, 7]  !< indices for N3, M3, and N4 in the cell_node_xy array
    integer(ik), dimension(3) :: edge_4_map = [7, 8, 1]  !< indices for N4, M4, and N1 in the cell_node_xy array
    integer(ik), dimension(3) :: edge_map  !< node/midpoint indices to use for a given edge index

    ! This mapping is go from edge number (for a quad cell) to the index location in the cell_node_xy array
    ! Refer to the implementation of regular_2d_grid_t
    select case(edge)
    case(1)
      edge_map = edge_1_map
    case(2)
      edge_map = edge_2_map
    case(3)
      edge_map = edge_3_map
    case(4)
      edge_map = edge_4_map
    case default
      error stop "Invalid edge number in local_evo_operator_t%evolve"
    end select

    select case(loc)
    case("k,1") ! 1st corner

      e1 = [N1, N4]
      e2 = []

      ! x = self%grid%get_cell_node_xy(i,j,edge_map(1), 1)
      ! y = self%grid%get_cell_node_xy(i,j,edge_map(1), 2)
    case("k,c") ! midpoint
      ! x = self%grid%get_cell_node_xy(i,j,edge_map(2), 1)
      ! y = self%grid%get_cell_node_xy(i,j,edge_map(2), 2)
    case("k,2") ! 2nd corner

      associate(x_head_1=>, y_head_1=>)
        e1 = []
      end associate
      ! x = self%grid%get_cell_node_xy(i,j,edge_map(3), 1)
      ! y = self%grid%get_cell_node_xy(i,j,edge_map(3), 2)
    case default
      error stop "Invalid location for the local evolution operator E0 "// &
        "(Needs to be one of the following ['k,1','k,c','k,2'])"
    end select

    U_bar = [rho, u, v, p]
  end function

  ! pure function construct_edge_vector()
  ! end function

  subroutine finalize(self)
    !< Cleanup the operator type
    type(local_evo_operator_t), intent(inout) :: self
    deallocate(self%mach_cone)
    nullify(self%grid)
  end subroutine

  pure function find_p_prime_location(self, x, y) result(ij)
    ! //TODO: implement P' finder!
    !< This is to find out where P' is so that density and pressure can be reconstructed for
    !> use in the get_density, or rho(P), formula (Eq. 42 in the text)
    real(rk), intent(in) :: x, y !< Location of P'
    integer(ik), dimension(2) :: ij !< Coordinates of the cell where P' lives
  end function

  pure function get_pressure(self) result(pressure)
    !< Implementation of the p(P) within the local evolution operator (Eq 45 in the text)

    class(local_evo_operator_t), intent(in) :: self
    real(rk) :: pressure
    integer(ik) :: i

    pressure = 0.0_rk

    do i = 1, self%mach_cone%n_intersections

      associate(theta_ie=>self%mach_cone%theta_ie(i), theta_ib=>self%mach_cone%theta_ib(i), &
                p_i=>self%reconstructed_primative_state[4], &
                u_i=>self%reconstructed_primative_state[2], v_i=>self%reconstructed_primative_state[1], &
                rho_tilde=>self%reference_state[1], a_tilde=>self%reference_state[4])

        pressure = pressure + p_i * (theta_ie - theta_ib) - &
                   rho_tilde * a_tilde * u_i * (sin(theta_ie) - sin(theta_ib)) + &
                   rho_tilde * a_tilde * v_i * (cos(theta_ie) - cos(theta_ib))
      end associate
    end do

    pressure = pressure / (2.0_rk * pi)

  end function get_pressure

  pure function get_density(self, reconstruction_operator) result(density)
    !< Implementation of the rho(P) within the local evolution operator (Eq. 42 in the text)
    class(local_evo_operator_t), intent(in) :: self
    class(abstract_reconstruction_t), pointer, intent(in) :: reconstruction_operator

    real(rk) :: density  !< rho(P)
    integer(ik) :: i, alloc_stat
    integer(ik), dimension(2) :: p_prime_ij
    real(rk), dimension(4) :: p_prime_u_bar !< [rho, u, v, p] at P'

    ! The first two terms in Eq. 42 require evaluating rho and p at P', which needs to be reconstructed
    p_prime_ij = self%find_p_prime_location(self%mach_cone%p_prime_xy)

    allocate(self%p_prime_reconstruction_operator, mold=reconstruction_operator, stat=alloc_stat)
    if(alloc_status /= 0) error stop "Unable to allocate local_evo_operator_t%p_prime_reconstruction_operator"

    self%p_prime_reconstruction_operator%select_cell_to_reconstruct(i=p_prime_ij[1], j=p_prime_ij[2])
    self%p_prime_reconstruction_operator%reconstruct(x=self%mach_cone%p_prime_xy[1], &
                                                     y=self%mach_cone%p_prime_xy[2]
    i = p_prime_ij[1], j = p_prime_ij[2])

    associate(rho_p_prime=>p_prime_u_bar[1], &
              pressure_p_prime=>p_prime_u_bar[4], a_tilde=>self%reference_state[4])
      density = rho_p_prime - pressure_p_prime / a_tilde**2
    end associate

    do i = 1, self%mach_cone%n_intersections

      associate(theta_ie=>self%mach_cone%theta_ie(i), theta_ib=>self%mach_cone%theta_ib(i), &
                p_i=>self%reconstructed_primative_state[4], &
                u_i=>self%reconstructed_primative_state[2], v_i=>self%reconstructed_primative_state[1], &
                rho_tilde=>self%reference_state[1], a_tilde=>self%reference_state[4])

        density = density + (p_i / a_tilde**2) * (theta_ie - theta_ib) - &
                  (rho_tilde / a_tilde) * u_i * (sin(theta_ie) - sin(theta_ib)) + &
                  (rho_tilde / a_tilde) * v_i * (cos(theta_ie) - cos(theta_ib))
      end associate
    end do

    pressure = pressure / (2.0_rk * pi)

    if(allocated(self%p_prime_reconstruction_operator)) then
      deallocate(self%p_prime_reconstruction_operator, stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to dallocate local_evo_operator_t%p_prime_reconstruction_operator"
    end if
  end function get_density

  pure function get_x_velocity(self) result(u)
    !< Implementation of the u(P) within the local evolution operator (Eq 43 in the text)

    class(local_evo_operator_t), intent(in) :: self
    real(rk) :: u
    integer(ik) :: i

    u = 0.0_rk

    do i = 1, self%mach_cone%n_intersections

      associate(theta_ie=>self%mach_cone%theta_ie(i), theta_ib=>self%mach_cone%theta_ib(i), &
                p_i=>self%reconstructed_primative_state[4], &
                u_i=>self%reconstructed_primative_state[2], v_i=>self%reconstructed_primative_state[1], &
                rho_tilde=>self%reference_state[1], a_tilde=>self%reference_state[4])

        u = u - (p_i / (rho_tilde * a_tilde))) * (sin(theta_ie) - sin(theta_ib)) + &
            u_i * (0.5_rk * (theta_ie - theta_ib) + 0.25_rk * (sin(2 * theta_ie) - sin(2 * theta_ib))) - &
            v_i * 0.25_rk * (cos(2 * theta_ie) - cos(2 * theta_ib))

      end associate
    end do

    u = u / pi

  end function get_x_velocity

  pure function get_y_velocity(self) result(v)
    !< Implementation of the v(P) within the local evolution operator (Eq 44 in the text)

    class(local_evo_operator_t), intent(in) :: self
    real(rk) :: v

    v = 0.0_rk

    do i = 1, self%mach_cone%n_intersections

      associate(theta_ie=>self%mach_cone%theta_ie(i), theta_ib=>self%mach_cone%theta_ib(i), &
                p_i=>self%reconstructed_primative_state[4], &
                u_i=>self%reconstructed_primative_state[2], v_i=>self%reconstructed_primative_state[1], &
                rho_tilde=>self%reference_state[1], a_tilde=>self%reference_state[4])

        v = v + (p_i / (rho_tilde * a_tilde))) * (cos(theta_ie) - cos(theta_ib)) + &
            u_i * 0.25_rk * (cos(2 * theta_ie) - cos(2 * theta_ib)) + &
            v_i * (0.5_rk * (theta_ie - theta_ib) + 0.25_rk * (sin(2 * theta_ie) - sin(2 * theta_ib))) - &
            end associate
        end do

        v = v / pi

        end function get_y_velocity
        end module mod_local_evo_operator
