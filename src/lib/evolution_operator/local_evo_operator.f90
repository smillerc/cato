module mod_local_evo_operator
  use iso_fortran_env, only: ik => int32, rk => real64
  use math_constants, only: pi
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_input, only: input_t
  use mod_mach_cone_geometry, only: mach_cone_geometry_t, new_cone
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
    procedure, private, nopass :: get_pressure
    procedure, private :: get_density
    procedure, private, nopass :: get_x_velocity
    procedure, private, nopass :: get_y_velocity
    procedure, private :: evolve_leftright_midpoints
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
    real(rk), dimension(:, :, :), intent(in), target :: conserved_vars
    ! real(rk), dimension(:, :, :), intent(in), target :: reference_state

    ! ! allocate(operator)
    operator%name = "FVLEG"
    call operator%set_tau(input%tau)
    operator%grid => grid
    operator%conserved_vars => conserved_vars
    ! operator%reference_state => reference_state

    ! p_prime_reconstruction_operator = recon_factory

  end function constructor

  subroutine finalize(self)
    !< Cleanup the operator type
    type(local_evo_operator_t), intent(inout) :: self
    ! deallocate(self%mach_cone)
    nullify(self%grid)
    nullify(self%reference_state)
    nullify(self%conserved_vars)
    nullify(self%reconstruction_operator)
  end subroutine finalize

!  ! Locals
!   integer(ik) :: i, j
!   integer(ik), dimension(2,2,2) :: midpoint_edge_vectors
!     !< ((x,y), head/tail, vector_id); set of vectors that define the midpoint

!   allocate(mach_cone_geometry_t :: self%mach_cone, stat=alloc_stat)
!   if (alloc_stat /= 0) error stop 'Unable to allocate local_evo_operator_t%mach_cone'
!   call self%mach_cone%initialize(tau=self%tau)

!   ! left/right midpoints
!   ! order is left then right
!   do i = 1, 2
!     do j = 1, 2

!       ! left
!       ! right
!       midpoint_neighbor_cells = reshape([[i-1,j], [i, j]], shape=[2,2])

!     end do
!   end do

!   ! up/down midpoints
!   ! order is down then up. I know this is opposite of the name, but who says down/up?

!   ! left/right
!   midpoint_edge_vectors = self%grid%get_midpoint_vectors(cell_ij=[i,j], edge='bottom')

!   call self%mach_cone%set_p_state(xy=, &
!                                   edge_vectors=midpoint_edge_vectors, &
!                                   reconstructed_state=self%, &
!                                   reference_state=)

!   ! midpoint_neighbor_cells = reshape([[i,j-1], [i, j]], shape=[2,2])

!   ! corners
!   ! the order of neighhor cell indices is important here,
!   ! it starts at the lower left-most and
!   ! goes counter-clockwise around the common corner node
!   corner_neighbor_cells = reshape([[i - 1,j - 1], [i, j-1], &
!                                    [i, j], [i - 1, j]], shape=[2,4])

  pure subroutine evolve_leftright_midpoints(self, reconstructed_state, &
                                             leftright_midpoints_reference_state, evolved_leftright_midpoints_state)
    !< Create Mach cones and evolve the state at all of the left/right midpoints in the domain. Left/right midpoints
    !< are the midpoints defined by vectors pointing to the left and right.

    class(local_evo_operator_t), intent(in) :: self

    real(rk), dimension(:, :, :, :, :), intent(in) :: reconstructed_state
    !< ((rho, u ,v, p), point, node/midpoint, i, j); The reconstructed state of each point P with respect to its parent cell

    real(rk), dimension(:, :, :), intent(in) :: leftright_midpoints_reference_state
    !< ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the left/right edges

    real(rk), dimension(:, :, :), intent(inout) :: evolved_leftright_midpoints_state
    !< ((rho,u,v,p), i, j); Reconstructed U at each midpoint on the left/right edges

    ! Locals
    type(mach_cone_geometry_t) :: mach_cone
    !< Mach cone used to provide angles theta_ib and theta_ie

    integer(ik) :: i, j, alloc_stat
    integer(ik) :: midpoint_idx  !< used to select the midpoint in the reconstructed_state
    integer(ik) :: edge_idx  !< used to select the edge in the reconstructed_state
    integer(ik), dimension(2, 2) :: neighbor_cell_indices
    real(rk), dimension(2, 2, 2) :: midpoint_edge_vectors
    !< ((x,y), head/tail, vector_id); set of vectors that define the midpoint

    real(rk), dimension(4, 2) :: reconstructed_midpoint_state
    !< ((rho, u, v, p), cell_id); the reconstructed state of the midpoint with respect to each cell

    ! allocate(mach_cone_geometry_t :: mach_cone, stat=alloc_stat)
    ! if (alloc_stat /= 0) error stop 'Unable to allocate mach_cone in local_evo_operator_t%evolve_leftright_midpoints'
    ! call mach_cone%initialize(tau=self%tau)

    midpoint_idx = 2

    ! Corner/midpoint index convention         Cell Indexing convention
    ! --------------------------------         ------------------------
    !
    !   C----M----C----M----C
    !   |         |         |                             E3
    !   O    x    O    x    O                      N4-----M3----N3
    !   |         |         |                      |            |
    !   C----M----C----M----C                  E4  M4     C     M2  E2
    !   |         |         |                      |            |
    !   O    x    O    x    O                      N1----M1----N2
    !   |         |         |                            E1
    !   C----M----C----M----C

    ! For left/right midpoints, the edge vectors go left then right.
    ! The neighboring cells are above (i,j) and below (i,j-1)
    ! For quad cells, N - corner, M - midpoint, E - edge

    do concurrent(i=1:ubound(evolved_leftright_midpoints_state, dim=2))
      do concurrent(j=1:ubound(evolved_leftright_midpoints_state, dim=3))

        neighbor_cell_indices = reshape([[i, j],[i, j + 1]], shape=[2, 2])

        ! <----M---->  (left and right vectors) these have (x_tail, y_tail) and (x_head, y_head)
        midpoint_edge_vectors = self%grid%get_midpoint_vectors(cell_ij=[i, j], edge='bottom')

        ! cell below the midpoint
        edge_idx = 3
        reconstructed_midpoint_state(:, 1) = reconstructed_state(:, edge_idx, midpoint_idx, i, j - 1)

        ! cell above the midpoint
        edge_idx = 1
        reconstructed_midpoint_state(:, 2) = reconstructed_state(:, edge_idx, midpoint_idx, i, j)

        mach_cone = new_cone(tau=self%tau, edge_vectors=midpoint_edge_vectors, &
                             reconstructed_state=reconstructed_midpoint_state, &
                             reference_state=leftright_midpoints_reference_state(:, i, j), &
                             cell_indices=neighbor_cell_indices)

        ! Set the evolved state at the midpoint
        evolved_leftright_midpoints_state(:, i, j) = [self%get_density(mach_cone), &
                                                      self%get_x_velocity(mach_cone), &
                                                      self%get_y_velocity(mach_cone), &
                                                      self%get_pressure(mach_cone)]
      end do
    end do
  end subroutine

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
    integer(ik), dimension(2, 4) :: cell_group

    ! cell_group = 0
    ! ! midpoints
    ! if (loc == 'k,c') then
    !   select case(edge)
    !   case(1)
    !     cell_group(:,1:2) = reshape([[i,j], [i,j-1]],[2,2])
    !   case(2)
    !     cell_group(:,1:2) = reshape([[i,j], [i+1,j]],[2,2])
    !   case(3)
    !     cell_group(:,1:2) = reshape([[i,j], [i,j+1]],[2,2])
    !   case(4)
    !     cell_group(:,1:2) = reshape([[i,j], [i-1,j]],[2,2])
    !   case default
    !     error stop 'Invalid edge for midpoint in local_evo_operator_t%evolve'
    !   end select
    ! end if

    ! ! corners
    ! if (loc == 'k,1' .or. loc == 'k,2') then
    !   select case(edge)
    !   case(1)
    !     cell_group(:,1:4) = reshape([[i,j], [i,j-1], &
    !                                  [i-1,j], [i-1,j-1]], [2,4])
    !   case(2)
    !     cell_group(:,1:4) = reshape([[i,j], [i,j-1], &
    !                                  [i+1,j], [i+1,j-1]], [2,4])
    !   case(3)
    !     cell_group(:,1:4) = reshape([[i,j], [i,j+1], &
    !                                  [i+1,j], [i+1,j+1]], [2,4])
    !   case(4)
    !     cell_group(:,1:4) = reshape([[i,j], [i,j+1], &
    !                                  [i-1,j], [i-1,j+1]], [2,4])
    !   case default
    !     error stop 'Invalid edge for midpoint in local_evo_operator_t%evolve'
    !   end select
    ! end if

    ! if (allocated(self%mach_cone)) deallocate(self%mach_cone, stat=alloc_stat)
    ! if (alloc_stat /= 0) error stop 'Unable to deallocate local_evo_operator_t%mach_cone inside of local_evo_operator_t%evolve'

    ! allocate(mach_cone_geometry_t :: self%mach_cone)

    ! select case(loc)
    ! case("k,1") ! 1st corner
    !   call self%mach_cone%initialize(cell_group=cell_group,reference_state=self%reference_state())
    ! ! case("k,c") ! midpoint
    ! !   self%mach_cone = mach_cone_geometry_t(edge_points=self%grid%cell_edge_vectors(:,2,edge,i,j))
    ! ! case("k,2") ! 2nd corner
    ! !   self%mach_cone = mach_cone_geometry_t(edge_points=self%grid%cell_edge_vectors(:,3,edge,i,j))
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

  pure function get_pressure(mach_cone) result(pressure)
    !< Implementation of p(P) within the local evolution operator (Eq 45 in the text)

    class(mach_cone_geometry_t), intent(in) :: mach_cone
    real(rk) :: pressure
    integer(ik) :: i, cell

    pressure = 0.0_rk

    do cell = 1, mach_cone%neighbor_cells
      do i = 1, mach_cone%n_intersections(cell)

        associate(theta_ie=>mach_cone%theta_ie(cell, i), &
                  theta_ib=>mach_cone%theta_ib(cell, i), &
                  p_i=>mach_cone%cell_conserved_vars(4, i, cell), &
                  u_i=>mach_cone%cell_conserved_vars(2, i, cell), &
                  v_i=>mach_cone%cell_conserved_vars(1, i, cell), &
                  rho_tilde=>mach_cone%reference_state(1), &
                  a_tilde=>mach_cone%reference_state(4))

          pressure = pressure + p_i * (theta_ie - theta_ib) - &
                     rho_tilde * a_tilde * u_i * (sin(theta_ie) - sin(theta_ib)) + &
                     rho_tilde * a_tilde * v_i * (cos(theta_ie) - cos(theta_ib))
        end associate
      end do
    end do

    pressure = pressure / (2.0_rk * pi)

  end function get_pressure

  pure function get_density(self, mach_cone) result(density)
    !< Implementation of rho(P) within the local evolution operator (Eq. 42 in the text)
    class(local_evo_operator_t), intent(in) :: self
    class(mach_cone_geometry_t), intent(in) :: mach_cone

    real(rk) :: density  !< rho(P)
    integer(ik) :: i, cell
    integer(ik), dimension(2) :: p_prime_ij
    real(rk), dimension(4) :: p_prime_u_bar !< [rho, u, v, p] at P'
    logical :: p_prime_found

    density = 0.0_rk

    do cell = 1, mach_cone%neighbor_cells
      do i = 1, mach_cone%n_intersections(cell)

        associate(theta_ie=>mach_cone%theta_ie(cell, i), &
                  theta_ib=>mach_cone%theta_ib(cell, i), &
                  p_i=>mach_cone%cell_conserved_vars(4, i, cell), &
                  u_i=>mach_cone%cell_conserved_vars(2, i, cell), &
                  v_i=>mach_cone%cell_conserved_vars(1, i, cell), &
                  rho_tilde=>mach_cone%reference_state(1), &
                  a_tilde=>mach_cone%reference_state(4))

          density = density + (p_i / a_tilde**2) * (theta_ie - theta_ib) - &
                    (rho_tilde / a_tilde) * u_i * (sin(theta_ie) - sin(theta_ib)) + &
                    (rho_tilde / a_tilde) * v_i * (cos(theta_ie) - cos(theta_ib))
        end associate
      end do
    end do

    density = density / (2.0_rk * pi)

    ! Reconstruct at P' to get rho(P') and p(P')
    p_prime_u_bar = self%reconstruction_operator%reconstruct_point(xy=mach_cone%p_prime_xy, &
                                                                   cell_ij=mach_cone%p_prime_ij)

    associate(rho_p_prime=>p_prime_u_bar(1), pressure_p_prime=>p_prime_u_bar(4), &
              a_tilde=>mach_cone%reference_state(4))
      density = density + rho_p_prime - (pressure_p_prime / a_tilde**2)
    end associate

  end function get_density

  pure function get_x_velocity(mach_cone) result(u)
    !< Implementation of u(P) within the local evolution operator (Eq 43 in the text)

    class(mach_cone_geometry_t), intent(in) :: mach_cone
    real(rk) :: u !< x velocity at P
    integer(ik) :: i
    integer(ik) :: cell  !< cell index within the mach cone

    u = 0.0_rk

    do cell = 1, mach_cone%neighbor_cells
      do i = 1, mach_cone%n_intersections(cell)

        associate(theta_ie=>mach_cone%theta_ie(cell, i), &
                  theta_ib=>mach_cone%theta_ib(cell, i), &
                  p_i=>mach_cone%cell_conserved_vars(4, i, cell), &
                  u_i=>mach_cone%cell_conserved_vars(2, i, cell), &
                  v_i=>mach_cone%cell_conserved_vars(1, i, cell), &
                  rho_tilde=>mach_cone%reference_state(1), &
                  a_tilde=>mach_cone%reference_state(4))

          u = u - (p_i / (rho_tilde * a_tilde)) * (sin(theta_ie) - sin(theta_ib)) + &
              u_i * (0.5_rk * (theta_ie - theta_ib) + 0.25_rk * (sin(2 * theta_ie) - sin(2 * theta_ib))) - &
              v_i * 0.25_rk * (cos(2 * theta_ie) - cos(2 * theta_ib))

        end associate
      end do
    end do

    u = u / pi

  end function get_x_velocity

  pure function get_y_velocity(mach_cone) result(v)
    !< Implementation of v(P) within the local evolution operator (Eq 44 in the text)

    class(mach_cone_geometry_t), intent(in) :: mach_cone
    real(rk) :: v !< y velocity at P
    integer(ik) :: i
    integer(ik) :: cell  !< cell index within the mach cone

    v = 0.0_rk

    do cell = 1, mach_cone%neighbor_cells
      do i = 1, mach_cone%n_intersections(cell)

        associate(theta_ie=>mach_cone%theta_ie(cell, i), &
                  theta_ib=>mach_cone%theta_ib(cell, i), &
                  p_i=>mach_cone%cell_conserved_vars(4, i, cell), &
                  u_i=>mach_cone%cell_conserved_vars(2, i, cell), &
                  v_i=>mach_cone%cell_conserved_vars(1, i, cell), &
                  rho_tilde=>mach_cone%reference_state(1), &
                  a_tilde=>mach_cone%reference_state(4))

          v = v + (p_i / (rho_tilde * a_tilde)) * (cos(theta_ie) - cos(theta_ib)) + &
              u_i * 0.25_rk * (cos(2 * theta_ie) - cos(2 * theta_ib)) + &
              v_i * (0.5_rk * (theta_ie - theta_ib) + 0.25_rk * (sin(2 * theta_ie) - sin(2 * theta_ib)))
        end associate
      end do
    end do

    v = v / pi
  end function get_y_velocity

end module mod_local_evo_operator
