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
    !< Local Evolution Operator (E0)
  contains
    procedure, private :: get_density
    procedure, private, nopass :: get_pressure
    procedure, private, nopass :: get_x_velocity
    procedure, private, nopass :: get_y_velocity
    procedure, public :: evolve_leftright_midpoints
    procedure, public :: evolve_downup_midpoints
    procedure, public :: evolve_corners
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
    nullify(self%grid)
    nullify(self%reference_state)
    nullify(self%conserved_vars)
    nullify(self%reconstruction_operator)
  end subroutine finalize

  pure subroutine evolve_leftright_midpoints(self, reconstructed_state, &
                                             reference_state, evolved_state)
    !< Create Mach cones and evolve the state at all of the left/right midpoints in the domain. Left/right midpoints
    !< are the midpoints defined by vectors pointing to the left and right.

    class(local_evo_operator_t), intent(in) :: self

    real(rk), dimension(:, :, :, :, :), intent(in) :: reconstructed_state
    !< ((rho, u ,v, p), point, node/midpoint, i, j); The reconstructed state of each point P with respect to its parent cell

    real(rk), dimension(:, :, :), intent(in) :: reference_state
    !< ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the left/right edges

    real(rk), dimension(:, :, :), intent(out) :: evolved_state
    !< ((rho,u,v,p), i, j); Reconstructed U at each midpoint on the left/right edges

    ! Locals
    type(mach_cone_geometry_t) :: mach_cone
    !< Mach cone used to provide angles theta_ib and theta_ie

    integer(ik) :: i, j
    integer(ik) :: midpoint_idx  !< used to select the midpoint in the reconstructed_state
    integer(ik) :: point_idx  !< used to select the edge in the reconstructed_state
    integer(ik), dimension(2, 2) :: neighbor_cell_indices
    real(rk), dimension(2, 2, 2) :: midpoint_edge_vectors
    !< ((x,y), head/tail, vector_id); set of vectors that define the midpoint

    real(rk), dimension(4, 2) :: reconstructed_midpoint_state
    !< ((rho, u, v, p), cell_id); the reconstructed state of the midpoint with respect to each cell

    midpoint_idx = 2

    ! Corner/midpoint index convention         Cell Indexing convention
    ! --------------------------------         ------------------------

    !   C----M----C----M----C
    !   |         |         |                             E3
    !   O    x    O    x    O                      N4-----M3----N3
    !   |         |         |                      |            |
    !   C----M----C----M----C                  E4  M4     C     M2  E2
    !   |         |         |                      |            |
    !   O    x    O    x    O                      N1----M1----N2
    !   |         |         |                            E1
    !   C----M----C----M----C

    !  |          |
    !  |  cell 1  |
    !  |  (i,j)   |
    !  C1----M----C2
    !  |  cell 2  |
    !  | (i,j-1)  |
    !  |          |

    ! For left/right midpoints, the edge vectors go left then right.
    ! The neighboring cells are above (i,j) and below (i,j-1)
    ! For quad cells, N - corner, M - midpoint, E - edge

    do concurrent(i=1:ubound(evolved_state, dim=2))
      do concurrent(j=1:ubound(evolved_state, dim=3))

        neighbor_cell_indices = reshape([[i, j], &  ! cell above
                                         [i, j - 1] &  ! cell below
                                         ], shape=[2, 2])

        ! <----M---->  (left and right vectors) these have (x_tail, y_tail) and (x_head, y_head)
        midpoint_edge_vectors = self%grid%get_midpoint_vectors(cell_ij=[i, j], edge='bottom')

        ! reconstructed_state indexing
        !< ((rho, u ,v, p), point, node/midpoint, i, j);

        ! Cell 1: cell above the midpoint -> midpoint point is on the bottom (M1) of the parent cell
        point_idx = 1
        reconstructed_midpoint_state(:, 1) = reconstructed_state(:, point_idx, midpoint_idx, i, j)

        ! Cell 2: cell below the midpoint -> midpoint is on the top (M3) of the parent cell
        point_idx = 3
        reconstructed_midpoint_state(:, 2) = reconstructed_state(:, point_idx, midpoint_idx, i, j - 1)

        mach_cone = new_cone(tau=self%tau, edge_vectors=midpoint_edge_vectors, &
                             reconstructed_state=reconstructed_midpoint_state, &
                             reference_state=reference_state(:, i, j), &
                             cell_indices=neighbor_cell_indices)

        ! Set the evolved state at the midpoint
        evolved_state(:, i, j) = [self%get_density(mach_cone), &    ! rho
                                  self%get_x_velocity(mach_cone), & ! u
                                  self%get_y_velocity(mach_cone), & ! v
                                  self%get_pressure(mach_cone) &    ! p
                                  ]
      end do
    end do
  end subroutine evolve_leftright_midpoints

  pure subroutine evolve_downup_midpoints(self, reconstructed_state, &
                                          reference_state, evolved_state)
    !< Create Mach cones and evolve the state at all of the down/up midpoints in the domain. Left/right midpoints
    !< are the midpoints defined by vectors pointing to the down and up.

    class(local_evo_operator_t), intent(in) :: self

    real(rk), dimension(:, :, :, :, :), intent(in) :: reconstructed_state
    !< ((rho, u ,v, p), point, node/midpoint, i, j); The reconstructed state of each point P with respect to its parent cell

    real(rk), dimension(:, :, :), intent(in) :: reference_state
    !< ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the left/right edges

    real(rk), dimension(:, :, :), intent(out) :: evolved_state
    !< ((rho,u,v,p), i, j); Evolved U at each midpoint on the left/right edges

    ! Locals
    type(mach_cone_geometry_t) :: mach_cone
    !< Mach cone used to provide angles theta_ib and theta_ie

    integer(ik) :: i, j
    integer(ik) :: midpoint_idx  !< used to select the midpoint in the reconstructed_state
    integer(ik) :: point_idx  !< used to select the edge in the reconstructed_state
    integer(ik), dimension(2, 2) :: neighbor_cell_indices
    real(rk), dimension(2, 2, 2) :: midpoint_edge_vectors
    !< ((x,y), head/tail, vector_id); set of vectors that define the midpoint

    real(rk), dimension(4, 2) :: reconstructed_midpoint_state
    !< ((rho, u, v, p), cell_id); the reconstructed state of the midpoint with respect to each cell

    midpoint_idx = 2

    ! Corner/midpoint index convention         Cell Indexing convention
    ! --------------------------------         ------------------------

    !   C----M----C----M----C
    !   |         |         |                             E3
    !   O    x    O    x    O                      N4-----M3----N3
    !   |         |         |                      |            |
    !   C----M----C----M----C                  E4  M4     C     M2  E2
    !   |         |         |                      |            |
    !   O    x    O    x    O                      N1----M1----N2
    !   |         |         |                            E1
    !   C----M----C----M----C

    !  |         C2         |
    !  |  cell 1 |  cell 2  |
    !  | (i-1,j) |  (i,j)   |
    !  |         M          |
    !  |         |          |
    !  |         |          |
    !  |         C1         |

    ! For down/up midpoints, the edge vectors go down then up.
    ! The neighboring cells are left (i, j) and right (i-1, j)
    ! For quad cells, N - corner, M - midpoint, E - edge

    do concurrent(i=1:ubound(evolved_state, dim=2))
      do concurrent(j=1:ubound(evolved_state, dim=3))

        ! cell ordering is 1) left, 2) right
        neighbor_cell_indices = reshape([[i - 1, j], &  ! cell left
                                         [i, j] &  ! cell right
                                         ], shape=[2, 2])

        midpoint_edge_vectors = self%grid%get_midpoint_vectors(cell_ij=[i, j], edge='left')

        ! Cell 1: cell to the left -> midpoint point is on the left (M4) of the parent cell
        point_idx = 4
        reconstructed_midpoint_state(:, 1) = reconstructed_state(:, point_idx, midpoint_idx, i - 1, j)

        ! Cell 2: cell to the right -> midpoint point is on the right (M2) of the parent cell
        point_idx = 2
        reconstructed_midpoint_state(:, 2) = reconstructed_state(:, point_idx, midpoint_idx, i, j)

        mach_cone = new_cone(tau=self%tau, edge_vectors=midpoint_edge_vectors, &
                             reconstructed_state=reconstructed_midpoint_state, &
                             reference_state=reference_state(:, i, j), &
                             cell_indices=neighbor_cell_indices)

        ! Set the evolved state at the midpoint
        evolved_state(:, i, j) = [self%get_density(mach_cone), &    ! rho
                                  self%get_x_velocity(mach_cone), & ! u
                                  self%get_y_velocity(mach_cone), & ! v
                                  self%get_pressure(mach_cone) &    ! p
                                  ]
      end do
    end do
  end subroutine evolve_downup_midpoints

  pure subroutine evolve_corners(self, reconstructed_state, &
                                 reference_state, evolved_state)
    !< Create Mach cones and evolve the state at all of the down/up midpoints in the domain. Left/right midpoints
    !< are the midpoints defined by vectors pointing to the down and up.

    class(local_evo_operator_t), intent(in) :: self

    real(rk), dimension(:, :, :, :, :), intent(in) :: reconstructed_state
    !< ((rho, u ,v, p), point, node/midpoint, i, j); The reconstructed state of each point P with respect to its parent cell

    real(rk), dimension(:, :, :), intent(in) :: reference_state
    !< ((rho,u,v,p), i, j); Reference state (tilde) at each corner

    real(rk), dimension(:, :, :), intent(out) :: evolved_state
    !< ((rho,u,v,p), i, j); Reconstructed U at each corner

    ! Locals
    type(mach_cone_geometry_t) :: mach_cone
    !< Mach cone used to provide angles theta_ib and theta_ie

    integer(ik) :: i, j
    integer(ik) :: corner_idx  !< used to select the corner in the reconstructed_state
    integer(ik) :: point_idx  !< used to select the point in the reconstructed_state
    integer(ik), dimension(2, 4) :: neighbor_cell_indices
    real(rk), dimension(2, 2, 4) :: corner_edge_vectors
    !< ((x,y), head/tail, vector_id); set of vectors that define the corner

    real(rk), dimension(4, 4) :: reconstructed_corner_state
    !< ((rho, u, v, p), cell_id); the reconstructed state of the corner with respect to each cell

    corner_idx = 1

    ! Corner/midpoint index convention         Cell Indexing convention
    ! --------------------------------         ------------------------

    !   C----M----C----M----C
    !   |         |         |                             E3
    !   O    x    O    x    O                      N4-----M3----N3
    !   |         |         |                      |            |
    !   C----M----C----M----C                  E4  M4     C     M2  E2
    !   |         |         |                      |            |
    !   O    x    O    x    O                      N1----M1----N2
    !   |         |         |                            E1
    !   C----M----C----M----C

    !           C3
    !   cell 4  |   cell 3
    !   (i-1,j) |   (i,j)
    !           |
    ! C4-------C0--------C2
    !           |
    ! (i-1,j-1) |   (i,j-1)
    !   cell 1  |   cell 2
    !          C1

    ! For down/up midpoints, the edge vectors go down then up.
    ! The neighboring cells are left (i, j) and right (i-1, j)
    ! For quad cells, N - corner, M - midpoint, E - edge

    do concurrent(i=1:ubound(evolved_state, dim=2))
      do concurrent(j=1:ubound(evolved_state, dim=3))

        ! cell ordering is 1) lower left, 2) lower right, 3) upper right, 4) upper left
        neighbor_cell_indices = reshape([[i - 1, j - 1], &  ! lower left
                                         [i, j - 1], &  ! lower right
                                         [i, j], &  ! upper right
                                         [i - 1, j] &  ! upper left
                                         ], shape=[2, 4])

        corner_edge_vectors = self%grid%get_corner_vectors(cell_ij=[i, j], corner='lowerleft')

        ! reconstructed_state indexing
        !< ((rho, u ,v, p), point, node/midpoint, i, j);

        ! Cell 1: lower left cell -> corner is in the upper right (N3) of its parent cell
        point_idx = 3
        reconstructed_corner_state(:, 1) = reconstructed_state(:, point_idx, corner_idx, i - 1, j - 1)

        ! Cell 2: lower right cell -> corner is in the upper left (N4) of its parent cell
        point_idx = 4
        reconstructed_corner_state(:, 2) = reconstructed_state(:, point_idx, corner_idx, i, j - 1)

        ! Cell 3: upper right cell-> corner is in the lower left (N1) of its parent cell
        point_idx = 1
        reconstructed_corner_state(:, 3) = reconstructed_state(:, point_idx, corner_idx, i, j)

        ! Cell 4: upper left cell -> corner is in the lower right (N2) of its parent cell
        point_idx = 2
        reconstructed_corner_state(:, 4) = reconstructed_state(:, point_idx, corner_idx, i - 1, j)

        mach_cone = new_cone(tau=self%tau, edge_vectors=corner_edge_vectors, &
                             reconstructed_state=reconstructed_corner_state, &
                             reference_state=reference_state(:, i, j), &
                             cell_indices=neighbor_cell_indices)

        ! Set the evolved state at the midpoint
        evolved_state(:, i, j) = [self%get_density(mach_cone), &    ! rho
                                  self%get_x_velocity(mach_cone), & ! u
                                  self%get_y_velocity(mach_cone), & ! v
                                  self%get_pressure(mach_cone) &    ! p
                                  ]
      end do
    end do
  end subroutine evolve_corners

  pure function get_pressure(mach_cone) result(pressure)
    !< Implementation of p(P) within the local evolution operator (Eq 45 in the text)

    class(mach_cone_geometry_t), intent(in) :: mach_cone
    real(rk) :: pressure
    integer(ik) :: i, cell

    pressure = 0.0_rk

    do cell = 1, mach_cone%n_neighbor_cells
      do i = 1, mach_cone%n_intersections(cell)

        associate(theta_ie=>mach_cone%theta_ie(cell, i), &
                  theta_ib=>mach_cone%theta_ib(cell, i), &
                  u_i=>mach_cone%cell_conserved_vars(2, i, cell), &
                  v_i=>mach_cone%cell_conserved_vars(3, i, cell), &
                  p_i=>mach_cone%cell_conserved_vars(4, i, cell), &
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
    real(rk), dimension(4) :: p_prime_u_bar !< [rho, u, v, p] at P'

    density = 0.0_rk

    do cell = 1, mach_cone%n_neighbor_cells
      do i = 1, mach_cone%n_intersections(cell)

        associate(theta_ie=>mach_cone%theta_ie(cell, i), &
                  theta_ib=>mach_cone%theta_ib(cell, i), &
                  u_i=>mach_cone%cell_conserved_vars(2, i, cell), &
                  v_i=>mach_cone%cell_conserved_vars(3, i, cell), &
                  p_i=>mach_cone%cell_conserved_vars(4, i, cell), &
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

    do cell = 1, mach_cone%n_neighbor_cells
      do i = 1, mach_cone%n_intersections(cell)

        associate(theta_ie=>mach_cone%theta_ie(cell, i), &
                  theta_ib=>mach_cone%theta_ib(cell, i), &
                  u_i=>mach_cone%cell_conserved_vars(2, i, cell), &
                  v_i=>mach_cone%cell_conserved_vars(3, i, cell), &
                  p_i=>mach_cone%cell_conserved_vars(4, i, cell), &
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

    do cell = 1, mach_cone%n_neighbor_cells
      do i = 1, mach_cone%n_intersections(cell)

        associate(theta_ie=>mach_cone%theta_ie(cell, i), &
                  theta_ib=>mach_cone%theta_ib(cell, i), &
                  u_i=>mach_cone%cell_conserved_vars(2, i, cell), &
                  v_i=>mach_cone%cell_conserved_vars(3, i, cell), &
                  p_i=>mach_cone%cell_conserved_vars(4, i, cell), &
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
