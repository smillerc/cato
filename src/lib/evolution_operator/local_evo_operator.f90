module mod_local_evo_operator
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_globals, only: debug_print

  use, intrinsic :: ieee_arithmetic
  use mod_floating_point_utils, only: near_zero, equal
  use math_constants, only: pi, rad2deg
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_input, only: input_t
  use mod_cone, only: cone_t, new_cone
  use mod_grid, only: grid_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_reconstruction_factory, only: reconstruction_factory

  implicit none

  private
  public :: local_evo_operator_t, get_density, get_pressure, get_x_velocity, get_y_velocity

  type, extends(abstract_evo_operator_t) :: local_evo_operator_t
    !< Local Evolution Operator (E0)
  contains
    procedure, public :: initialize
    procedure, private, nopass :: get_y_velocity
    procedure, public :: evolve_leftright_midpoints
    procedure, public :: evolve_downup_midpoints
    procedure, public :: evolve_corners
    procedure, public :: copy
    final :: finalize
  end type local_evo_operator_t

contains

  subroutine initialize(self, input, grid_target, recon_operator_target, reconstructed_state_target, lbounds)
    !< Constructor for the FVLEG operator
    class(local_evo_operator_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    class(grid_t), intent(in), target :: grid_target
    class(abstract_reconstruction_t), intent(in), target :: recon_operator_target
    integer(ik), dimension(5), intent(in) :: lbounds
    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):, &
                        lbounds(4):, lbounds(5):), intent(in), target :: reconstructed_state_target

    self%name = "FVLEG"
    call self%set_tau(input%tau)
    self%grid => grid_target
    self%reconstructed_state => reconstructed_state_target
    self%reconstruction_operator => recon_operator_target

  end subroutine

  subroutine finalize(self)
    !< Cleanup the operator type
    type(local_evo_operator_t), intent(inout) :: self

    call debug_print('Calling local_evo_operator_t%finalize()', __FILE__, __LINE__)

    if(allocated(self%name)) deallocate(self%name)
    if(associated(self%grid)) nullify(self%grid)
    if(associated(self%reconstructed_state)) nullify(self%reconstructed_state)
    if(associated(self%reconstruction_operator)) nullify(self%reconstruction_operator)
  end subroutine finalize

  subroutine copy(out_evo, in_evo)
    class(abstract_evo_operator_t), intent(in) :: in_evo
    class(local_evo_operator_t), intent(inout) :: out_evo

    call debug_print('Calling local_evo_operator_t%copy()', __FILE__, __LINE__)

    if(associated(out_evo%grid)) nullify(out_evo%grid)
    ! allocate(out_evo%grid, source=in_evo%grid)
    out_evo%grid => in_evo%grid

    if(associated(out_evo%reconstructed_state)) nullify(out_evo%reconstructed_state)
    ! allocate(out_evo%reconstructed_state, source=in_evo%reconstructed_state)
    out_evo%reconstructed_state => in_evo%reconstructed_state

    if(associated(out_evo%reconstruction_operator)) nullify(out_evo%reconstruction_operator)
    ! allocate(out_evo%reconstruction_operator, source=in_evo%reconstruction_operator)
    out_evo%reconstruction_operator => in_evo%reconstruction_operator

    if(allocated(out_evo%name)) deallocate(out_evo%name)
    allocate(out_evo%name, source=in_evo%name)

  end subroutine

  pure subroutine evolve_leftright_midpoints(self, reference_state, evolved_state)
    !< Create Mach cones and evolve the state at all of the left/right midpoints in the domain. Left/right midpoints
    !< are the midpoints defined by vectors pointing to the left and right.

    class(local_evo_operator_t), intent(in) :: self

    ! real(rk), dimension(:, :, :, 0:, 0:), intent(in) :: reconstructed_state
    !< ((rho, u ,v, p), point, node/midpoint, i, j); The reconstructed state of each point P with respect to its parent cell

    real(rk), dimension(:, :, :), intent(in) :: reference_state
    ! < ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the left/right edges

    real(rk), dimension(:, :, :), intent(out) :: evolved_state
    !< ((rho,u,v,p), i, j); Reconstructed U at each midpoint on the left/right edges

    ! real(rk), dimension(:, 0:, 0:), intent(in) :: conserved_vars

    ! Locals
    type(cone_t) :: mach_cone
    !< Mach cone used to provide angles theta_ib and theta_ie

    integer(ik) :: i, j
    integer(ik) :: midpoint_idx  !< used to select the midpoint in the reconstructed_state
    integer(ik) :: point_idx  !< used to select the edge in the reconstructed_state
    integer(ik), dimension(2, 2) :: neighbor_cell_indices
    real(rk), dimension(2, 2, 2) :: midpoint_edge_vectors
    !< ((x,y), head/tail, vector_id); set of vectors that define the midpoint

    real(rk), dimension(4, 2) :: reconstructed_midpoint_state
    !< ((rho, u, v, p), cell_id); the reconstructed state of the midpoint with respect to each cell

    integer(ik) :: ilo, ihi, jlo, jhi

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

    ! We only reconstruct the real domain
    ilo = lbound(reference_state, dim=2) ! in the i-direction, # left/right midpoints = # cells
    ihi = ubound(reference_state, dim=2) ! in the i-direction, # left/right midpoints = # cells
    jlo = lbound(reference_state, dim=3) ! in the j-direction, # left/right midpoints = # nodes
    jhi = ubound(reference_state, dim=3) ! in the j-direction, # left/right midpoints = # nodes

    do concurrent(i=ilo:ihi)
      do concurrent(j=jlo:jhi)

        neighbor_cell_indices = reshape([[i, j], &  ! cell above
                                         [i, j - 1] &  ! cell below
                                         ], shape=[2, 2])

        ! <----M---->  (left and right vectors) these have (x_tail, y_tail) and (x_head, y_head)
        midpoint_edge_vectors = self%grid%get_midpoint_vectors(cell_ij=[i, j], edge='bottom')

        ! reconstructed_state indexing
        !< ((rho, u ,v, p), point, node/midpoint, i, j);

        ! Cell 1: cell above the midpoint -> midpoint point is on the bottom (M1) of the parent cell
        point_idx = 1
        reconstructed_midpoint_state(:, 1) = self%reconstructed_state(:, point_idx, midpoint_idx, i, j)

        ! Cell 2: cell below the midpoint -> midpoint is on the top (M3) of the parent cell
        point_idx = 3
        reconstructed_midpoint_state(:, 2) = self%reconstructed_state(:, point_idx, midpoint_idx, i, j - 1)
        mach_cone = new_cone(tau=self%tau, edge_vectors=midpoint_edge_vectors, &
                             reconstructed_state=reconstructed_midpoint_state, &
                             reference_state=reference_state(:, i, j), &
                             cell_indices=neighbor_cell_indices)

        ! Set the evolved state at the midpoint
        evolved_state(:, i, j) = [get_density(self, mach_cone), &    ! rho
                                  get_x_velocity(mach_cone), & ! u
                                  get_y_velocity(mach_cone), & ! v
                                  get_pressure(mach_cone) &    ! p
                                  ]

        ! if (i == 50 .and. j == 50) then
        ! print *, mach_cone

        ! error stop
        ! end if
      end do
    end do
  end subroutine evolve_leftright_midpoints

  pure subroutine evolve_downup_midpoints(self, reference_state, evolved_state)
    !< Create Mach cones and evolve the state at all of the down/up midpoints in the domain. Left/right midpoints
    !< are the midpoints defined by vectors pointing to the down and up.

    class(local_evo_operator_t), intent(in) :: self

    ! real(rk), dimension(:, :, :, 0:, 0:), intent(in) :: reconstructed_state
    !< ((rho, u ,v, p), point, node/midpoint, i, j);
    !< The reconstructed state of each point P with respect to its parent cell

    real(rk), dimension(:, :, :), intent(in) :: reference_state
    ! < ((rho,u,v,p), i, j); Reference state (tilde) at each midpoint on the left/right edges

    real(rk), dimension(:, :, :), intent(out) :: evolved_state
    !< ((rho,u,v,p), i, j); Evolved U at each midpoint on the left/right edges

    ! real(rk), dimension(:, 0:, 0:), intent(in) :: conserved_vars

    ! Locals
    type(cone_t) :: mach_cone
    !< Mach cone used to provide angles theta_ib and theta_ie

    integer(ik) :: i, j
    integer(ik) :: midpoint_idx  !< used to select the midpoint in the reconstructed_state
    integer(ik) :: point_idx  !< used to select the edge in the reconstructed_state
    integer(ik), dimension(2, 2) :: neighbor_cell_indices
    real(rk), dimension(2, 2, 2) :: midpoint_edge_vectors !< ((x,y), (tail,head), (vector1:vector2))

    real(rk), dimension(4, 2) :: reconstructed_midpoint_state
    !< ((rho, u, v, p), cell_id); the reconstructed state of the midpoint with respect to each cell

    integer(ik) :: ilo, ihi, jlo, jhi

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

    ilo = lbound(reference_state, dim=2)  ! in the i-direction, # down/up midpoints = # nodes
    ihi = ubound(reference_state, dim=2)  ! in the i-direction, # down/up midpoints = # nodes
    jlo = lbound(reference_state, dim=3)  ! in the j-direction, # down/up midpoints = # cells
    jhi = ubound(reference_state, dim=3)  ! in the j-direction, # down/up midpoints = # cells

    do concurrent(i=ilo:ihi)
      do concurrent(j=jlo:jhi)

        ! cell ordering is 1) left, 2) right
        neighbor_cell_indices = reshape([[i - 1, j], &  ! cell left
                                         [i, j] &  ! cell right
                                         ], shape=[2, 2])

        midpoint_edge_vectors = self%grid%get_midpoint_vectors(cell_ij=[i, j], edge='left')

        ! Cell 1: cell to the left -> midpoint point is on the left (M4) of the parent cell
        point_idx = 4
        reconstructed_midpoint_state(:, 1) = self%reconstructed_state(:, point_idx, midpoint_idx, i - 1, j)

        ! Cell 2: cell to the right -> midpoint point is on the right (M2) of the parent cell
        point_idx = 2
        reconstructed_midpoint_state(:, 2) = self%reconstructed_state(:, point_idx, midpoint_idx, i, j)
        mach_cone = new_cone(tau=self%tau, edge_vectors=midpoint_edge_vectors, &
                             reconstructed_state=reconstructed_midpoint_state, &
                             reference_state=reference_state(:, i, j), &
                             cell_indices=neighbor_cell_indices)

        ! Set the evolved state at the midpoint
        evolved_state(:, i, j) = [get_density(self, mach_cone), &    ! rho
                                  get_x_velocity(mach_cone), & ! u
                                  get_y_velocity(mach_cone), & ! v
                                  get_pressure(mach_cone) &    ! p
                                  ]

      end do
    end do
  end subroutine evolve_downup_midpoints

  pure subroutine evolve_corners(self, reference_state, evolved_state)
    !< Create Mach cones and evolve the state at all of the down/up midpoints in the domain. Left/right midpoints
    !< are the midpoints defined by vectors pointing to the down and up.

    class(local_evo_operator_t), intent(in) :: self

    ! real(rk), dimension(:, :, :, 0:, 0:), intent(in) :: reconstructed_state
    ! !< ((rho, u ,v, p), point, node/midpoint, i, j); The reconstructed state of each point P with respect to its parent cell

    real(rk), dimension(:, :, :), intent(in) :: reference_state
    !< ((rho,u,v,p), i, j); Reference state (tilde) at each corner

    real(rk), dimension(:, :, :), intent(out) :: evolved_state
    !< ((rho,u,v,p), i, j); Reconstructed U at each corner

    ! real(rk), dimension(:, 0:, 0:), intent(in) :: conserved_vars

    ! Locals
    type(cone_t) :: mach_cone
    !< Mach cone used to provide angles theta_ib and theta_ie

    integer(ik) :: i, j
    integer(ik) :: corner_idx  !< used to select the corner in the reconstructed_state
    integer(ik) :: point_idx  !< used to select the point in the reconstructed_state
    integer(ik), dimension(2, 4) :: neighbor_cell_indices
    real(rk), dimension(2, 2, 4) :: corner_edge_vectors !< ((x,y), (tail,head), (vector1:vector4))
    !< ((x,y), head/tail, vector_id); set of vectors that define the corner

    real(rk), dimension(4, 4) :: reconstructed_corner_state
    !< ((rho, u, v, p), cell_id); the reconstructed state of the corner with respect to each cell

    integer(ik) :: ilo, ihi, jlo, jhi

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

    ! Unlike the midpoints, the corner ref state dimensions are the same in i,j as the # nodes
    ilo = lbound(reference_state, dim=2)
    ihi = ubound(reference_state, dim=2)
    jlo = lbound(reference_state, dim=3)
    jhi = ubound(reference_state, dim=3)

    do concurrent(i=ilo:ihi)
      do concurrent(j=jlo:jhi)

        ! cell ordering is 1) lower left, 2) lower right, 3) upper right, 4) upper left
        neighbor_cell_indices = reshape([[i - 1, j - 1], &  ! lower left
                                         [i, j - 1], &  ! lower right
                                         [i, j], &  ! upper right
                                         [i - 1, j] &  ! upper left
                                         ], shape=[2, 4])

        corner_edge_vectors = self%grid%get_corner_vectors(cell_ij=[i, j], corner='lower-left')

        ! reconstructed_state indexing
        !< ((rho, u ,v, p), point, node/midpoint, i, j);

        ! Cell 1: lower left cell -> corner is in the upper right (N3) of its parent cell
        point_idx = 3
        reconstructed_corner_state(:, 1) = self%reconstructed_state(:, point_idx, corner_idx, i - 1, j - 1)

        ! Cell 2: lower right cell -> corner is in the upper left (N4) of its parent cell
        point_idx = 4
        reconstructed_corner_state(:, 2) = self%reconstructed_state(:, point_idx, corner_idx, i, j - 1)

        ! Cell 3: upper right cell-> corner is in the lower left (N1) of its parent cell
        point_idx = 1
        reconstructed_corner_state(:, 3) = self%reconstructed_state(:, point_idx, corner_idx, i, j)

        ! Cell 4: upper left cell -> corner is in the lower right (N2) of its parent cell
        point_idx = 2
        reconstructed_corner_state(:, 4) = self%reconstructed_state(:, point_idx, corner_idx, i - 1, j)

        mach_cone = new_cone(tau=self%tau, edge_vectors=corner_edge_vectors, &
                             reconstructed_state=reconstructed_corner_state, &
                             reference_state=reference_state(:, i, j), &
                             cell_indices=neighbor_cell_indices)

        ! Set the evolved state at the midpoint
        evolved_state(:, i, j) = [get_density(self, mach_cone), &    ! rho
                                  get_x_velocity(mach_cone), & ! u
                                  get_y_velocity(mach_cone), & ! v
                                  get_pressure(mach_cone) &    ! p
                                  ]
      end do
    end do
  end subroutine evolve_corners

  pure function get_pressure(mach_cone) result(pressure)
    !< Implementation of p(P) within the local evolution operator (Eq 45 in the text)

    class(cone_t), intent(in) :: mach_cone
    real(rk) :: pressure
    integer(ik) :: arc, cell

    pressure = 0.0_rk

    do cell = 1, mach_cone%n_neighbor_cells
      if(mach_cone%n_arcs(cell) > 0) then
        do arc = 1, mach_cone%n_arcs(cell)

          associate(theta_ie=>mach_cone%theta_ie(cell, arc), &
                    theta_ib=>mach_cone%theta_ib(cell, arc), &
                    u_i=>mach_cone%cell_conserved_vars(2, arc, cell), &
                    v_i=>mach_cone%cell_conserved_vars(3, arc, cell), &
                    p_i=>mach_cone%cell_conserved_vars(4, arc, cell), &
                    rho_tilde=>mach_cone%reference_state(1), &
                    a_tilde=>mach_cone%reference_state(4))

            pressure = pressure + (p_i * (theta_ie - theta_ib) - &
                                   rho_tilde * a_tilde * u_i * (sin(theta_ie) - sin(theta_ib)) + &
                                   rho_tilde * a_tilde * v_i * (cos(theta_ie) - cos(theta_ib)))
          end associate
        end do
      end if
    end do

    if(ieee_is_nan(pressure)) error stop "NaNs getting generated by local_evo_operator_t%get_pressure"

    pressure = pressure / (2.0_rk * pi)
    if(near_zero(pressure) .or. pressure < 0.0_rk) error stop "Pressure <= 0 in local_evo_operator_t%get_pressure"
  end function get_pressure

  pure function get_density(self, mach_cone) result(density)
    !< Implementation of rho(P) within the local evolution operator (Eq. 42 in the text)
    class(local_evo_operator_t), intent(in) :: self
    class(cone_t), intent(in) :: mach_cone
    ! real(rk), dimension(:, 0:, 0:), intent(in) :: conserved_vars

    real(rk) :: density  !< rho(P)
    integer(ik) :: arc, cell
    real(rk), dimension(4) :: p_prime_u_bar !< [rho, u, v, p] at P'
    real(rk) :: a
    density = 0.0_rk

    do cell = 1, mach_cone%n_neighbor_cells
      if(mach_cone%n_arcs(cell) > 0) then
        do arc = 1, mach_cone%n_arcs(cell)

          associate(theta_ie=>mach_cone%theta_ie(cell, arc), &
                    theta_ib=>mach_cone%theta_ib(cell, arc), &
                    u_i=>mach_cone%cell_conserved_vars(2, arc, cell), &
                    v_i=>mach_cone%cell_conserved_vars(3, arc, cell), &
                    p_i=>mach_cone%cell_conserved_vars(4, arc, cell), &
                    rho_tilde=>mach_cone%reference_state(1), &
                    a_tilde=>mach_cone%reference_state(4))

            density = density + ((p_i / a_tilde**2) * (theta_ie - theta_ib) &
                                 - (rho_tilde / a_tilde) * u_i * (sin(theta_ie) - sin(theta_ib)) &
                                 + (rho_tilde / a_tilde) * v_i * (cos(theta_ie) - cos(theta_ib)))

            if(ieee_is_nan(density)) then
              error stop "NaNs getting generated by local_evo_operator_t%get_density"
            end if
          end associate
        end do
      end if
    end do

    density = density / (2.0_rk * pi)

    ! Reconstruct at P' to get rho(P') and p(P')
    p_prime_u_bar = self%reconstruction_operator%reconstruct_point(xy=mach_cone%p_prime_xy, &
                                                                   cell_ij=mach_cone%p_prime_ij)

    associate(rho_p_prime=>p_prime_u_bar(1), pressure_p_prime=>p_prime_u_bar(4), &
              a_tilde=>mach_cone%reference_state(4))
      density = density + rho_p_prime - (pressure_p_prime / a_tilde**2)
    end associate

    if(near_zero(density) .or. density < 0.0_rk) then
      error stop "Density <= 0 in local_evo_operator_t%get_density"
    end if
  end function get_density

  pure function get_x_velocity(mach_cone) result(u)
    !< Implementation of u(P) within the local evolution operator (Eq 43 in the text)

    class(cone_t), intent(in) :: mach_cone
    real(rk) :: u !< x velocity at P
    integer(ik) :: arc
    integer(ik) :: cell  !< cell index within the mach cone

    u = 0.0_rk

    do cell = 1, mach_cone%n_neighbor_cells
      if(mach_cone%n_arcs(cell) > 0) then
        do arc = 1, mach_cone%n_arcs(cell)
          ! print*, cell, arc
          associate(theta_ie=>mach_cone%theta_ie(cell, arc), &
                    theta_ib=>mach_cone%theta_ib(cell, arc), &
                    u_i=>mach_cone%cell_conserved_vars(2, arc, cell), &
                    v_i=>mach_cone%cell_conserved_vars(3, arc, cell), &
                    p_i=>mach_cone%cell_conserved_vars(4, arc, cell), &
                    rho_tilde=>mach_cone%reference_state(1), &
                    a_tilde=>mach_cone%reference_state(4))

            u = u + ((-p_i / (rho_tilde * a_tilde)) * (sin(theta_ie) - sin(theta_ib)) + &
                     u_i * (0.5_rk * (theta_ie - theta_ib) + &
                            0.25_rk * (sin(2 * theta_ie) - sin(2 * theta_ib)) &
                            ) - &
                     v_i * 0.25_rk * (cos(2 * theta_ie) - cos(2 * theta_ib)))
          end associate
        end do
      end if
    end do

    if(ieee_is_nan(u)) error stop "NaNs getting generated by local_evo_operator_t%get_x_velocity"
    u = u / pi
  end function get_x_velocity

  pure function get_y_velocity(mach_cone) result(v)
    !< Implementation of v(P) within the local evolution operator (Eq 44 in the text)

    class(cone_t), intent(in) :: mach_cone
    real(rk) :: v !< y velocity at P
    integer(ik) :: arc
    integer(ik) :: cell  !< cell index within the mach cone

    v = 0.0_rk

    do cell = 1, mach_cone%n_neighbor_cells
      if(mach_cone%n_arcs(cell) > 0) then
        do arc = 1, mach_cone%n_arcs(cell)

          associate(theta_ie=>mach_cone%theta_ie(cell, arc), &
                    theta_ib=>mach_cone%theta_ib(cell, arc), &
                    u_i=>mach_cone%cell_conserved_vars(2, arc, cell), &
                    v_i=>mach_cone%cell_conserved_vars(3, arc, cell), &
                    p_i=>mach_cone%cell_conserved_vars(4, arc, cell), &
                    rho_tilde=>mach_cone%reference_state(1), &
                    a_tilde=>mach_cone%reference_state(4))

            v = v + ((p_i / (rho_tilde * a_tilde)) * (cos(theta_ie) - cos(theta_ib)) - &
                     u_i * 0.25_rk * (cos(2 * theta_ie) - cos(2 * theta_ib)) + &
                     v_i * (0.5_rk * (theta_ie - theta_ib) + &
                            0.25_rk * (sin(2 * theta_ie) - sin(2 * theta_ib))))
          end associate
        end do
      end if
    end do

    if(ieee_is_nan(v)) error stop "NaNs getting generated by local_evo_operator_t%get_y_velocity"
    v = v / pi
  end function get_y_velocity

end module mod_local_evo_operator
