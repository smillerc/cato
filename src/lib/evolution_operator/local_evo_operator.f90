module mod_local_evo_operator
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use, intrinsic :: ieee_arithmetic
  use mod_floating_point_utils, only: near_zero, equal
  use math_constants, only: pi
  use mod_abstract_evo_operator, only: abstract_evo_operator_t
  use mod_input, only: input_t
  use mod_mach_cone_geometry, only: mach_cone_geometry_t, new_cone
  use mod_grid, only: grid_t
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_reconstruction_factory, only: reconstruction_factory

  implicit none

  private
  public :: local_evo_operator_t

  ! type(reconstruction_factory_t) :: recon_factory

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

  interface local_evo_operator_t
    module procedure :: constructor
  end interface
contains

  type(local_evo_operator_t) function constructor(input, grid, reconstruction_operator) result(operator)
    !< Constructor for the FVLEG operator
    class(input_t), intent(in) :: input
    class(grid_t), intent(in), target :: grid
    class(abstract_reconstruction_t), intent(in), target :: reconstruction_operator

    ! print*, 'a'
    operator%name = "FVLEG"
    ! print*, 'a'
    ! allocate(operator%grid, source=grid)

    ! print*, 'a'
    ! allocate(operator%reconstruction_operator, source=reconstruction_operator)
    ! print*, 'a'
    call operator%set_tau(input%tau)
    ! print*, 'a'
  end function constructor

  subroutine finalize(self)
    !< Cleanup the operator type
    type(local_evo_operator_t), intent(inout) :: self
    nullify(self%grid)
    ! nullify(self%conserved_vars)
    ! nullify(self%reference_state)
    nullify(self%reconstruction_operator)
  end subroutine finalize

  subroutine evolve_leftright_midpoints(self, reference_state, evolved_state)
    ! subroutine evolve_leftright_midpoints(self, conserved_vars, reconstructed_state, &
    !                                       reference_state, evolved_state)
    ! //TODO: make pure
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

    ! do concurrent(i=ilo:ihi)
    !   do concurrent(j=jlo:jhi)
    do i = ilo, ihi
      do j = jlo, jhi

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
        write(*, '(a, 2(i0,1x), 4(f0.3, 1x))') 'reference_state(:, i, j) ', i, j, reference_state(:, i, j)
        mach_cone = new_cone(tau=self%tau, edge_vectors=midpoint_edge_vectors, &
                             reconstructed_state=reconstructed_midpoint_state, &
                             reference_state=reference_state(:, i, j), &
                             cell_indices=neighbor_cell_indices)
        print *, mach_cone

        ! error stop
        ! Set the evolved state at the midpoint
        evolved_state(:, i, j) = [self%get_density(mach_cone), &    ! rho
                                  self%get_x_velocity(mach_cone), & ! u
                                  self%get_y_velocity(mach_cone), & ! v
                                  self%get_pressure(mach_cone) &    ! p
                                  ]
      end do
    end do
  end subroutine evolve_leftright_midpoints

  subroutine evolve_downup_midpoints(self, reference_state, evolved_state)
    ! subroutine evolve_downup_midpoints(self, conserved_vars, reconstructed_state, &
    !                                    reference_state, evolved_state)
    ! //TODO: make pure
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
    type(mach_cone_geometry_t) :: mach_cone
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

    ! do concurrent(i=ilo:ihi)
    !   do concurrent(j=jlo:jhi)
    do i = ilo, ihi
      do j = jlo, jhi

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
        write(*, '(a, 2(i0,1x), 4(f0.3, 1x))') 'reference_state(:, i, j) ', i, j, reference_state(:, i, j)
        write(*, '(a, 2(i0,1x), 4(f0.3, 1x))') 'reconstructed_midpoint_state(:, 1) ', i, j, reconstructed_midpoint_state(:, 1)
        write(*, '(a, 2(i0,1x), 4(f0.3, 1x))') 'reconstructed_midpoint_state(:, 2) ', i, j, reconstructed_midpoint_state(:, 2)
        write(*, '(a, 2(f0.3,1x))') 'midpoint_edge_vectors 1: ', midpoint_edge_vectors(:, 1, 1)
        write(*, '(a, 2(f0.3,1x))') 'midpoint_edge_vectors 1: ', midpoint_edge_vectors(:, 2, 1)
        write(*, '(a, 2(f0.3,1x))') 'midpoint_edge_vectors 2: ', midpoint_edge_vectors(:, 1, 2)
        write(*, '(a, 2(f0.3,1x))') 'midpoint_edge_vectors 2: ', midpoint_edge_vectors(:, 2, 2)
        mach_cone = new_cone(tau=self%tau, edge_vectors=midpoint_edge_vectors, &
                             reconstructed_state=reconstructed_midpoint_state, &
                             reference_state=reference_state(:, i, j), &
                             cell_indices=neighbor_cell_indices)
        write(*, *) 'any(theta_ib==theta_ie)? ', any(equal(mach_cone%theta_ib, mach_cone%theta_ie))
        write(*, *) mach_cone
        ! error stop
        ! write(*,'(a, 2(i0,1x), 4(f0.3, 1x))') 'mach_cone%reference_state', i,j, mach_cone%reference_state

        ! Set the evolved state at the midpoint
        evolved_state(:, i, j) = [self%get_density(mach_cone), &    ! rho
                                  self%get_x_velocity(mach_cone), & ! u
                                  self%get_y_velocity(mach_cone), & ! v
                                  self%get_pressure(mach_cone) &    ! p
                                  ]
        write(*, '(a, 2(i0,1x), 4(f0.3, 1x))') '  evolved_state(:, i, j) ', i, j, evolved_state(:, i, j)
        print *
      end do
    end do
  end subroutine evolve_downup_midpoints

  subroutine evolve_corners(self, reference_state, evolved_state)
    ! subroutine evolve_corners(self, conserved_vars, reconstructed_state, &
    !                           reference_state, evolved_state)
    ! //TODO: make pure
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
    type(mach_cone_geometry_t) :: mach_cone
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

    do i = ilo, ihi
      do j = jlo, jhi

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
        evolved_state(:, i, j) = [self%get_density(mach_cone), &    ! rho
                                  self%get_x_velocity(mach_cone), & ! u
                                  self%get_y_velocity(mach_cone), & ! v
                                  self%get_pressure(mach_cone) &    ! p
                                  ]
      end do
    end do
  end subroutine evolve_corners

  function get_pressure(mach_cone) result(pressure)
    !< Implementation of p(P) within the local evolution operator (Eq 45 in the text)

    class(mach_cone_geometry_t), intent(in) :: mach_cone
    real(rk) :: pressure
    integer(ik) :: arc, cell

    pressure = 0.0_rk

    do cell = 1, mach_cone%n_neighbor_cells
      do arc = 1, mach_cone%n_arcs(cell)

        associate(theta_ie=>mach_cone%theta_ie(cell, arc), &
                  theta_ib=>mach_cone%theta_ib(cell, arc), &
                  u_i=>mach_cone%cell_conserved_vars(2, arc, cell), &
                  v_i=>mach_cone%cell_conserved_vars(3, arc, cell), &
                  p_i=>mach_cone%cell_conserved_vars(4, arc, cell), &
                  rho_tilde=>mach_cone%reference_state(1), &
                  a_tilde=>mach_cone%reference_state(4))

          pressure = pressure + p_i * (theta_ie - theta_ib) - &
                     rho_tilde * a_tilde * u_i * (sin(theta_ie) - sin(theta_ib)) + &
                     rho_tilde * a_tilde * v_i * (cos(theta_ie) - cos(theta_ib))
          ! write(*, '(8(f0.3,1x))') pressure, p_i, theta_ie, theta_ib, rho_tilde, a_tilde, u_i, v_i
        end associate
      end do
    end do

    if(ieee_is_nan(pressure)) error stop "NaNs getting generated by local_evo_operator_t%get_pressure"
    if(near_zero(pressure) .or. pressure < 0.0_rk) error stop "Pressure <= 0 in local_evo_operator_t%get_pressure"
    pressure = pressure / (2.0_rk * pi)

  end function get_pressure

  function get_density(self, mach_cone) result(density)
    ! //TODO: make pure
    !< Implementation of rho(P) within the local evolution operator (Eq. 42 in the text)
    class(local_evo_operator_t), intent(in) :: self
    class(mach_cone_geometry_t), intent(in) :: mach_cone
    ! real(rk), dimension(:, 0:, 0:), intent(in) :: conserved_vars

    real(rk) :: density  !< rho(P)
    integer(ik) :: arc, cell
    real(rk), dimension(4) :: p_prime_u_bar !< [rho, u, v, p] at P'

    density = 0.0_rk

    do cell = 1, mach_cone%n_neighbor_cells
      do arc = 1, mach_cone%n_arcs(cell)

        associate(theta_ie=>mach_cone%theta_ie(cell, arc), &
                  theta_ib=>mach_cone%theta_ib(cell, arc), &
                  u_i=>mach_cone%cell_conserved_vars(2, arc, cell), &
                  v_i=>mach_cone%cell_conserved_vars(3, arc, cell), &
                  p_i=>mach_cone%cell_conserved_vars(4, arc, cell), &
                  rho_tilde=>mach_cone%reference_state(1), &
                  a_tilde=>mach_cone%reference_state(4))

          density = density + (p_i / a_tilde**2) * (theta_ie - theta_ib) - &
                    (rho_tilde / a_tilde) * u_i * (sin(theta_ie) - sin(theta_ib)) + &
                    (rho_tilde / a_tilde) * v_i * (cos(theta_ie) - cos(theta_ib))
          if(ieee_is_nan(density)) then
            write(*, '(8(f0.3,1x))') density, theta_ie, theta_ib, p_i, a_tilde, rho_tilde, u_i, v_i
            error stop "NaNs getting generated by local_evo_operator_t%get_density"
          end if
        end associate
      end do
    end do

    if(near_zero(density) .or. density < 0.0_rk) error stop "Density <= 0 in local_evo_operator_t%get_density"
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
    integer(ik) :: arc
    integer(ik) :: cell  !< cell index within the mach cone

    u = 0.0_rk

    do cell = 1, mach_cone%n_neighbor_cells
      do arc = 1, mach_cone%n_arcs(cell)

        associate(theta_ie=>mach_cone%theta_ie(cell, arc), &
                  theta_ib=>mach_cone%theta_ib(cell, arc), &
                  u_i=>mach_cone%cell_conserved_vars(2, arc, cell), &
                  v_i=>mach_cone%cell_conserved_vars(3, arc, cell), &
                  p_i=>mach_cone%cell_conserved_vars(4, arc, cell), &
                  rho_tilde=>mach_cone%reference_state(1), &
                  a_tilde=>mach_cone%reference_state(4))

          u = u - (p_i / (rho_tilde * a_tilde)) * (sin(theta_ie) - sin(theta_ib)) + &
              u_i * (0.5_rk * (theta_ie - theta_ib) + 0.25_rk * (sin(2 * theta_ie) - sin(2 * theta_ib))) - &
              v_i * 0.25_rk * (cos(2 * theta_ie) - cos(2 * theta_ib))

        end associate
      end do
    end do

    if(ieee_is_nan(u)) error stop "NaNs getting generated by local_evo_operator_t%get_x_velocity"
    u = u / pi
  end function get_x_velocity

  pure function get_y_velocity(mach_cone) result(v)
    !< Implementation of v(P) within the local evolution operator (Eq 44 in the text)

    class(mach_cone_geometry_t), intent(in) :: mach_cone
    real(rk) :: v !< y velocity at P
    integer(ik) :: arc
    integer(ik) :: cell  !< cell index within the mach cone

    v = 0.0_rk

    do cell = 1, mach_cone%n_neighbor_cells
      do arc = 1, mach_cone%n_arcs(cell)

        associate(theta_ie=>mach_cone%theta_ie(cell, arc), &
                  theta_ib=>mach_cone%theta_ib(cell, arc), &
                  u_i=>mach_cone%cell_conserved_vars(2, arc, cell), &
                  v_i=>mach_cone%cell_conserved_vars(3, arc, cell), &
                  p_i=>mach_cone%cell_conserved_vars(4, arc, cell), &
                  rho_tilde=>mach_cone%reference_state(1), &
                  a_tilde=>mach_cone%reference_state(4))

          v = v + (p_i / (rho_tilde * a_tilde)) * (cos(theta_ie) - cos(theta_ib)) + &
              u_i * 0.25_rk * (cos(2 * theta_ie) - cos(2 * theta_ib)) + &
              v_i * (0.5_rk * (theta_ie - theta_ib) + 0.25_rk * (sin(2 * theta_ie) - sin(2 * theta_ib)))
        end associate
      end do
    end do

    if(ieee_is_nan(v)) error stop "NaNs getting generated by local_evo_operator_t%get_y_velocity"
    v = v / pi
  end function get_y_velocity

end module mod_local_evo_operator
