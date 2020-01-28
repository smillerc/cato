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
  public :: local_evo_operator_t

  type, extends(abstract_evo_operator_t) :: local_evo_operator_t
    !< Local Evolution Operator (E0)
  contains
    procedure, public :: initialize
    procedure, public :: evolve_leftright_midpoints
    procedure, public :: evolve_downup_midpoints
    procedure, public :: evolve_corners
    procedure, public :: e0_operator
    procedure, public :: copy
    final :: finalize
  end type local_evo_operator_t

contains

  subroutine initialize(self, input, grid_target, recon_operator_target)
    !< Constructor for the FVLEG operator
    class(local_evo_operator_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    class(grid_t), intent(in), target :: grid_target
    class(abstract_reconstruction_t), intent(in), target :: recon_operator_target

    self%name = "FVLEG"
    call self%set_tau(input%tau)
    self%grid => grid_target
    self%reconstruction_operator => recon_operator_target

  end subroutine

  subroutine finalize(self)
    !< Cleanup the operator type
    type(local_evo_operator_t), intent(inout) :: self

    call debug_print('Running local_evo_operator_t%finalize()', __FILE__, __LINE__)

    if(allocated(self%name)) deallocate(self%name)
    if(associated(self%grid)) nullify(self%grid)
    if(associated(self%reconstructed_state)) nullify(self%reconstructed_state)
    if(associated(self%reconstruction_operator)) nullify(self%reconstruction_operator)
  end subroutine finalize

  subroutine copy(out_evo, in_evo)
    class(abstract_evo_operator_t), intent(in) :: in_evo
    class(local_evo_operator_t), intent(inout) :: out_evo

    call debug_print('Running local_evo_operator_t%copy()', __FILE__, __LINE__)

    if(associated(out_evo%grid)) nullify(out_evo%grid)
    out_evo%grid => in_evo%grid

    if(associated(out_evo%reconstructed_state)) nullify(out_evo%reconstructed_state)
    out_evo%reconstructed_state => in_evo%reconstructed_state

    if(associated(out_evo%reconstruction_operator)) nullify(out_evo%reconstruction_operator)
    out_evo%reconstruction_operator => in_evo%reconstruction_operator

    if(allocated(out_evo%name)) deallocate(out_evo%name)
    allocate(out_evo%name, source=in_evo%name)

  end subroutine

  subroutine evolve_leftright_midpoints(self, reference_state, evolved_state, lbounds)
    !< Create Mach cones and evolve the state at all of the left/right midpoints in the domain. Left/right midpoints
    !< are the midpoints defined by vectors pointing to the left and right.

    class(local_evo_operator_t), intent(in) :: self

    integer(ik), dimension(3), intent(in) :: lbounds

    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(in) :: reference_state
    !< ((rho, u, v, p), i, j); Reference state (tilde) at each midpoint on the left/right edges

    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(out) :: evolved_state
    !< ((rho, u, v, p), i, j); Reconstructed U at each midpoint on the left/right edges

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

    ! print*, 'left/right -> ', ilo, ihi, jlo, jhi
    do j = jlo, jhi
      do i = ilo, ihi
        ! do concurrent(j=jlo:jhi) local(neighbor_cell_indices, midpoint_edge_vectors, reconstructed_midpoint_state)
        !   do concurrent(i=ilo:ihi) local(neighbor_cell_indices, midpoint_edge_vectors, reconstructed_midpoint_state)

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
        evolved_state(:, i, j) = self%e0_operator(mach_cone)

      end do
    end do
  end subroutine evolve_leftright_midpoints

  subroutine evolve_downup_midpoints(self, reference_state, evolved_state, lbounds)
    !< Create Mach cones and evolve the state at all of the down/up midpoints in the domain. Left/right midpoints
    !< are the midpoints defined by vectors pointing to the down and up.

    class(local_evo_operator_t), intent(in) :: self

    integer(ik), dimension(3), intent(in) :: lbounds

    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(in) :: reference_state
    ! < ((rho, u, v, p), i, j); Reference state (tilde) at each midpoint on the left/right edges

    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(out) :: evolved_state
    !< ((rho, u, v, p), i, j); Evolved U at each midpoint on the left/right edges

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

    ! print*, 'down/up -> ', ilo, ihi, jlo, jhi
    do j = jlo, jhi
      do i = ilo, ihi
        ! do concurrent(j=jlo:jhi)
        !   do concurrent(i=ilo:ihi)

        ! cell ordering is 1) left, 2) right
        neighbor_cell_indices = reshape([[i - 1, j], &  ! cell left
                                         [i, j] &  ! cell right
                                         ], shape=[2, 2])

        midpoint_edge_vectors = self%grid%get_midpoint_vectors(cell_ij=[i, j], edge='left')
        ! debug_write(*,*) i, j
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
        evolved_state(:, i, j) = self%e0_operator(mach_cone)
      end do
    end do
  end subroutine evolve_downup_midpoints

  subroutine evolve_corners(self, reference_state, evolved_state, lbounds)
    !< Create Mach cones and evolve the state at all of the down/up midpoints in the domain. Left/right midpoints
    !< are the midpoints defined by vectors pointing to the down and up.

    class(local_evo_operator_t), intent(in) :: self

    integer(ik), dimension(3), intent(in) :: lbounds

    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(in) :: reference_state
    !< ((rho, u, v, p), i, j); Reference state (tilde) at each corner

    real(rk), dimension(lbounds(1):, lbounds(2):, lbounds(3):), intent(out) :: evolved_state
    !< ((rho, u, v, p), i, j); Reconstructed U at each corner

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

    ! print*, 'corner -> ', ilo, ihi, jlo, jhi
    do j = jlo, jhi
      do i = ilo, ihi

        ! do concurrent(j=jlo:jhi) local(neighbor_cell_indices, corner_edge_vectors)
        !   do concurrent(i=ilo:ihi) local(neighbor_cell_indices, corner_edge_vectors)

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
        evolved_state(:, i, j) = self%e0_operator(mach_cone)
      end do
    end do
  end subroutine evolve_corners

  function e0_operator(self, mach_cone) result(primitive_variables)

    class(local_evo_operator_t), intent(in) :: self
    class(cone_t), intent(in) :: mach_cone
    real(rk), dimension(4) :: primitive_variables
    real(rk), dimension(4) :: p_prime_prim_vars !< [rho, u, v, p] at P'

    real(rk) :: cone_density_term, pressure, density, u, v

    primitive_variables = 0.0_rk

    p_prime_prim_vars = self%reconstruction_operator%reconstruct_point(xy=mach_cone%p_prime_xy, &
                                                                       cell_ij=mach_cone%p_prime_ij)

    associate(theta_ie=>mach_cone%theta_ie, &
              theta_ib=>mach_cone%theta_ib, &
              rho_p_prime=>p_prime_prim_vars(1), &
              pressure_p_prime=>p_prime_prim_vars(4), &
              u_i=>mach_cone%arc_primitive_vars(2, :, :), &
              v_i=>mach_cone%arc_primitive_vars(3, :, :), &
              p_i=>mach_cone%arc_primitive_vars(4, :, :), &
              rho_tilde=>mach_cone%reference_state(1), &
              a_tilde=>mach_cone%reference_state(4))

      pressure = sum(p_i * (theta_ie - theta_ib) &
                     - rho_tilde * a_tilde * u_i * (sin(theta_ie) - sin(theta_ib)) &
                     + rho_tilde * a_tilde * v_i * (cos(theta_ie) - cos(theta_ib))) / (2.0_rk * pi)

      cone_density_term = sum((p_i / a_tilde**2) * (theta_ie - theta_ib) &
                              - (rho_tilde / a_tilde) * u_i * (sin(theta_ie) - sin(theta_ib)) &
                              + (rho_tilde / a_tilde) * v_i * (cos(theta_ie) - cos(theta_ib))) / (2.0_rk * pi)

      density = rho_p_prime - (pressure_p_prime / a_tilde**2) + cone_density_term
      ! if(density < 0.0_rk) then
      !   print *, '===================================='
      !   print *, ' density < 0.0, trying pressure fix @ (x,y) = ', mach_cone%p_xy
      !   print *, '===================================='
      !   density = rho_p_prime - (pressure / a_tilde**2) + cone_density_term
      ! end if

      u = sum((-p_i / (rho_tilde * a_tilde)) * (sin(theta_ie) - sin(theta_ib)) &
              + u_i * ( &
              ((theta_ie - theta_ib) / 2.0_rk) &
              + ((sin(2.0_rk * theta_ie) - sin(2.0_rk * theta_ib)) / 4.0_rk) &
              ) &
              - v_i * ((cos(2.0_rk * theta_ie) - cos(2.0_rk * theta_ib)) / 4.0_rk) &
              ) / pi

      v = sum((p_i / (rho_tilde * a_tilde)) * (cos(theta_ie) - cos(theta_ib)) &
              - u_i * ((cos(2.0_rk * theta_ie) - cos(2.0_rk * theta_ib)) / 4.0_rk) &
              + v_i * ( &
              ((theta_ie - theta_ib) / 2.0_rk) &
              - ((sin(2.0_rk * theta_ie) - sin(2.0_rk * theta_ib)) / 4.0_rk) &
              ) &
              ) / pi

      if(density < 0.0_rk) then
        print *, mach_cone
        error stop "Density < 0 in e0_operator"
      end if

      if(pressure < 0.0_rk) then
        print *, mach_cone
        error stop "Pressure < 0 in e0_operator"
      end if
    end associate

    primitive_variables = [density, u, v, pressure]

  end function

end module mod_local_evo_operator
