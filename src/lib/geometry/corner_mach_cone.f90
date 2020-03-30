module mod_corner_mach_cone
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use math_constants, only: pi, rad2deg
  use mod_globals, only: TINY_VEL, grid_is_orthogonal
  use mod_floating_point_utils, only: near_zero, equal
  use mod_eos, only: eos
  use mod_vector, only: vector_t
  use mod_geometry, only: super_circle, get_arc_segments_new
  use mod_mach_cone_utilties, only: get_cone_extents, determine_if_p_prime_is_in_cell

  implicit none

  private
  public :: corner_mach_cone_t, new_corner_cone

  integer(ik), parameter :: N_CELLS = 4

  type :: corner_mach_cone_t
    real(rk), dimension(4, N_CELLS) :: recon_state = 0.0_rk
    !< ((rho, u, v, p), (cell 1:N_CELLS)); Reconstructed value of U at P for each neighbor cell

    logical, dimension(N_CELLS) :: cell_is_supersonic = .false.
    !< (cell 1:N_CELLS); Is each neighbor cell supersonic or not? This is needed for transonic mach cones in
    !< order to apply an entropy fix for transonic rarefaction regions.

    ! real(rk), dimension(2, N_CELLS) :: theta_ie = 0.0_rk
    ! !< ((arc 1:2), (cell 1:N_CELLS)); arc end angle

    ! real(rk), dimension(2, N_CELLS) :: theta_ib = 0.0_rk
    ! !< ((arc 1:2), (cell 1:N_CELLS)); arc begin angle

    real(rk), dimension(2, N_CELLS) :: dtheta = 0.0_rk
    !< ((arc 1:2), (cell 1:N_CELLS)); theta_ie - theta_ib [radians]

    real(rk), dimension(2, N_CELLS) :: sin_dtheta = 0.0_rk
    !< ((arc 1:2), (cell 1:N_CELLS)); sin(theta_ie) - sin(theta_ib)

    real(rk), dimension(2, N_CELLS) :: cos_dtheta = 0.0_rk
    !< ((arc 1:2), (cell 1:N_CELLS)); cos(theta_ie) - cos(theta_ib)

    real(rk), dimension(2, N_CELLS) :: sin_d2theta = 0.0_rk
    !< ((arc 1:2), (cell 1:N_CELLS)); sin(2*theta_ie) - sine(2*theta_ib)

    real(rk), dimension(2, N_CELLS) :: cos_d2theta = 0.0_rk
    !< ((arc 1:2), (cell 1:N_CELLS)); cos(2*theta_ie) - cos(2*theta_ib)

    integer(ik), dimension(2, N_CELLS) :: p_prime_ij = 0
    !< ((i,j), (cell 1:N_CELLS)) i,j location of P' or the apex of the Mach cone (global, not relative to P0)

    real(rk), dimension(2, 2, N_CELLS) :: edge_vectors = 0.0_rk
    !< ((x,y), (tail,head), (vector 1:N_CELLS)); set of vectors that define the cell edges

    integer(ik), dimension(N_CELLS) :: n_arcs_per_cell = 0
    !< (cell 1:N_CELLS); Number of arcs that each cell contributes (0, 1, or 2)

    logical, dimension(N_CELLS) :: p_prime_in_cell = .false.
    !< (cell 1:N_CELLS); is P' inside the control volume?

    real(rk), dimension(4, 2, N_CELLS) :: arc_primitive_vars = 0.0_rk
    !< ((rho,u,v,p), (arc 1:2), (cell 1:N_CELLS))

    real(rk), dimension(N_CELLS) :: sound_speed = 0.0_rk
    !< (cell 1:N_CELLS); speed of sound for each neighbor cell

    real(rk), dimension(2) :: p_prime_xy = 0.0_rk !< (x,y); Location of P'
    real(rk), dimension(2) :: p_xy = 0.0_rk       !< (x,y); Location of P
    integer(ik) :: n_neighbor_cells = N_CELLS           !< Number of neighbor cells
    real(rk) :: radius = -1.0_rk                  !< Radius of the cone
    real(rk) :: tau = 1.0e-10_rk                  !< Time evolution increment
    real(rk) :: reference_density = 0.0_rk        !< Reference density (e.g. neighbor averaged)
    real(rk) :: reference_u = 0.0_rk              !< Reference x velocity (e.g. neighbor averaged)
    real(rk) :: reference_v = 0.0_rk              !< Reference y velocity (e.g. neighbor averaged)
    real(rk) :: reference_sound_speed = 0.0_rk    !< Reference sound speed (e.g. neighbor averaged)
    real(rk) :: reference_mach_number = 0.0_rk    !< Reference Mach number (e.g. neighbor averaged)
    logical :: cone_is_transonic = .false.        !< Flag to enable special treatment for transonic cones
    logical :: cone_is_centered = .false.
    !< Flag signifying if P and P' are collocated -> allows for simplified angles/trig
    !< so a bunch of functions can be skipped for speed (no need to find line/circle intersections)

    character(len=32) :: cone_location = 'corner'       !< Corner or midpoint cone?
  contains
    procedure, private :: get_reference_state
    procedure, private :: determine_p_prime_cell
    procedure, private :: sanity_checks
    procedure, private :: find_arc_angles
    procedure, private :: precompute_trig_angles
    procedure, private :: get_transonic_cone_extents
    procedure, pass :: write => write_cone
    generic, public :: write(formatted) => write
  end type corner_mach_cone_t

contains

  function new_corner_cone(tau, edge_vectors, reconstructed_state, cell_indices, cone_location) result(cone)
    !< Constructor for the corner cone type

    type(corner_mach_cone_t) :: cone
    real(rk), intent(in) :: tau  !< time increment, tau -> 0 (very small number)
    character(len=*), intent(in) :: cone_location !< corner, down/up midpoint, left/right midpoint

    real(rk), dimension(2, 2, N_CELLS), intent(in) :: edge_vectors
    !< ((x,y), (tail,head), (vector 1:N_CELLS)); set of vectors that define the cell edges

    integer(ik), dimension(2, N_CELLS), intent(in) :: cell_indices
    !< ((i,j), cell 1:N_CELLS); set of indices for the neighboring cells -> needed to find P' i,j index

    real(rk), dimension(4, N_CELLS), intent(in) :: reconstructed_state
    !< ((rho,u,v,p), cell); reconstructed state for point P.

    ! Locals
    integer(ik) :: l, i, arc

    cone%tau = tau
    cone%p_xy = [edge_vectors(1, 1, 1), edge_vectors(2, 1, 1)]
    cone%cone_location = trim(cone_location)
    cone%edge_vectors = edge_vectors

    if(any(reconstructed_state(1, :) < 0.0_rk)) then
      write(std_err, '(a, 4(es10.3,1x))') 'Reconstructed density states (cell 1:4)', reconstructed_state(1, :)
      write(std_err, '(a, 4("[",i0,", ", i0,"] "))') 'Midpoint [i, j] cell indices: ', cell_indices
      error stop "Error in corner_mach_cone_t initialization, density in the reconstructed state is < 0"

    end if

    if(any(reconstructed_state(4, :) < 0.0_rk)) then
      write(std_err, '(a, 4(es10.3,1x))') 'Reconstructed pressure states (cell 1:4): ', reconstructed_state(4, :)
      write(std_err, '(a, 4("[",i0,", ", i0,"] "))') 'Midpoint [i, j] cell indices: ', cell_indices
      error stop "Error in corner_mach_cone_t initialization, pressure in the reconstructed state is < 0"
    end if

    do i = 1, N_CELLS
      do l = 1, 4
        cone%recon_state(l, i) = reconstructed_state(l, i)
      end do
      ! Tiny velocity fix to reduce noise
      if(abs(cone%recon_state(2, i)) < TINY_VEL) cone%recon_state(2, i) = 0.0_rk
      if(abs(cone%recon_state(3, i)) < TINY_VEL) cone%recon_state(3, i) = 0.0_rk
    end do

    call cone%get_reference_state(reconstructed_state)

    if(cone%cone_is_transonic) then
      call cone%get_transonic_cone_extents(origin=cone%p_prime_xy, &
                                           radius=cone%radius)
      cone%reference_sound_speed = cone%radius / cone%tau
    else
      call get_cone_extents(tau=cone%tau, xy=cone%p_xy, &
                            vel=[cone%reference_u, cone%reference_v], &
                            sound_speed=cone%reference_sound_speed, &
                            origin=cone%p_prime_xy, &
                            radius=cone%radius)
    end if

    call cone%precompute_trig_angles(cell_indices, edge_vectors)

    ! Assign the primitive state variables to each arc contained in each cell
    do i = 1, N_CELLS
      do arc = 1, cone%n_arcs_per_cell(i)
        cone%arc_primitive_vars(:, arc, i) = cone%recon_state(:, i)
      end do
    end do

    call cone%determine_p_prime_cell()
    call cone%sanity_checks()

  end function new_corner_cone

  subroutine get_reference_state(self, reconstructed_state)
    !< Find the reference state of the cone only based on the cells that are touched by the cone.
    !< Sometimes when a cone is near a large density jump, a skewed reference state can cause
    !< negative densities and pressures. This aims to alleviate that problem...

    class(corner_mach_cone_t), intent(inout) :: self
    real(rk), dimension(4, N_CELLS), intent(in) :: reconstructed_state !< [rho, u, v, a]

    real(rk) :: ave_p
    real(rk), dimension(N_CELLS) :: mach_number, sound_speed
    integer(ik) :: i, n_supersonic_cells
    real(rk) :: gamma

    gamma = eos%get_gamma()

    n_supersonic_cells = 0
    do i = 1, N_CELLS
      self%sound_speed(i) = sqrt(gamma * self%recon_state(4, i) / self%recon_state(1, i))
    end do

    do i = 1, N_CELLS
      mach_number(i) = sqrt(self%recon_state(2, i)**2 + &
                            self%recon_state(3, i)**2) / self%sound_speed(i)
      if(mach_number(i) > 1.0_rk) then
        self%cell_is_supersonic(i) = .true.
        n_supersonic_cells = n_supersonic_cells + 1
      end if
    end do

    ! n_supersonic_cells = count(self%cell_is_supersonic)
    if(n_supersonic_cells == N_CELLS .or. n_supersonic_cells == 0) then
      self%cone_is_transonic = .false.
    else
      self%cone_is_transonic = .true.
    end if

    self%reference_density = sum(self%recon_state(1, :)) / real(N_CELLS, rk)
    self%reference_u = sum(self%recon_state(2, :)) / real(N_CELLS, rk)
    self%reference_v = sum(self%recon_state(3, :)) / real(N_CELLS, rk)

    if(abs(self%reference_u) < TINY_VEL .and. abs(self%reference_v) < TINY_VEL) then
      self%reference_u = 0.0_rk
      self%reference_v = 0.0_rk
      self%cone_is_centered = .true.
    else
      self%cone_is_centered = .false.
    end if

    ave_p = sum(self%recon_state(4, :)) / real(N_CELLS, rk)
    self%reference_sound_speed = sqrt(gamma * ave_p / self%reference_density)
    self%reference_mach_number = sqrt(self%reference_u**2 + self%reference_v**2) / &
                                 self%reference_sound_speed
  end subroutine get_reference_state

  subroutine find_arc_angles(self, cell_indices, edge_vectors, theta_ib, theta_ie, n_arcs_per_cell)
    !< Find the starting and ending angle of each arc that is held within each cell

    class(corner_mach_cone_t), intent(inout) :: self

    real(rk), dimension(2, 2, N_CELLS), intent(in) :: edge_vectors
    integer(ik), dimension(2, N_CELLS), intent(in) :: cell_indices
    !< ((i,j), cell 1:N_CELLS); set of indices for the neighboring cells -> needed to find P' i,j index

    real(rk), dimension(2, N_CELLS), intent(out) :: theta_ib
    !< ((arc1, arc2), (cell 1:N_CELLS)); starting angle [rad] for the arc contained in each cell

    real(rk), dimension(2, N_CELLS), intent(out) :: theta_ie
    !< ((arc1, arc2), (cell 1:N_CELLS)); ending angle [rad] for the arc contained in each cell

    integer(ik), dimension(N_CELLS), intent(out) :: n_arcs_per_cell !< # of arcs in a given cell

    real(rk), dimension(2, 2) :: theta_ib_ie
    !< ((start,end), (arc1, arc2)); start/end angle for a given cell

    real(rk), dimension(2) :: vector_1_head !< (x,y) location of vector 1's head
    real(rk), dimension(2) :: vector_2_head !< (x,y) location of vector 2's head
    real(rk), dimension(2) :: vector_3_head !< (x,y) location of vector 3's head
    real(rk), dimension(2) :: vector_4_head !< (x,y) location of vector 4's head
    real(rk), dimension(2) :: origin        !< (x,y) location of the origin (common for all vectors)
    real(rk), dimension(2) :: second_vector, first_vector
    type(vector_t) :: p_prime_vector
    integer(ik) :: i
    logical :: p_prime_in_cell
    integer(ik), dimension(2) :: cell_ij
    integer(ik) :: arc  !< index for looping through the # arcs in each neighbor cell
    integer(ik) :: n_arcs

    cell_ij = 0
    ! corner cone (for quadrilateral cells)
    !       cell 4                   cell 3
    !      (i-1,j)                   (i,j)
    !  N4----M3----N3    P3   N4----M3----N3
    !  |            |    |    |            |
    !  M4    C4    M2    |    M4    C3    M2
    !  |            |    |    |            |
    !  N1----M1----N2    |    N1----M1----N2
    !                    |
    !  P4----------------O-----------------P2
    !                    |
    !  N4----M3----N3    |    N4----M3----N3
    !  |            |    |    |            |
    !  M4    C1    M2    |    M4    C2    M2
    !  |            |    |    |            |
    !  N1----M1----N2    P1   N1----M1----N2
    !      cell 1                  cell 2
    !     (i-1,j-1)               (i,j-1)

    ! For corners, the edge vectors go down (0-P1), right(O-P2), up(O-P3), right (O-P4)
    ! O is the corner point shared by the 4 cells

    ! the order of the vectors matters, mainly for the cross product to determine
    ! if P' is in the neighboring cell or not
    ! edge_vector_ordering = [4, 1],[1, 2],[2, 3],[3, 4]

    vector_1_head = [edge_vectors(1, 2, 1), edge_vectors(2, 2, 1)] ! (x,y)
    vector_2_head = [edge_vectors(1, 2, 2), edge_vectors(2, 2, 2)] ! (x,y)
    vector_3_head = [edge_vectors(1, 2, 3), edge_vectors(2, 2, 3)] ! (x,y)
    vector_4_head = [edge_vectors(1, 2, 4), edge_vectors(2, 2, 4)] ! (x,y)
    origin = [edge_vectors(1, 1, 1), edge_vectors(2, 1, 1)] ! (x,y)

    p_prime_vector = vector_t(x=[edge_vectors(1, 1, 1), self%p_prime_xy(1)], &
                              y=[edge_vectors(2, 1, 1), self%p_prime_xy(2)])

    do i = 1, N_CELLS
      select case(i)
      case(1) ! Cell 1
        first_vector = vector_4_head
        second_vector = vector_1_head
      case(2) ! Cell 2
        first_vector = vector_1_head
        second_vector = vector_2_head
      case(3) ! Cell 3
        first_vector = vector_2_head
        second_vector = vector_3_head
      case(4) ! Cell 4
        first_vector = vector_3_head
        second_vector = vector_4_head
      end select

      p_prime_in_cell = determine_if_p_prime_is_in_cell(origin=origin, &
                                                        vec_1_head=first_vector, &
                                                        vec_2_head=second_vector, &
                                                        p_prime_vector=p_prime_vector)
      self%p_prime_in_cell(i) = p_prime_in_cell
      if(p_prime_in_cell) self%p_prime_ij(:, i) = cell_indices(:, i)
      call get_arc_segments_new(origin=origin, &
                                vec_1_head=first_vector, &
                                vec_2_head=second_vector, &
                                origin_in_cell=p_prime_in_cell, &
                                circle_xy=self%p_prime_xy, circle_radius=self%radius, &
                                arc_segments=theta_ib_ie, n_arcs=n_arcs)
      n_arcs_per_cell(i) = n_arcs
      theta_ib(:, i) = theta_ib_ie(1, :)
      theta_ie(:, i) = theta_ib_ie(2, :)
    end do

  end subroutine find_arc_angles

  subroutine precompute_trig_angles(self, cell_indices, edge_vectors)
    !< Precompute some of the trig values, like sin(theta_ie), etc. to save compute time. If
    !< the grid is orthogonal and the references state is at rest, then the math
    !< gets even easier

    class(corner_mach_cone_t), intent(inout) :: self

    real(rk), dimension(2, N_CELLS) :: theta_ib
    !< ((arc1, arc2), (cell 1:N_CELLS)); starting angle [rad] for the arc contained in each cell

    real(rk), dimension(2, N_CELLS) :: theta_ie
    !< ((arc1, arc2), (cell 1:N_CELLS)); ending angle [rad] for the arc contained in each cell

    real(rk), dimension(2, 2, N_CELLS), intent(in) :: edge_vectors
    integer(ik), dimension(2, N_CELLS), intent(in) :: cell_indices
    !< ((i,j), cell 1:N_CELLS); set of indices for the neighboring cells -> needed to find P' i,j index

    real(rk), dimension(2, N_CELLS) :: sin_theta_ib, cos_theta_ib
    real(rk), dimension(2, N_CELLS) :: sin_theta_ie, cos_theta_ie
    integer(ik) :: i, arc, max_n_arcs
    integer(ik), dimension(N_CELLS) :: n_arcs_per_cell

    ! If we know the cone is at the origin and the grid is orthogonal, then
    ! the angles are easily known ahead of time, so we can skip the additional
    ! math needed to find the intersection angles
    if(self%cone_is_centered .and. grid_is_orthogonal) then

      self%dtheta(1, :) = pi / 2.0_rk ! all of the 1st arcs are a quarter of the circle
      self%n_arcs_per_cell = 1
      self%p_prime_in_cell = [.true., .false., .false., .false.]
      ! Cell 1 (lower left):  180 -> 270
      self%sin_dtheta(1, 1) = -1.0_rk
      self%cos_dtheta(1, 1) = 1.0_rk
      self%sin_d2theta(1, 1) = 0.0_rk
      self%cos_d2theta(1, 1) = -2.0_rk

      ! Cell 2 (lower right): 270 -> 360
      self%sin_dtheta(1, 2) = 1.0_rk
      self%cos_dtheta(1, 2) = 1.0_rk
      self%sin_d2theta(1, 2) = 0.0_rk
      self%cos_d2theta(1, 2) = 2.0_rk

      ! Cell 3 (upper right): 0 -> 90
      self%sin_dtheta(1, 3) = 1.0_rk
      self%cos_dtheta(1, 3) = -1.0_rk
      self%sin_d2theta(1, 3) = 0.0_rk
      self%cos_d2theta(1, 3) = -2.0_rk

      ! Cell 4 (upper left):  90 -> 180
      self%sin_dtheta(1, 4) = -1.0_rk
      self%cos_dtheta(1, 4) = -1.0_rk
      self%sin_d2theta(1, 4) = 0.0_rk
      self%cos_d2theta(1, 4) = 2.0_rk

    else
      call self%find_arc_angles(cell_indices, edge_vectors, theta_ib, theta_ie, n_arcs_per_cell)

      ! Precompute for a bit of speed
      self%n_arcs_per_cell = n_arcs_per_cell
      max_n_arcs = maxval(n_arcs_per_cell)

      if(max_n_arcs == 1) then ! most cells will only have 1 arc in each
        do i = 1, N_CELLS

          self%dtheta(1, i) = abs(theta_ie(1, i) - theta_ib(1, i))
          sin_theta_ib(1, i) = sin(theta_ib(1, i))
          cos_theta_ib(1, i) = cos(theta_ib(1, i))
          sin_theta_ie(1, i) = sin(theta_ie(1, i))
          cos_theta_ie(1, i) = cos(theta_ie(1, i))
        end do
      else ! occasionally, some will have 2 arcs in a cell, so we loop through them all
        do i = 1, N_CELLS
          do arc = 1, max_n_arcs
            self%dtheta(arc, i) = abs(theta_ie(arc, i) - theta_ib(arc, i))
            sin_theta_ib(arc, i) = sin(theta_ib(arc, i))
            cos_theta_ib(arc, i) = cos(theta_ib(arc, i))
            sin_theta_ie(arc, i) = sin(theta_ie(arc, i))
            cos_theta_ie(arc, i) = cos(theta_ie(arc, i))
          end do
        end do
      end if

      self%sin_dtheta = sin_theta_ie - sin_theta_ib
      self%cos_dtheta = cos_theta_ie - cos_theta_ib
      self%sin_d2theta = 2.0_rk * sin_theta_ie * cos_theta_ie - 2.0_rk * sin_theta_ib * cos_theta_ib
      self%cos_d2theta = (2.0_rk * cos_theta_ie**2 - 1.0_rk) - (2.0_rk * cos_theta_ib**2 - 1.0_rk)
    end if

  end subroutine precompute_trig_angles

  subroutine get_transonic_cone_extents(self, origin, radius)
    !< If the set of neighbor cells are transonic (some supersonic, some subsonic),
    !< then the cone needs to be extended so that it encorporates both supersonic and subsonic cells.

    class(corner_mach_cone_t), intent(inout) :: self
    real(rk), dimension(2), intent(out) :: origin  !< origin of the new cone
    real(rk), intent(out) :: radius                !< radius of the new cone

    real(rk), dimension(2) :: origin_1v3, origin_2v4
    real(rk) :: radius_1v3, radius_2v4
    real(rk), dimension(2, 2) :: origins_to_compare
    real(rk), dimension(2) :: radii_to_compare

    integer(ik) :: i
    real(rk) :: sound_speed
    logical :: is_inside
    real(rk), dimension(2) :: vel
    integer(ik) :: largest_cone_idx
    !< index of the cone with the largets radius

    real(rk), dimension(N_CELLS) :: cone_radii
    !< (cell); radius of each cone based on the neighbor cell state

    real(rk), dimension(2, N_CELLS) :: cone_origins
    !< ((x,y), cell); origin of each cone based on the neighbor cell state

    associate(u=>self%recon_state(2, :), &
              v=>self%recon_state(3, :), &
              cs=>self%sound_speed)

      ! Get the cone size based on each neighbor cell's state
      do i = 1, N_CELLS
        vel = [u(i), v(i)]
        call get_cone_extents(tau=self%tau, &
                              xy=self%p_xy, vel=vel, sound_speed=cs(i), &
                              origin=cone_origins(:, i), radius=cone_radii(i))
      end do
    end associate

    ! Since the circle comparison is only 2 at a time, this
    ! does the diagonal cells first (1v3 and 2v4) then then
    ! compares the results from those for the final circle

    ! Compare cells 1 and 3
    origins_to_compare(:, 1) = cone_origins(:, 1)
    origins_to_compare(:, 2) = cone_origins(:, 3)
    radii_to_compare = [cone_radii(1), cone_radii(3)]
    call super_circle(origins=origins_to_compare, radii=radii_to_compare, &
                      new_origin=origin_1v3, new_radius=radius_1v3)

    ! Compare cells 2 and 4
    origins_to_compare(:, 1) = cone_origins(:, 2)
    origins_to_compare(:, 2) = cone_origins(:, 4)
    radii_to_compare = [cone_radii(2), cone_radii(4)]
    call super_circle(origins=origins_to_compare, radii=radii_to_compare, &
                      new_origin=origin_2v4, new_radius=radius_2v4)

    ! Compare result from 1 vs 3 and 2 vs 4
    origins_to_compare(:, 1) = origin_1v3
    origins_to_compare(:, 2) = origin_2v4
    radii_to_compare = [radius_1v3, radius_2v4]
    call super_circle(origins=origins_to_compare, radii=radii_to_compare, &
                      new_origin=origin, new_radius=radius)
  end subroutine get_transonic_cone_extents

  subroutine determine_p_prime_cell(self)
    !< Sometimes when u and v tilde (reference velocities) are 0, P' and P are collocated. When this
    !< is the case, P' needs to be chosen from one of the neighbor cells. This subroutine scans all
    !< the neighbor cells and picks the one with the highest pressure (if any). If not, then it juse
    !< choses the first cell from the list of cells that "contain" P'

    class(corner_mach_cone_t), intent(inout) :: self
    integer(ik) :: p_prime_cell, i
    real(rk) :: max_p
    integer(ik), dimension(2) :: p_prime_ij

    max_p = 0.0_rk
    p_prime_ij = [0, 0]

    if(count(self%p_prime_in_cell) > 1) then
      do i = 1, N_CELLS
        if(self%p_prime_in_cell(i)) then
          if(self%recon_state(4, i) > max_p) then
            max_p = self%recon_state(4, i)
            p_prime_ij = self%p_prime_ij(:, i)
            p_prime_cell = i
          end if
        end if
      end do

      self%p_prime_in_cell = .false.
      self%p_prime_in_cell(p_prime_cell) = .true.
    end if
  end subroutine determine_p_prime_cell

  subroutine sanity_checks(self)
    !< Do some sanity checks to make sure the mach cone is valid
    class(corner_mach_cone_t), intent(in) :: self
    real(rk) :: arc_sum

    arc_sum = sum(self%dtheta)

    if(count(self%n_arcs_per_cell >= 2) > 1) then
      error stop "Too many arcs in the mach cone (count(cone%n_arcs >= 2) > 1)"
    end if

    if(.not. equal(arc_sum, 2 * pi, 1e-6_rk)) then
      write(std_err, *) self
      error stop "Cone arcs do not add up to 2pi"
    end if

    if(self%radius < 0.0_rk) then
      error stop "Error: cone radius < 0"
    end if

  end subroutine sanity_checks

  subroutine write_cone(self, unit, iotype, v_list, iostat, iomsg)
    !< Implementation of `write(*,*) mach_cone_base_t`

    class(corner_mach_cone_t), intent(in) :: self     !< cone class
    integer, intent(in) :: unit           !< input/output unit
    character(*), intent(in) :: iotype    !< LISTDIRECTED or DTxxx
    integer, intent(in) :: v_list(:)      !< parameters from fmt spec.
    integer, intent(out) :: iostat        !< non zero on error, etc.
    character(*), intent(inout) :: iomsg  !< define if iostat non zero.

    integer(ik) :: i

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) 'Mach Cone Details'//new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) '================='//new_line('a')

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) "location = '"//trim(self%cone_location)//"'"//new_line('a')
    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "tau = ", self%tau, new_line('a')
    write(unit, '(a, es14.6, a)', iostat=iostat, iomsg=iomsg) "radius = ", self%radius, new_line('a')
    write(unit, '(a, es14.6, ",", es14.6, a)', iostat=iostat, iomsg=iomsg) 'cone["P (x,y)"] = [', self%p_xy, "]"//new_line('a')
    write(unit, '(a, es14.6, ",", es14.6, a)', iostat=iostat, iomsg=iomsg) &
      'cone["P'//"'(x,y)"//'"] = [', self%p_prime_xy, "]"//new_line('a')
    write(unit, '(a, es14.6, ",", es14.6, a)', iostat=iostat, iomsg=iomsg) &
      "dist = [", self%p_prime_xy - self%p_xy, "]"//new_line('a')

    do i = 1, N_CELLS
      write(unit, '(a, i7,",",i7, a, i0, a)', iostat=iostat, iomsg=iomsg) &
        'cone["P'//"'(i,j)"//'"] = [', self%p_prime_ij(:, i), "] # Cell: ", i, new_line('a')
    end do

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) "# Edge Vectors: tail (x,y) -> head (x,y)"//new_line('a')

    do i = 1, N_CELLS
      write(unit, '(a, i0, a, 2(es10.3, ", ", es10.3, a))', iostat=iostat, iomsg=iomsg) "vector_", i, " = [[", &
        self%edge_vectors(:, 1, i), "], [", self%edge_vectors(:, 2, i), "] ]"//new_line('a')
    end do

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "reference_density =     ", self%reference_density, new_line('a')
    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "reference_x_velocity =  ", self%reference_u, new_line('a')
    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "reference_y_velocity =  ", self%reference_v, new_line('a')
    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "reference_sound_speed = ", self%reference_sound_speed, new_line('a')
    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "reference_mach_number = ", self%reference_mach_number, new_line('a')

    write(unit, '(a)', iostat=iostat) new_line('a')
    write(unit, '(a, 4(i0, 1x),a)', iostat=iostat) '# Valid arcs in each cell: [', self%n_arcs_per_cell, ']'
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    ! write(unit, '(a)', iostat=iostat, iomsg=iomsg) '# Angles:       Theta_ib [deg]         Theta_ie [deg]'//new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) '# Delta Theta (for each arc) [deg]'//new_line('a')
    do i = 1, N_CELLS
      ! write(unit, '(a, i0, 4(a,f9.5, ",",f9.5), a )', iostat=iostat, iomsg=iomsg) &
      !   'angles[:,', i, &
      !   '] = [[ ', rad2deg(self%theta_ib(:, i)), &
      !   '], [ ', rad2deg(self%theta_ie(:, i)), &
      !   ']]'//new_line('a')
      write(unit, '(a, i0,a,f9.5,",",f9.5, a)', iostat=iostat, iomsg=iomsg) &
        'dtheta[', i, '] = [', rad2deg(self%dtheta(:, i)), ']'//new_line('a')
    end do
    write(unit, '(a, f0.2, a)', iostat=iostat, iomsg=iomsg) 'arc_length_sum = ', &
      rad2deg(sum(self%dtheta)), ' '//new_line('a')

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) 'Cell Primitive Vars State [rho,u,v,p]'//new_line('a')
    do i = 1, N_CELLS
      write(unit, '(a, i0, a, 4(es10.3, 1x), a)', iostat=iostat, iomsg=iomsg) &
        'Cell: ', i, ' [ ', self%arc_primitive_vars(:, 1, i), '] (arc 1)'//new_line('a')
      write(unit, '(a, 4(es10.3, 1x), a)', iostat=iostat, iomsg=iomsg) &
        '        [ ', self%arc_primitive_vars(:, 2, i), '] (arc 2)'//new_line('a')
    end do

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) "# P' in cell?"//new_line('a')
    do i = 1, N_CELLS
      write(unit, '(a, i0, a, l2, a)', iostat=iostat, iomsg=iomsg) 'Cell: ', i, ' [', self%p_prime_in_cell(i), ' ]'//new_line('a')
    end do

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) '# Neighbor cells contributing to the mach cone '//new_line('a')
    write(unit, '(a, i0, a)', iostat=iostat, iomsg=iomsg) 'neighbor_cells = ', self%n_neighbor_cells, new_line('a')
    do i = 1, N_CELLS
      write(unit, '(a, i0, " = [", 3(es10.3, ","), es10.3,"]", a)', iostat=iostat, iomsg=iomsg) &
        'cell_', i, self%recon_state(:, i), new_line('a')
    end do

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
  end subroutine write_cone
end module mod_corner_mach_cone
