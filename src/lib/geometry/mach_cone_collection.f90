module mod_mach_cone_collection
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use, intrinsic :: ieee_arithmetic
  use math_constants, only: pi, rad2deg
  use mod_globals, only: TINY_MACH, TINY_DIST, grid_is_orthogonal
  use mod_floating_point_utils, only: near_zero, equal
  use mod_eos, only: eos
  use mod_vector, only: vector_t
  use mod_geometry
  use mod_mach_cone_utilties, only: tau_scaling, enable_tau_scaling

  implicit none

  private
  public :: mach_cone_collection_t

  real(rk) :: tau_param = 1e-16_rk

  type :: mach_cone_collection_t
    ! TODO: midpoint cones can't have more than 1 arc!!
    real(rk), dimension(:, :, :, :), allocatable :: recon_state
    !< ((rho, u, v, p), (cell 1:N), i, j); Reconstructed value of U at P for each neighbor cell

    ! logical, dimension(:, :), allocatable :: cell_is_supersonic
    !< (i, j); Is each neighbor cell supersonic or not? This is needed for transonic mach cones in
    !< order to apply an entropy fix for transonic rarefaction regions.

    ! Note: The cell/arc indexing is a bit tricky, but there is a good reason why it is so...
    ! Take self%dtheta for example; it's indexed as ((cell_1_arc_1:cell_N_arc_2), i, j). The primary
    ! reason for making it this way is so that the E0 operator can directly sum along dim=1 and make
    ! the math/code easy to read as well as have nice memory access patterns (all data is next to each other in sum)

    ! So, for a corner mach cone with 4 cells, the indexing goes like the following
    ! cell: 1, arc: 1 => array index = (1, i, j)
    ! cell: 1, arc: 2 => array index = (2, i, j)
    ! cell: 2, arc: 1 => array index = (3, i, j)
    ! cell: 2, arc: 2 => array index = (4, i, j)
    ! cell: 3, arc: 1 => array index = (5, i, j)
    ! cell: 3, arc: 2 => array index = (6, i, j)
    ! cell: 4, arc: 1 => array index = (7, i, j)
    ! cell: 4, arc: 2 => array index = (8, i, j)

    ! When looping through these arrays use the following:

    ! when the arc index is first
    ! do arc = 1, 2
    !   do cell = 1, 4
    !     idx = cell + (arc-1) * ncells ! ncells = 4 or 2
    !     array(idx, i, j) = something...
    !   end do
    ! end do

    ! when the cell index is first
    ! do cell = 1, 4
    !   do arc = 1, 2
    !     idx = arc + (cell-1)*2
    !     array(idx, i, j) = something...
    !   end do
    ! end do

    real(rk), dimension(:, :, :), allocatable :: dtheta
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); theta_ie - theta_ib [radians]

    real(rk), dimension(:, :, :), allocatable :: sin_dtheta
    !< (cell_1_arc_1:cell_N_arc_2); sin(theta_ie) - sin(theta_ib)

    real(rk), dimension(:, :, :), allocatable :: cos_dtheta
    !< (cell_1_arc_1:cell_N_arc_2); cos(theta_ie) - cos(theta_ib)

    real(rk), dimension(:, :, :), allocatable :: sin_d2theta
    !< (cell_1_arc_1:cell_N_arc_2); sin(2*theta_ie) - sine(2*theta_ib)

    real(rk), dimension(:, :, :), allocatable :: cos_d2theta
    !< (cell_1_arc_1:cell_N_arc_2); cos(2*theta_ie) - cos(2*theta_ib)

    integer(ik), dimension(:, :, :), allocatable :: p_prime_ij
    !< ((P' i, P' j), i, j) i,j location of P' or the apex of the Mach cone (global, not relative to P0)

    real(rk), dimension(:, :, :), allocatable :: p_prime_vector
    !< ((x,y), i, j) vector from P0 to P' (P0 is shifted to (0,0))

    integer(ik), dimension(:, :, :), allocatable :: n_arcs_per_cell
    !< (cell 1:N); Number of arcs that each cell contributes (0, 1, or 2)

    logical, dimension(:, :, :), allocatable :: p_prime_in_cell
    !< ((cell 1:N), i, j); is P' inside the control volume?

    real(rk), dimension(:, :), allocatable :: pressure_p_prime
    !< (i, j); pressure from the P' cell (due to FVLEG assumptions, this is the reconstructed value at P0 of the cell that contains P')

    real(rk), dimension(:, :), allocatable :: density_p_prime
    !< density from the P' cell (due to FVLEG assumptions, this is the reconstructed value at P0 of the cell that contains P')

    real(rk), dimension(:, :, :), allocatable :: rho
    !< (cell_1_arc_1:cell_N_arc_2); pressure of the cell along the are; Note the dimension (cell 1:N)*2 , where 2 means # arcs
    real(rk), dimension(:, :, :), allocatable :: u
    !< (cell_1_arc_1:cell_N_arc_2); u velocity of the cell along the are; Note the dimension (cell 1:N)*2 , where 2 means # arcs
    real(rk), dimension(:, :, :), allocatable :: v
    !< (cell_1_arc_1:cell_N_arc_2); v velocity of the cell along the are; Note the dimension (cell 1:N)*2 , where 2 means # arcs
    real(rk), dimension(:, :, :), allocatable :: p
    !< (cell_1_arc_1:cell_N_arc_2); pressure of the cell along the are; Note the dimension (cell 1:N)*2 , where 2 means # arcs

    integer(ik), dimension(:, :, :, :), allocatable :: cell_indices
    !< ((i,j), cell 1:N, i, j); indices of the neighboring cells

    integer(ik) :: n_neighbor_cells    !< Number of neighbor cells
    real(rk), dimension(:, :), allocatable :: radius                  !< Radius of the cone
    real(rk), dimension(:, :), allocatable :: tau                     !< Time evolution increment
    real(rk), dimension(:, :), allocatable :: reference_density       !< Reference density (e.g. neighbor averaged)
    real(rk), dimension(:, :), allocatable :: reference_u             !< Reference x velocity (e.g. neighbor averaged)
    real(rk), dimension(:, :), allocatable :: reference_v             !< Reference y velocity (e.g. neighbor averaged)
    real(rk), dimension(:, :), allocatable :: reference_sound_speed   !< Reference sound speed (e.g. neighbor averaged)
    logical, dimension(:, :), allocatable:: cone_is_transonic        !< Flag to enable special treatment for transonic cones
    logical, dimension(:, :), allocatable:: cone_is_centered         !< is the cone located at P0?
    character(len=32) :: cone_location  !< Corner or midpoint cone?
    integer(ik) :: ni !< # of cones in the i direction
    integer(ik) :: nj !< # of cones in the j direction
  contains
    procedure, public :: initialize
    procedure, private :: get_reference_state
    ! procedure, private :: determine_p_prime_cell
    procedure, private :: sanity_checks
    procedure, private :: find_arc_angles
    procedure, private :: find_arc_segments
    procedure, private :: compute_trig_angles
    ! procedure, private :: get_transonic_cone_extents
    ! procedure, pass :: write => write_cone
    ! generic, public :: write(formatted) => write
  end type mach_cone_collection_t

contains

  subroutine initialize(self, edge_vectors, vector_scaling, reconstructed_state, cell_indices, cone_location, tau)
    !< Class constructor

    class(mach_cone_collection_t), intent(inout) :: self

    character(len=*), intent(in) :: cone_location !< corner, down/up midpoint, left/right midpoint

    real(rk), intent(in), optional :: tau

    real(rk), dimension(:, :, :, :), intent(in) :: edge_vectors
    !< ((x,y), (vector_1:N), i, j) ; set of vectors that define the corner

    real(rk), dimension(:, :), intent(in) :: vector_scaling
    !< (i, j); scale factor for each set of edge vectors (this is the min length of the set). This is used to
    !< scale the P' vector down to the same frame as the edge vectors

    integer(ik), dimension(:, :, :, :), intent(in) :: cell_indices
    !< ((i,j), (cell_1:N), i, j); set of indices for the neighboring cells -> needed to find P' i,j index

    real(rk), dimension(:, :, :, :), intent(in) :: reconstructed_state
    !< ((rho,u,v,p), (cell_1:N), i, j); reconstructed state for point P.

    integer(ik) :: idx, arc, c, i, j, ni, nj

    self%ni = size(reconstructed_state, dim=3)
    self%nj = size(reconstructed_state, dim=4)
    self%cone_location = trim(cone_location)

    if(present(tau)) then
      tau_param = tau
    end if

    if(any(vector_scaling < 0.0_rk)) error stop "vector_scaling < 0"

    print *, 'min/max vector scale', minval(vector_scaling), maxval(vector_scaling)

    select case(trim(cone_location))
    case('corner')
      self%n_neighbor_cells = 4
    case('down/up midpoint', 'left/right midpoint')
      self%n_neighbor_cells = 2
    case default
      error stop "Error in mach_cone_collection_t%initialize(), unsupported cone location"
    end select

    do j = 1, self%nj
      do i = 1, self%ni
        do c = 1, self%n_neighbor_cells
          if(reconstructed_state(1, c, i, j) < 0.0_rk) then
            write(std_err, '(a, 4(es10.3,1x))') 'Reconstructed density states [' // trim(self%cone_location) //'] (cell 1:N)', reconstructed_state(1, :, i, j)
            write(std_err, '(a, i0, 1x, i0, a)') 'Cone index [i, j]: ', i, j
            write(std_err, '(a, 8(i0, 1x))') 'Cone neighbor cell indices [i, j]: ', cell_indices(:, :, i, j)
            error stop "Error in corner_mach_cone_t initialization, density in the reconstructed state is < 0"
          end if

          if(reconstructed_state(4, c, i, j) < 0.0_rk) then
            write(std_err, '(a, 4(es10.3,1x))') 'Reconstructed pressure states [' // trim(self%cone_location) //'] (cell 1:N): ', reconstructed_state(4, :, i, j)
            write(std_err, '(a, i0, 1x, i0, a)') 'Cone index [i, j]: [', i, j, ']'
            write(std_err, '(a, 8(i0, 1x))') 'Cone neighbor cell indices [i, j]: ', cell_indices(:, :, i, j)
            error stop "Error in corner_mach_cone_t initialization, pressure in the reconstructed state is < 0"
          end if
        end do
      end do
    end do

    associate(nc=>self%n_neighbor_cells, ni=>self%ni, nj=>self%nj)
      ! allocate(self%cell_is_supersonic(ni, nj))
      ! self%cell_is_supersonic = .false.

      if(.not. allocated(self%p_prime_in_cell)) allocate(self%p_prime_in_cell(nc, ni, nj))
      self%p_prime_in_cell = .false.

      if(.not. allocated(self%recon_state)) allocate(self%recon_state(4, nc, ni, nj))
      self%recon_state = reconstructed_state

      if(.not. allocated(self%dtheta)) allocate(self%dtheta(nc * 2, ni, nj))
      self%dtheta = 0.0_rk

      if(.not. allocated(self%sin_dtheta)) allocate(self%sin_dtheta(nc * 2, ni, nj))
      self%sin_dtheta = 0.0_rk

      if(.not. allocated(self%cos_dtheta)) allocate(self%cos_dtheta(nc * 2, ni, nj))
      self%cos_dtheta = 0.0_rk

      if(.not. allocated(self%sin_d2theta)) allocate(self%sin_d2theta(nc * 2, ni, nj))
      self%sin_d2theta = 0.0_rk

      if(.not. allocated(self%cos_d2theta)) allocate(self%cos_d2theta(nc * 2, ni, nj))
      self%cos_d2theta = 0.0_rk

      if(.not. allocated(self%p_prime_vector)) allocate(self%p_prime_vector(2, ni, nj))
      self%p_prime_vector = 0.0_rk

      if(.not. allocated(self%pressure_p_prime)) allocate(self%pressure_p_prime(ni, nj))
      self%pressure_p_prime = 0.0_rk

      if(.not. allocated(self%density_p_prime)) allocate(self%density_p_prime(ni, nj))
      self%density_p_prime = 0.0_rk

      if(.not. allocated(self%rho)) allocate(self%rho(nc * 2, ni, nj))
      self%rho = 0.0_rk
      if(.not. allocated(self%u)) allocate(self%u(nc * 2, ni, nj))
      self%u = 0.0_rk
      if(.not. allocated(self%v)) allocate(self%v(nc * 2, ni, nj))
      self%v = 0.0_rk
      if(.not. allocated(self%p)) allocate(self%p(nc * 2, ni, nj))
      self%p = 0.0_rk

      if(.not. allocated(self%p_prime_ij)) allocate(self%p_prime_ij(2, ni, nj))

      if(.not. allocated(self%n_arcs_per_cell)) allocate(self%n_arcs_per_cell(nc, ni, nj))
      self%n_arcs_per_cell = 0

      if(.not. allocated(self%cell_indices)) allocate(self%cell_indices(2, nc, ni, nj))
      self%cell_indices = cell_indices

      if(.not. allocated(self%radius)) allocate(self%radius(ni, nj))
      if(.not. allocated(self%tau)) allocate(self%tau(ni, nj))
      if(.not. allocated(self%reference_density)) allocate(self%reference_density(ni, nj))
      if(.not. allocated(self%reference_u)) allocate(self%reference_u(ni, nj))
      if(.not. allocated(self%reference_v)) allocate(self%reference_v(ni, nj))
      if(.not. allocated(self%reference_sound_speed)) allocate(self%reference_sound_speed(ni, nj))

      if(.not. allocated(self%cone_is_transonic)) allocate(self%cone_is_transonic(ni, nj))
      self%cone_is_transonic = .false.

      if(.not. allocated(self%cone_is_centered)) allocate(self%cone_is_centered(ni, nj))
      self%cone_is_centered = .false.
    end associate

    call self%get_reference_state()

    ! Set tau (infinitesimal time increment)
    if(enable_tau_scaling .and. .not. present(tau)) then
      self%tau = tau_scaling / self%reference_sound_speed
    else
      self%tau = tau_param
    end if

    print *, 'min/max tau', minval(self%tau), maxval(self%tau)

    ! Set the P' vector. Remember that all edge vector sets are moved to the origin, so P0 is at (0,0)
    do j = 1, self%nj
      do i = 1, self%ni
        ! In the FVLEG paper, P'(x,y) = [x - tau * u_tilde, y - tau * v_tilde], but P0 is at (0,0) below
        self%p_prime_vector(:, i, j) = -self%tau(i, j) * vector_scaling(i, j)* &
                                       [self%reference_u(i, j), self%reference_v(i, j)]
      end do
    end do

    self%radius = self%tau * self%reference_sound_speed

    ! For very small distances, just make P' coincide with P0 at the origin
    do j = 1, self%nj
      do i = 1, self%ni
        if(abs(self%p_prime_vector(1, i, j)) < TINY_DIST) self%p_prime_vector(1, i, j) = 0.0_rk
        if(abs(self%p_prime_vector(2, i, j)) < TINY_DIST) self%p_prime_vector(2, i, j) = 0.0_rk
      end do
    end do

    call self%compute_trig_angles(edge_vectors=edge_vectors)

    ! Assign the reconstructed quantities to the primitive var values for each arc of the mach cone
    do j = 1, self%nj
      do i = 1, self%ni
        do arc = 1, 2
          do c = 1, self%n_neighbor_cells

            idx = c + (arc - 1) * self%n_neighbor_cells ! convert indexing for better storage

            ! Normalize for numbers closer to 1 (better for floating point storage)
            self%rho(idx, i, j) = self%recon_state(1, c, i, j)
            self%u(idx, i, j) = self%recon_state(2, c, i, j) / self%reference_sound_speed(i, j)
            self%v(idx, i, j) = self%recon_state(3, c, i, j) / self%reference_sound_speed(i, j)
            self%p(idx, i, j) = self%recon_state(4, c, i, j) / &
                                (self%reference_density(i, j) * self%reference_sound_speed(i, j)**2)
          end do
        end do
      end do
    end do

    do j = 1, self%nj
      do i = 1, self%ni
        do c = 1, self%n_neighbor_cells
          if(self%p_prime_in_cell(c, i, j)) then
            self%density_p_prime(i, j) = self%recon_state(1, c, i, j)
            self%pressure_p_prime(i, j) = self%recon_state(4, c, i, j)
          end if
        end do
      end do
    end do

    ! call self%determine_p_prime_cell()
    call self%sanity_checks()
  end subroutine initialize

  subroutine get_reference_state(self)
    !< Find the reference state of the cone only based on the cells that are touched by the cone
    class(mach_cone_collection_t), intent(inout) :: self

    real(rk), dimension(:, :, :), allocatable :: mach_number, sound_speed
    integer(ik), dimension(:, :), allocatable :: n_supersonic_cells
    integer(ik) :: i, j, c
    real(rk) :: gamma, mach_u, mach_v

    gamma = eos%get_gamma()

    allocate(mach_number(self%n_neighbor_cells, self%ni, self%nj))
    allocate(sound_speed(self%n_neighbor_cells, self%ni, self%nj))
    allocate(n_supersonic_cells(self%ni, self%nj))
    n_supersonic_cells = 0

    sound_speed = sqrt(gamma * self%recon_state(4, :, :, :) / self%recon_state(1, :, :, :))
    mach_number = sqrt(self%recon_state(2, :, :, :)**2 + self%recon_state(3, :, :, :)**2) / sound_speed

    do j = 1, self%nj
      do i = 1, self%ni
        n_supersonic_cells(i, j) = 0
        do c = 1, self%n_neighbor_cells
          if(mach_number(c, i, j) > 1.0_rk) then
            n_supersonic_cells(i, j) = n_supersonic_cells(i, j) + 1
          end if
        end do
      end do
    end do

    where(n_supersonic_cells == self%n_neighbor_cells .or. n_supersonic_cells == 0)
    self%cone_is_transonic = .false.
    else where
    self%cone_is_transonic = .true.
    end where

    ! Use the mach number to eliminate small velocities
    ! (this scales with the problem better than an absolute velocity value)
    ! do j = 1, self%nj
    !   do i = 1, self%ni
    !     do c = 1, self%n_neighbor_cells
    !       mach_u = abs(self%recon_state(2, c, i, j) / sound_speed(c, i, j))
    !       mach_v = abs(self%recon_state(3, c, i, j) / sound_speed(c, i, j))
    !       if(mach_u < TINY_MACH) self%recon_state(2, c, i, j) = 0.0_rk
    !       if(mach_v < TINY_MACH) self%recon_state(3, c, i, j) = 0.0_rk
    !       if(mach_u < TINY_MACH .and. mach_v < TINY_MACH) then
    !         self%cone_is_centered(i, j) = .true.
    !       else
    !         self%cone_is_centered(i, j) = .false.
    !       end if
    !     end do
    !   end do
    ! end do

    self%reference_density = sum(self%recon_state(1, :, :, :), dim=1) / real(self%n_neighbor_cells, rk)
    self%reference_u = sum(self%recon_state(2, :, :, :), dim=1) / real(self%n_neighbor_cells, rk)
    self%reference_v = sum(self%recon_state(3, :, :, :), dim=1) / real(self%n_neighbor_cells, rk)

    self%reference_sound_speed = sqrt(gamma * (sum(self%recon_state(4, :, :, :), dim=1) / real(self%n_neighbor_cells, rk)) &
                                      / self%reference_density)

    deallocate(mach_number)
    deallocate(sound_speed)
    deallocate(n_supersonic_cells)
  end subroutine get_reference_state

  subroutine compute_trig_angles(self, edge_vectors)
    !< Precompute some of the trig values, like sin(theta_ie), etc. to save compute time. If
    !< the grid is orthogonal and the references state is at rest, then the math
    !< gets even easier

    class(mach_cone_collection_t), intent(inout) :: self

    real(rk), dimension(:, :, :, :), intent(in) :: edge_vectors

    real(rk), dimension(:, :, :), allocatable :: theta_ib
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); starting angle [rad] for the arc contained in each cell

    real(rk), dimension(:, :, :), allocatable:: theta_ie
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); ending angle [rad] for the arc contained in each cell

    real(rk), dimension(:, :, :), allocatable :: sin_theta_ib
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); sin(theta_ib) for the arc contained in each cell

    real(rk), dimension(:, :, :), allocatable :: cos_theta_ib
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); cos(theta_ib) for the arc contained in each cell

    real(rk), dimension(:, :, :), allocatable :: sin_theta_ie
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); sin(theta_ie) for the arc contained in each cell

    real(rk), dimension(:, :, :), allocatable :: cos_theta_ie
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); cos(theta_ie) for the arc contained in each cell

    integer(ik) :: i, j, c, arc, max_n_arcs
    integer(ik) :: idx
    !< index conversion to go from ((arc1, arc2), (cell 1:N_CELLS), i, j) to ((cell_1_arc_1, cell_N_arc_1:cell_1_arc_N, cell_N_arc_N), i, j)

    integer(ik), dimension(:, :, :), allocatable :: n_arcs_per_cell

    ! allocate(theta_ib(2*self%n_neighbor_cells, self%ni, self%nj))
    ! theta_ib = 0.0_rk

    ! allocate(theta_ie(2*self%n_neighbor_cells, self%ni, self%nj))
    ! theta_ie = 0.0_rk

    allocate(sin_theta_ib(2 * self%n_neighbor_cells, self%ni, self%nj))
    sin_theta_ib = 0.0_rk

    allocate(cos_theta_ib(2 * self%n_neighbor_cells, self%ni, self%nj))
    cos_theta_ib = 0.0_rk

    allocate(sin_theta_ie(2 * self%n_neighbor_cells, self%ni, self%nj))
    sin_theta_ie = 0.0_rk

    allocate(cos_theta_ie(2 * self%n_neighbor_cells, self%ni, self%nj))
    cos_theta_ie = 0.0_rk

    allocate(n_arcs_per_cell(self%n_neighbor_cells, self%ni, self%nj))
    n_arcs_per_cell = 0

    ! If we know the cone is at the origin and the grid is orthogonal, then
    ! the angles are easily known ahead of time, so we can skip the additional
    ! math needed to find the intersection angles
    ! if(grid_is_orthogonal) then
    !   do j = 1, self%nj
    !     do i = 1, self%ni
    !       if(self%cone_is_centered(i, j)) then

    !         self%dtheta(1:4, i, j) = pi / 2.0_rk ! all of the 1st arcs are a quarter of the circle
    !         self%n_arcs_per_cell(:, i, j) = 1
    !         self%p_prime_in_cell(:, i, j) = [.true., .false., .false., .false.]
    !         ! Cell 1 (lower left):  180 -> 270
    !         self%sin_dtheta(1, i, j) = -1.0_rk
    !         self%cos_dtheta(1, i, j) = 1.0_rk
    !         self%sin_d2theta(1, i, j) = 0.0_rk
    !         self%cos_d2theta(1, i, j) = -2.0_rk

    !         ! Cell 2 (lower right): 270 -> 360
    !         self%sin_dtheta(2, i, j) = 1.0_rk
    !         self%cos_dtheta(2, i, j) = 1.0_rk
    !         self%sin_d2theta(2, i, j) = 0.0_rk
    !         self%cos_d2theta(2, i, j) = 2.0_rk

    !         ! Cell 3 (upper right): 0 -> 90
    !         self%sin_dtheta(3, i, j) = 1.0_rk
    !         self%cos_dtheta(3, i, j) = -1.0_rk
    !         self%sin_d2theta(3, i, j) = 0.0_rk
    !         self%cos_d2theta(3, i, j) = -2.0_rk

    !         ! Cell 4 (upper left):  90 -> 180
    !         self%sin_dtheta(4, i, j) = -1.0_rk
    !         self%cos_dtheta(4, i, j) = -1.0_rk
    !         self%sin_d2theta(4, i, j) = 0.0_rk
    !         self%cos_d2theta(4, i, j) = 2.0_rk
    !       end if
    !     end do
    !   end do
    ! end if

    call self%find_arc_angles(edge_vectors, theta_ib, theta_ie, n_arcs_per_cell)

    ! Precompute for a bit of speed
    self%n_arcs_per_cell = n_arcs_per_cell
    max_n_arcs = maxval(n_arcs_per_cell)
    do j = 1, self%nj
      do i = 1, self%ni
        do arc = 1, 2
          do c = 1, self%n_neighbor_cells
            idx = c + (arc - 1) * self%n_neighbor_cells ! convert indexing for better storage
            self%dtheta(idx, i, j) = abs(theta_ie(idx, i, j) - theta_ib(idx, i, j))
            sin_theta_ib(idx, i, j) = sin(theta_ib(idx, i, j))
            cos_theta_ib(idx, i, j) = cos(theta_ib(idx, i, j))
            sin_theta_ie(idx, i, j) = sin(theta_ie(idx, i, j))
            cos_theta_ie(idx, i, j) = cos(theta_ie(idx, i, j))
          end do
        end do
      end do
    end do

    where(ieee_is_nan(sin_theta_ib)) sin_theta_ib = 0.0_rk
    where(ieee_is_nan(cos_theta_ib)) cos_theta_ib = 0.0_rk
    where(ieee_is_nan(sin_theta_ie)) sin_theta_ie = 0.0_rk
    where(ieee_is_nan(cos_theta_ie)) cos_theta_ie = 0.0_rk

    self%sin_dtheta = sin_theta_ie - sin_theta_ib
    self%cos_dtheta = cos_theta_ie - cos_theta_ib
    where(abs(self%sin_dtheta) < 1e-10_rk) self%sin_dtheta = 0.0_rk
    where(abs(self%cos_dtheta) < 1e-10_rk) self%cos_dtheta = 0.0_rk

    self%sin_d2theta = 2.0_rk * sin_theta_ie * cos_theta_ie - 2.0_rk * sin_theta_ib * cos_theta_ib
    self%cos_d2theta = (2.0_rk * cos_theta_ie**2 - 1.0_rk) - (2.0_rk * cos_theta_ib**2 - 1.0_rk)
    where(abs(self%sin_d2theta) < 1e-10_rk) self%sin_d2theta = 0.0_rk
    where(abs(self%cos_d2theta) < 1e-10_rk) self%cos_d2theta = 0.0_rk

    deallocate(sin_theta_ib)
    deallocate(cos_theta_ib)
    deallocate(sin_theta_ie)
    deallocate(cos_theta_ie)
    deallocate(n_arcs_per_cell)

  end subroutine compute_trig_angles

  subroutine find_arc_angles(self, edge_vectors, theta_ib, theta_ie, n_arcs_per_cell)
    !< Find the starting and ending angle of each arc that is held within each cell

    class(mach_cone_collection_t), intent(inout) :: self

    real(rk), dimension(:, :, :, :), intent(in) :: edge_vectors
    !< ((x,y), (vector_1:N), i, j) ; set of vectors that define the corner

    real(rk), dimension(:, :, :), allocatable :: theta_ib
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); starting angle [rad] for the arc contained in each cell

    real(rk), dimension(:, :, :), allocatable:: theta_ie
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); ending angle [rad] for the arc contained in each cell

    integer(ik), dimension(:, :, :), allocatable, intent(out) :: n_arcs_per_cell !< # of arcs in a given cell

    real(rk), dimension(2, 2) :: theta_ib_ie
    !< ((start,end), (arc1, arc2)); start/end angle for a given cell

    real(rk), dimension(2) :: second_vector, first_vector
    real(rk), dimension(2, 4) :: corner_vector_set
    real(rk), dimension(2, 2) :: midpoint_vector_set
    integer(ik) :: i, j, c, idx
    integer(ik) :: arc  !< index for looping through the # arcs in each neighbor cell
    integer(ik) :: n_arcs

    real(rk) :: p_prime_cross_vec_1
    real(rk) :: vec_2_cross_p_prime

    allocate(theta_ib(2 * self%n_neighbor_cells, self%ni, self%nj))
    theta_ib = 0.0_rk
    allocate(theta_ie(2 * self%n_neighbor_cells, self%ni, self%nj))
    theta_ie = 0.0_rk
    allocate(n_arcs_per_cell(self%n_neighbor_cells, self%ni, self%nj))
    n_arcs_per_cell = 0

    p_prime_cross_vec_1 = 0.0_rk
    vec_2_cross_p_prime = 0.0_rk

    second_vector = 0.0_rk
    first_vector = 0.0_rk
    corner_vector_set = 0.0_rk
    midpoint_vector_set = 0.0_rk

    select case(trim(self%cone_location))
    case('corner')
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
      do j = 1, self%nj
        do i = 1, self%ni
          corner_vector_set = edge_vectors(:, 1:4, i, j)

          do c = 1, 4
            select case(c)
            case(1) ! Cell 1
              first_vector = corner_vector_set(:, 4)
              second_vector = corner_vector_set(:, 1)
            case(2) ! Cell 2
              first_vector = corner_vector_set(:, 1)
              second_vector = corner_vector_set(:, 2)
            case(3) ! Cell 3
              first_vector = corner_vector_set(:, 2)
              second_vector = corner_vector_set(:, 3)
            case(4) ! Cell 4
              first_vector = corner_vector_set(:, 3)
              second_vector = corner_vector_set(:, 4)
            end select

            ! Determine if P' is in this cell
            p_prime_cross_vec_1 = self%p_prime_vector(1, i, j) * first_vector(2) - self%p_prime_vector(2, i, j) * first_vector(1)
            vec_2_cross_p_prime = second_vector(1) * self%p_prime_vector(2, i, j) - second_vector(2) * self%p_prime_vector(1, i, j)

            if(p_prime_cross_vec_1 <= 0.0_rk .and. vec_2_cross_p_prime <= 0.0_rk) then
              self%p_prime_in_cell(c, i, j) = .true.
              self%p_prime_ij(:, i, j) = [i, j]
            end if

            call self%find_arc_segments(vector_1=first_vector, vector_2=second_vector, &
                                        origin_in_cell=self%p_prime_in_cell(c, i, j), &
                                        circle_xy=self%p_prime_vector(:, i, j), &
                                        circle_radius=self%radius(i, j), &
                                        arc_segments=theta_ib_ie, n_arcs=n_arcs)

            n_arcs_per_cell(c, i, j) = n_arcs

            do arc = 1, 2
              idx = arc + (c - 1) * 2
              theta_ib(idx, i, j) = theta_ib_ie(1, arc)
              theta_ie(idx, i, j) = theta_ib_ie(2, arc)
            end do

          end do
        end do
      end do
    case('down/up midpoint', 'left/right midpoint')

      ! down/up midpoint
      !       cell 1          edge          cell 2
      !      (i-1,j)         vector          (i,j)
      !  N4----M3----N3       P2       N4----M3----N3
      !  |           |        |       |            |
      !  M4    C1    M2       O       M4    C2     M2
      !  |           |        |       |            |
      !  N1----M1----N2      P1       N1----M1----N2

      ! left/right midpoint
      !       cell 1
      !      (i-1,j)
      !  N4----M3----N3
      !  |            |
      !  M4    C1     M2
      !  |            |
      !  N1----M1----N2
      !
      !  P1----O-----P2  edge vector
      !
      !  N4----M3----N3
      !  |           |
      !  M4    C2    M2
      !  |           |
      !  N1----M1----N2
      !       cell 2
      !      (i,j-1)

      ! the order of the vectors matters, mainly for the cross product to determine
      ! if P' is in the neighboring cell or not
      ! edge_vector_ordering = [2, 1] [1, 2]

      do j = 1, self%nj
        do i = 1, self%ni
          midpoint_vector_set = edge_vectors(:, 1:2, i, j)

          do c = 1, 2
            select case(c)
            case(1) ! Cell 1
              first_vector = midpoint_vector_set(:, 2)
              second_vector = midpoint_vector_set(:, 1)
            case(2) ! Cell 2
              first_vector = midpoint_vector_set(:, 1)
              second_vector = midpoint_vector_set(:, 2)
            end select

            ! Determine if P' is in this cell
            p_prime_cross_vec_1 = self%p_prime_vector(1, i, j) * first_vector(2) - &
                                  self%p_prime_vector(2, i, j) * first_vector(1)
            vec_2_cross_p_prime = second_vector(1) * self%p_prime_vector(2, i, j) - &
                                  second_vector(2) * self%p_prime_vector(1, i, j)

            if(p_prime_cross_vec_1 <= 0.0_rk .and. vec_2_cross_p_prime <= 0.0_rk) then
              self%p_prime_in_cell(c, i, j) = .true.
              self%p_prime_ij(:, i, j) = [i, j]
            end if

            call self%find_arc_segments(vector_1=first_vector, vector_2=second_vector, &
                                        origin_in_cell=self%p_prime_in_cell(c, i, j), &
                                        circle_xy=self%p_prime_vector(:, i, j), &
                                        circle_radius=self%radius(i, j), &
                                        arc_segments=theta_ib_ie, n_arcs=n_arcs)

            n_arcs_per_cell(c, i, j) = n_arcs

            do arc = 1, 2
              idx = arc + (c - 1) * 2
              theta_ib(idx, i, j) = theta_ib_ie(1, arc)
              theta_ie(idx, i, j) = theta_ib_ie(2, arc)
            end do

          end do
        end do
      end do

    case default
      error stop 'Error in mach_cone_collection_t%find_arc_angles(), invalid cone location'
    end select

  end subroutine find_arc_angles

  pure subroutine find_arc_segments(self, vector_1, vector_2, origin_in_cell, &
                                    circle_xy, circle_radius, arc_segments, n_arcs)
    !< Given 2 lines and a circle, find their intersections and starting/ending angles for each arc

    ! Input
    ! real(rk), dimension(2, 2, 2), intent(in) :: lines
    !< ((x,y), (tail,head), (line_1:line_2)); set of vectors that define the lines to intersect with the circle
    class(mach_cone_collection_t), intent(inout) :: self

    real(rk), dimension(2), intent(in) :: vector_1 !< (x,y) 1st vector
    real(rk), dimension(2), intent(in) :: vector_2 !< (x,y) 2nd vector

    logical, intent(in) :: origin_in_cell  !< is the circle origin in the cell defined by the two lines?
    real(rk), dimension(2), intent(in) :: circle_xy
    real(rk), intent(in) :: circle_radius                 !< radius of the circle

    ! Output
    real(rk), dimension(2, 2), intent(out):: arc_segments !< ((start,end), (arc1, arc2))
    integer(ik), intent(out) :: n_arcs !< Number of arcs with a valid start and end angle

    ! Dummy
    integer(ik) :: n_intersections_per_line !< # intersections for a single line
    integer(ik), dimension(2) :: n_intersections !< # intersections for all lines
    real(rk), dimension(2) :: intersection_angles_per_line !< intersection angles
    real(rk), dimension(2, 2) :: intersection_angles !< intersection angles for all lines
    integer(ik) :: i
    logical, dimension(2) :: valid_intersections
    logical, dimension(2, 2) :: total_valid_intersections !< ((intersection 1, intersection 2), (line 1, line 2)); valid intersections for each line

    valid_intersections = .false.
    total_valid_intersections = .false.
    n_intersections_per_line = 0
    n_intersections = 0
    intersection_angles_per_line = 0.0_rk
    intersection_angles = 0.0_rk

    ! Line 1
    ! For a given line & circle intersection, find the angle with respect to the x-axis for each intersection point
    call get_origin_intersection_angles(line_xy=vector_1, circle_xy=circle_xy, circle_radius=circle_radius, &
                                        arc_angles=intersection_angles_per_line, valid_intersections=valid_intersections)
    intersection_angles(1, :) = intersection_angles_per_line
    n_intersections(1) = count(valid_intersections)
    total_valid_intersections(:, 1) = valid_intersections

    ! Line 2
    ! For a given line & circle intersection, find the angle with respect to the x-axis for each intersection point
    call get_origin_intersection_angles(line_xy=vector_2, circle_xy=circle_xy, circle_radius=circle_radius, &
                                        arc_angles=intersection_angles_per_line, valid_intersections=valid_intersections)
    intersection_angles(2, :) = intersection_angles_per_line
    n_intersections(2) = count(valid_intersections)
    total_valid_intersections(:, 2) = valid_intersections

    ! Find arc starting and ending angles
    call get_theta_start_end(thetas=intersection_angles, origin_in_cell=origin_in_cell, &
                             valid_intersections=total_valid_intersections, &
                             n_intersections=n_intersections, &
                             theta_start_end=arc_segments, n_arcs=n_arcs)
  end subroutine find_arc_segments

  ! subroutine get_transonic_cone_extents(self, origin, radius)
  !   !< If the set of neighbor cells are transonic (some supersonic, some subsonic),
  !   !< then the cone needs to be extended so that it encorporates both supersonic and subsonic cells.

  !   class(corner_mach_cone_t), intent(inout) :: self
  !   real(rk), dimension(2), intent(out) :: origin  !< origin of the new cone
  !   real(rk), intent(out) :: radius                !< radius of the new cone

  !   real(rk), dimension(2) :: origin_1v3, origin_2v4
  !   real(rk) :: radius_1v3, radius_2v4
  !   real(rk), dimension(2, 2) :: origins_to_compare
  !   real(rk), dimension(2) :: radii_to_compare

  !   integer(ik) :: i
  !   real(rk) :: sound_speed
  !   logical :: is_inside
  !   real(rk), dimension(2) :: vel
  !   integer(ik) :: largest_cone_idx
  !   !< index of the cone with the largets radius

  !   real(rk), dimension(:) :: cone_radii
  !   !< (cell); radius of each cone based on the neighbor cell state

  !   real(rk), dimension(2, :) :: cone_origins
  !   !< ((x,y), cell); origin of each cone based on the neighbor cell state

  !   associate(u=>self%recon_state(2, :), &
  !             v=>self%recon_state(3, :), &
  !             cs=>self%sound_speed)

  !     ! Get the cone size based on each neighbor cell's state
  !     do i = 1, :
  !       vel = [u(i), v(i)]
  !       call get_cone_extents(tau=self%tau, &
  !                             xy=self%p_xy, vel=vel, sound_speed=cs(i), &
  !                             origin=cone_origins(:, i), radius=cone_radii(i))
  !     end do
  !   end associate

  !   ! Since the circle comparison is only 2 at a time, this
  !   ! does the diagonal cells first (1v3 and 2v4) then then
  !   ! compares the results from those for the final circle

  !   ! Compare cells 1 and 3
  !   origins_to_compare(:, 1) = cone_origins(:, 1)
  !   origins_to_compare(:, 2) = cone_origins(:, 3)
  !   radii_to_compare = [cone_radii(1), cone_radii(3)]
  !   call super_circle(origins=origins_to_compare, radii=radii_to_compare, &
  !                     new_origin=origin_1v3, new_radius=radius_1v3)

  !   ! Compare cells 2 and 4
  !   origins_to_compare(:, 1) = cone_origins(:, 2)
  !   origins_to_compare(:, 2) = cone_origins(:, 4)
  !   radii_to_compare = [cone_radii(2), cone_radii(4)]
  !   call super_circle(origins=origins_to_compare, radii=radii_to_compare, &
  !                     new_origin=origin_2v4, new_radius=radius_2v4)

  !   ! Compare result from 1 vs 3 and 2 vs 4
  !   origins_to_compare(:, 1) = origin_1v3
  !   origins_to_compare(:, 2) = origin_2v4
  !   radii_to_compare = [radius_1v3, radius_2v4]
  !   call super_circle(origins=origins_to_compare, radii=radii_to_compare, &
  !                     new_origin=origin, new_radius=radius)
  ! end subroutine get_transonic_cone_extents

  subroutine sanity_checks(self)
    !< Do some sanity checks to make sure the mach cone is valid
    class(mach_cone_collection_t), intent(in) :: self
    real(rk), dimension(self%ni, self%nj) :: arc_sum
    integer(ik) :: i, j

    arc_sum = sum(self%dtheta, dim=1)

    do j = 1, self%nj
      do i = 1, self%ni
        if(count(self%n_arcs_per_cell(:, i, j) >= 2) > 1) then
          error stop "Too many arcs in the mach cone (count(cone%n_arcs >= 2) > 1)"
        end if
      end do
    end do

    do j = 1, self%nj
      do i = 1, self%ni
        if(.not. equal(arc_sum(i, j), 2 * pi, 1e-6_rk)) then
          error stop "Cone arcs do not add up to 2pi"
        end if
      end do
    end do

  end subroutine sanity_checks

  ! subroutine write_cone(self, unit, iotype, v_list, iostat, iomsg)
  !   !< Implementation of `write(*,*) mach_cone_base_t`

  !   class(corner_mach_cone_t), intent(in) :: self     !< cone class
  !   integer, intent(in) :: unit           !< input/output unit
  !   character(*), intent(in) :: iotype    !< LISTDIRECTED or DTxxx
  !   integer, intent(in) :: v_list(:)      !< parameters from fmt spec.
  !   integer, intent(out) :: iostat        !< non zero on error, etc.
  !   character(*), intent(inout) :: iomsg  !< define if iostat non zero.

  !   integer(ik) :: i

  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) 'Mach Cone Details'//new_line('a')
  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) '================='//new_line('a')

  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) "location = '"//trim(self%cone_location)//"'"//new_line('a')
  !   write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "tau = ", self%tau, new_line('a')
  !   write(unit, '(a, es14.6, a)', iostat=iostat, iomsg=iomsg) "radius = ", self%radius, new_line('a')
  !   write(unit, '(a, es14.6, ",", es14.6, a)', iostat=iostat, iomsg=iomsg) 'cone["P (x,y)"] = [', self%p_xy, "]"//new_line('a')
  !   write(unit, '(a, es14.6, ",", es14.6, a)', iostat=iostat, iomsg=iomsg) &
  !     'cone["P'//"'(x,y)"//'"] = [', self%p_prime_xy, "]"//new_line('a')
  !   write(unit, '(a, es14.6, ",", es14.6, a)', iostat=iostat, iomsg=iomsg) &
  !     "dist = [", self%p_prime_xy - self%p_xy, "]"//new_line('a')

  !   do i = 1, :
  !     write(unit, '(a, i7,",",i7, a, i0, a)', iostat=iostat, iomsg=iomsg) &
  !       'cone["P'//"'(i,j)"//'"] = [', self%p_prime_ij(:, i), "] # Cell: ", i, new_line('a')
  !   end do

  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) "# Edge Vectors: tail (x,y) -> head (x,y)"//new_line('a')

  !   do i = 1, :
  !     write(unit, '(a, i0, a, 2(es10.3, ", ", es10.3, a))', iostat=iostat, iomsg=iomsg) "vector_", i, " = [[", &
  !       self%edge_vectors(:, 1, i), "], [", self%edge_vectors(:, 2, i), "] ]"//new_line('a')
  !   end do

  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
  !   write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "reference_density =     ", self%reference_density, new_line('a')
  !   write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "reference_x_velocity =  ", self%reference_u, new_line('a')
  !   write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "reference_y_velocity =  ", self%reference_v, new_line('a')
  !   write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "reference_sound_speed = ", self%reference_sound_speed, new_line('a')
  !   ! write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "reference_mach_number = ", self%reference_mach_number, new_line('a')

  !   write(unit, '(a)', iostat=iostat) new_line('a')
  !   write(unit, '(a, 4(i0, 1x),a)', iostat=iostat) '# Valid arcs in each cell: [', self%n_arcs_per_cell, ']'
  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
  !   ! write(unit, '(a)', iostat=iostat, iomsg=iomsg) '# Angles:       Theta_ib [deg]         Theta_ie [deg]'//new_line('a')
  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) '# Delta Theta (for each arc) [deg]'//new_line('a')
  !   do i = 1, :
  !     write(unit, '(a, i0,a,f9.5,",",f9.5, a)', iostat=iostat, iomsg=iomsg) &
  !       'dtheta[', i, '] = [', rad2deg(self%dtheta(:, i)), ']'//new_line('a')
  !   end do

  !   do i = 1, :
  !     write(unit, '(a, i0,a,f9.5,",",f9.5, a)', iostat=iostat, iomsg=iomsg) &
  !       'sin_dtheta[', i, '] = [', self%sin_dtheta(:, i), ']'//new_line('a')
  !   end do
  !   do i = 1, :
  !     write(unit, '(a, i0,a,f9.5,",",f9.5, a)', iostat=iostat, iomsg=iomsg) &
  !       'cos_dtheta[', i, '] = [', self%cos_dtheta(:, i), ']'//new_line('a')
  !   end do
  !   do i = 1, :
  !     write(unit, '(a, i0,a,f9.5,",",f9.5, a)', iostat=iostat, iomsg=iomsg) &
  !       'sin_d2theta[', i, '] = [', self%sin_d2theta(:, i), ']'//new_line('a')
  !   end do
  !   do i = 1, :
  !     write(unit, '(a, i0,a,f9.5,",",f9.5, a)', iostat=iostat, iomsg=iomsg) &
  !       'cos_d2theta[', i, '] = [', self%cos_d2theta(:, i), ']'//new_line('a')
  !   end do

  !   write(unit, '(a, f0.2, a)', iostat=iostat, iomsg=iomsg) 'arc_length_sum = ', &
  !     rad2deg(sum(self%dtheta)), ' '//new_line('a')

  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) 'Cell Primitive Vars State [rho,u,v,p]'//new_line('a')
  !   do i = 1, :
  !     write(unit, '(a, i0, a, 4(es10.3, 1x), a)', iostat=iostat, iomsg=iomsg) &
  !       'Cell: ', i, ' [ ', self%arc_primitive_vars(:, 1, i), '] (arc 1)'//new_line('a')
  !     write(unit, '(a, 4(es10.3, 1x), a)', iostat=iostat, iomsg=iomsg) &
  !       '        [ ', self%arc_primitive_vars(:, 2, i), '] (arc 2)'//new_line('a')
  !   end do

  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) "# P' in cell?"//new_line('a')
  !   do i = 1, :
  !     write(unit, '(a, i0, a, l2, a)', iostat=iostat, iomsg=iomsg) 'Cell: ', i, ' [', self%p_prime_in_cell(i), ' ]'//new_line('a')
  !   end do

  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')

  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) '# Neighbor cells contributing to the mach cone '//new_line('a')
  !   write(unit, '(a, i0, a)', iostat=iostat, iomsg=iomsg) 'neighbor_cells = ', self%n_neighbor_cells, new_line('a')
  !   do i = 1, :
  !     write(unit, '(a, i0, " = [", 3(es10.3, ","), es10.3,"]", a)', iostat=iostat, iomsg=iomsg) &
  !       'cell_', i, self%recon_state(:, i), new_line('a')
  !   end do

  !   write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
  ! end subroutine write_cone

end module mod_mach_cone_collection
