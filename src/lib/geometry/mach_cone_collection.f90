module mod_mach_cone_collection
#ifdef USE_OPENMP
  use omp_lib
#endif /* USE_OPENMP */
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use, intrinsic :: ieee_arithmetic
  use math_constants, only: pi, rad2deg
  use mod_globals, only: grid_is_orthogonal
  use mod_floating_point_utils, only: near_zero, equal
  use mod_eos, only: eos
  use mod_vector, only: vector_t
  use mod_geometry, only: find_arc_segments

  implicit none

  private
  public :: mach_cone_collection_t

  real(rk), parameter :: TINY_DIST_RATIO = 1.0e-8_rk

  type :: mach_cone_collection_t
    real(rk), dimension(:, :, :), allocatable :: recon_rho
    !< ((cell 1:N), i, j); Reconstructed value of density at P for each neighbor cell

    real(rk), dimension(:, :, :), allocatable :: recon_u
    !< ((cell 1:N), i, j); Reconstructed value of x-velocity at P for each neighbor cell

    real(rk), dimension(:, :, :), allocatable :: recon_v
    !< ((cell 1:N), i, j); Reconstructed value of y-velocity at P for each neighbor cell

    real(rk), dimension(:, :, :), allocatable :: recon_p
    !< ((cell 1:N), i, j); Reconstructed value of pressure at P for each neighbor cell

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
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); sin(theta_ie) - sin(theta_ib)

    real(rk), dimension(:, :, :), allocatable :: cos_dtheta
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); cos(theta_ie) - cos(theta_ib)

    real(rk), dimension(:, :, :), allocatable :: sin_d2theta
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); sin(2*theta_ie) - sine(2*theta_ib)

    real(rk), dimension(:, :, :), allocatable :: cos_d2theta
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); cos(2*theta_ie) - cos(2*theta_ib)

    integer(ik), dimension(:, :, :), allocatable :: p_prime_ij
    !< ((P' i, P' j), i, j) i,j location of P' or the apex of the Mach cone (global, not relative to P0)

    real(rk), dimension(:, :), allocatable :: p_prime_x  !< (i, j); x location of P'
    real(rk), dimension(:, :), allocatable :: p_prime_y  !< (i, j); y location of P'
    real(rk), dimension(:, :), allocatable :: p0_x       !< (i, j); x location of P0
    real(rk), dimension(:, :), allocatable :: p0_y       !< (i, j); y location of P0

    integer(ik), dimension(:, :, :), allocatable :: n_arcs_per_cell
    !< (cell 1:N); Number of arcs that each cell contributes (0, 1, or 2)

    logical, dimension(:, :, :), allocatable :: p_prime_in_cell
    !< ((cell 1:N), i, j); is P' inside the control volume?

    real(rk), dimension(:, :), allocatable :: pressure_p_prime
    !< (i, j); pressure from the P' cell (due to FVLEG assumptions, this is the reconstructed value at P0 of the cell that contains P')

    real(rk), dimension(:, :), allocatable :: density_p_prime
    !< density from the P' cell (due to FVLEG assumptions, this is the reconstructed value at P0 of the cell that contains P')

    real(rk), dimension(:, :, :), allocatable :: u
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); u velocity of the cell along the are; Note the dimension (cell 1:N)*2 , where 2 means # arcs
    real(rk), dimension(:, :, :), allocatable :: v
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); v velocity of the cell along the are; Note the dimension (cell 1:N)*2 , where 2 means # arcs
    real(rk), dimension(:, :, :), allocatable :: p
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); pressure of the cell along the are; Note the dimension (cell 1:N)*2 , where 2 means # arcs

    real(rk), dimension(:, :, :, :), allocatable:: edge_vectors
    !< ((x,y), (origin, cell 1:N) , i, j); vectors that define the cell edges

    integer(ik), dimension(:, :, :, :), allocatable :: cell_indices
    !< ((i,j), cell 1:N, i, j); indices of the neighboring cells

    integer(ik) :: n_neighbor_cells    !< Number of neighbor cells
    real(rk) :: tau                     !< Time evolution increment
    real(rk), dimension(:, :), allocatable :: radius                  !< Radius of the cone
    real(rk), dimension(:, :), allocatable :: reference_density       !< Reference density (e.g. neighbor averaged)
    real(rk), dimension(:, :), allocatable :: reference_u             !< Reference x velocity (e.g. neighbor averaged)
    real(rk), dimension(:, :), allocatable :: reference_v             !< Reference y velocity (e.g. neighbor averaged)
    real(rk), dimension(:, :), allocatable :: reference_sound_speed   !< Reference sound speed (e.g. neighbor averaged)
    ! logical, dimension(:, :), allocatable:: cone_is_transonic        !< Flag to enable special treatment for transonic cones
    ! logical, dimension(:, :), allocatable:: cone_is_centered         !< is the cone located at P0?
    character(len=32) :: cone_location  !< Corner or midpoint cone?
    integer(ik) :: ni !< # of cones in the i direction
    integer(ik) :: nj !< # of cones in the j direction
  contains
    procedure, public :: initialize
    procedure, public :: print
    procedure, private :: get_reference_state
    procedure, private :: sanity_checks
    procedure, private :: find_arc_angles
    procedure, private :: compute_trig_angles
  end type mach_cone_collection_t

contains

  subroutine initialize(self, tau, edge_vectors, reconstructed_state, cell_indices, cone_location)
    !< Class constructor

    class(mach_cone_collection_t), intent(inout) :: self

    character(len=*), intent(in) :: cone_location !< corner, down/up midpoint, left/right midpoint

    real(rk), intent(in) :: tau

    real(rk), dimension(:, 0:, :, :), intent(in) :: edge_vectors
    !< ((x,y), (origin, vector_1:N), i, j) ; set of vectors that define the corner

    integer(ik), dimension(:, :, :, :), intent(in) :: cell_indices
    !< ((i,j), (cell_1:N), i, j); set of indices for the neighboring cells -> needed to find P' i,j index

    real(rk), dimension(:, :, :, :), intent(in) :: reconstructed_state
    !< ((rho,u,v,p), (cell_1:N), i, j); reconstructed state for point P.

    integer(ik) :: idx, arc, c, i, j, ni, nj
    real(rk) :: recon_u, recon_v, recon_p

    self%ni = size(reconstructed_state, dim=3)
    self%nj = size(reconstructed_state, dim=4)
    self%cone_location = trim(cone_location)
    self%tau = tau

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
            write(std_err, '(a, 4(es16.8,1x))') 'Reconstructed density states [' // trim(self%cone_location) //'] (cell 1:N)', reconstructed_state(1, :, i, j)
            write(std_err, '(a, i0, 1x, i0, a)') 'Cone index [i, j]: ', i, j
            write(std_err, '(a, 8(i0, 1x))') 'Cone neighbor cell indices [i, j]: ', cell_indices(:, :, i, j)
            error stop "Error in mach_cone_collection_t%initialize(), density in the reconstructed state is < 0"
          end if

          if(reconstructed_state(4, c, i, j) < 0.0_rk) then
            write(std_err, '(a, 4(es16.8,1x))') 'Reconstructed pressure states [' // trim(self%cone_location) //'] (cell 1:N): ', reconstructed_state(4, :, i, j)
            write(std_err, '(a, i0, 1x, i0, a)') 'Cone index [i, j]: [', i, j, ']'
            write(std_err, '(a, 8(i0, 1x))') 'Cone neighbor cell indices [i, j]: ', cell_indices(:, :, i, j)
            error stop "Error in mach_cone_collection_t%initialize(), pressure in the reconstructed state is < 0"
          end if
        end do
      end do
    end do

    associate(nc=>self%n_neighbor_cells, ni=>self%ni, nj=>self%nj)
      ! allocate(self%cell_is_supersonic(ni, nj))
      ! self%cell_is_supersonic = .false.

      if(.not. allocated(self%edge_vectors)) then
        allocate(self%edge_vectors(2, 0:nc, ni, nj))
        self%edge_vectors = edge_vectors
      end if

      if(.not. allocated(self%p_prime_in_cell)) allocate(self%p_prime_in_cell(nc, ni, nj))
      self%p_prime_in_cell = .false.

      if(.not. allocated(self%recon_rho)) allocate(self%recon_rho(nc, ni, nj))
      if(.not. allocated(self%recon_u)) allocate(self%recon_u(nc, ni, nj))
      if(.not. allocated(self%recon_v)) allocate(self%recon_v(nc, ni, nj))
      if(.not. allocated(self%recon_p)) allocate(self%recon_p(nc, ni, nj))
      self%recon_rho = reconstructed_state(1, :, :, :)
      self%recon_u = reconstructed_state(2, :, :, :)
      self%recon_v = reconstructed_state(3, :, :, :)
      self%recon_p = reconstructed_state(4, :, :, :)

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

      if(.not. allocated(self%p_prime_x)) allocate(self%p_prime_x(ni, nj))
      self%p_prime_x = 0.0_rk

      if(.not. allocated(self%p_prime_y)) allocate(self%p_prime_y(ni, nj))
      self%p_prime_y = 0.0_rk

      if(.not. allocated(self%p0_x)) allocate(self%p0_x(ni, nj))
      self%p0_x = 0.0_rk

      if(.not. allocated(self%p0_y)) allocate(self%p0_y(ni, nj))
      self%p0_y = 0.0_rk

      if(.not. allocated(self%pressure_p_prime)) allocate(self%pressure_p_prime(ni, nj))
      self%pressure_p_prime = 0.0_rk

      if(.not. allocated(self%density_p_prime)) allocate(self%density_p_prime(ni, nj))
      self%density_p_prime = 0.0_rk

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
      if(.not. allocated(self%reference_density)) allocate(self%reference_density(ni, nj))
      if(.not. allocated(self%reference_u)) allocate(self%reference_u(ni, nj))
      if(.not. allocated(self%reference_v)) allocate(self%reference_v(ni, nj))
      if(.not. allocated(self%reference_sound_speed)) allocate(self%reference_sound_speed(ni, nj))

      ! if(.not. allocated(self%cone_is_transonic)) allocate(self%cone_is_transonic(ni, nj))
      ! self%cone_is_transonic = .false.

      ! if(.not. allocated(self%cone_is_centered)) allocate(self%cone_is_centered(ni, nj))
      ! self%cone_is_centered = .false.
    end associate

    call self%get_reference_state()

    ! Set the P' vector. Remember that all edge vector sets are moved to the origin, so P0 is at (0,0)
    do j = 1, self%nj
      do i = 1, self%ni
        self%p0_x(i, j) = self%edge_vectors(1, 0, i, j)
        self%p0_y(i, j) = self%edge_vectors(2, 0, i, j)
        associate(x=>self%p0_x(i, j), y=>self%p0_y(i, j), &
                  u_tilde=>self%reference_u(i, j), v_tilde=>self%reference_v(i, j))
          self%p_prime_x(i, j) = x - self%tau * u_tilde
          self%p_prime_y(i, j) = y - self%tau * v_tilde
        end associate
      end do
    end do

    self%radius = self%tau * self%reference_sound_speed
    ! where(abs(self%p0_x - self%p_prime_x) / self%radius < TINY_DIST_RATIO) self%p_prime_x = self%p0_x
    ! where(abs(self%p0_y - self%p_prime_y) / self%radius < TINY_DIST_RATIO) self%p_prime_y = self%p0_y

    call self%compute_trig_angles()

    ! Assign the reconstructed quantities to the primitive var values for each arc of the mach cone
    !$omp parallel default(shared) private(i,j,c, arc, idx,recon_u, recon_v, recon_p)
    !$omp do
    do j = 1, self%nj
      do i = 1, self%ni
        do c = 1, self%n_neighbor_cells
          recon_u = self%recon_u(c, i, j)
          recon_v = self%recon_v(c, i, j)
          recon_p = self%recon_p(c, i, j)
          do arc = 1, 2
            idx = arc + (c - 1) * 2 ! convert indexing for better storage
            self%u(idx, i, j) = recon_u
            self%v(idx, i, j) = recon_v
            self%p(idx, i, j) = recon_p
          end do
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do j = 1, self%nj
      do i = 1, self%ni
        do c = 1, self%n_neighbor_cells
          if(self%p_prime_in_cell(c, i, j)) then
            self%density_p_prime(i, j) = self%recon_rho(c, i, j)
            self%pressure_p_prime(i, j) = self%recon_p(c, i, j)
          end if
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

    ! call self%determine_p_prime_cell()
    call self%sanity_checks()
  end subroutine initialize

  subroutine get_reference_state(self)
    !< Find the reference state of the cone only based on the cells that are touched by the cone
    class(mach_cone_collection_t), intent(inout) :: self

    real(rk), dimension(:, :, :), allocatable :: mach_number, sound_speed
    integer(ik), dimension(:, :), allocatable :: n_supersonic_cells
    integer(ik) :: i, j, c
    real(rk) :: gamma, mach_u, mach_v, n_cells_real

    gamma = eos%get_gamma()

    n_cells_real = real(self%n_neighbor_cells, rk)

    !$omp parallel default(shared) private(i,j,c)
    !$omp do
    do j = 1, self%nj
      do i = 1, self%ni
        self%reference_density(i, j) = sum(self%recon_rho(:, i, j)) / n_cells_real
        self%reference_u(i, j) = sum(self%recon_u(:, i, j)) / n_cells_real
        self%reference_v(i, j) = sum(self%recon_v(:, i, j)) / n_cells_real
        self%reference_sound_speed(i, j) = sqrt(gamma * (sum(self%recon_p(:, i, j)) / n_cells_real) / self%reference_density(i, j))
      end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine get_reference_state

  subroutine compute_trig_angles(self)
    !< Precompute some of the trig values, like sin(theta_ie), etc. to save compute time. If
    !< the grid is orthogonal and the references state is at rest, then the math
    !< gets even easier

    class(mach_cone_collection_t), intent(inout) :: self

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
    integer(ik) :: idx, idx_max
    !< index conversion to go from ((arc1, arc2), (cell 1:N_CELLS), i, j) to ((cell_1_arc_1, cell_N_arc_1:cell_1_arc_N, cell_N_arc_N), i, j)

    integer(ik), dimension(:, :, :), allocatable :: n_arcs_per_cell

    idx_max = 2 * self%n_neighbor_cells
    allocate(sin_theta_ib(idx_max, self%ni, self%nj))
    sin_theta_ib = 0.0_rk

    allocate(cos_theta_ib(idx_max, self%ni, self%nj))
    cos_theta_ib = 0.0_rk

    allocate(sin_theta_ie(idx_max, self%ni, self%nj))
    sin_theta_ie = 0.0_rk

    allocate(cos_theta_ie(idx_max, self%ni, self%nj))
    cos_theta_ie = 0.0_rk

    allocate(n_arcs_per_cell(self%n_neighbor_cells, self%ni, self%nj))
    n_arcs_per_cell = 0

    call self%find_arc_angles(theta_ib, theta_ie, n_arcs_per_cell)

    self%n_arcs_per_cell = n_arcs_per_cell

    !$omp parallel default(shared) private(i,j,c,arc,idx)
    !$omp do
    do j = 1, self%nj
      do i = 1, self%ni
        do c = 1, self%n_neighbor_cells
          do arc = 1, 2
            idx = arc + (c - 1) * 2
            self%dtheta(idx, i, j) = abs(theta_ie(idx, i, j) - theta_ib(idx, i, j))
            sin_theta_ib(idx, i, j) = sin(theta_ib(idx, i, j))
            cos_theta_ib(idx, i, j) = cos(theta_ib(idx, i, j))
            sin_theta_ie(idx, i, j) = sin(theta_ie(idx, i, j))
            cos_theta_ie(idx, i, j) = cos(theta_ie(idx, i, j))
          end do
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do j = 1, self%nj
      do i = 1, self%ni
        do idx = 1, idx_max
          if(abs(sin_theta_ib(idx, i, j)) < 1e-15_rk) sin_theta_ib(idx, i, j) = 0.0_rk
          if(abs(cos_theta_ib(idx, i, j)) < 1e-15_rk) cos_theta_ib(idx, i, j) = 0.0_rk
          if(abs(sin_theta_ie(idx, i, j)) < 1e-15_rk) sin_theta_ie(idx, i, j) = 0.0_rk
          if(abs(cos_theta_ie(idx, i, j)) < 1e-15_rk) cos_theta_ie(idx, i, j) = 0.0_rk
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do j = 1, self%nj
      do i = 1, self%ni
        do idx = 1, idx_max
          self%sin_dtheta(idx, i, j) = sin_theta_ie(idx, i, j) - sin_theta_ib(idx, i, j)
          self%cos_dtheta(idx, i, j) = cos_theta_ie(idx, i, j) - cos_theta_ib(idx, i, j)
          self%sin_d2theta(idx, i, j) = 2.0_rk * sin_theta_ie(idx, i, j) * cos_theta_ie(idx, i, j) - &
                                        2.0_rk * sin_theta_ib(idx, i, j) * cos_theta_ib(idx, i, j)
          self%cos_d2theta(idx, i, j) = (2.0_rk * cos_theta_ie(idx, i, j)**2 - 1.0_rk) - &
                                        (2.0_rk * cos_theta_ib(idx, i, j)**2 - 1.0_rk)
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

    deallocate(sin_theta_ib)
    deallocate(cos_theta_ib)
    deallocate(sin_theta_ie)
    deallocate(cos_theta_ie)
    deallocate(n_arcs_per_cell)

  end subroutine compute_trig_angles

  subroutine find_arc_angles(self, theta_ib, theta_ie, n_arcs_per_cell)
    !< Find the starting and ending angle of each arc that is held within each cell

    class(mach_cone_collection_t), intent(inout) :: self

    real(rk), dimension(:, :, :), allocatable, intent(out) :: theta_ib
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); starting angle [rad] for the arc contained in each cell

    real(rk), dimension(:, :, :), allocatable, intent(out):: theta_ie
    !< ((cell_1_arc_1:cell_N_arc_2), i, j); ending angle [rad] for the arc contained in each cell

    integer(ik), dimension(:, :, :), allocatable, intent(out) :: n_arcs_per_cell !< # of arcs in a given cell
    ! logical, dimension(:, :, :), intent(out) :: p_prime_in_cell
    ! integer(ik), dimension(:, :, :), intent(out) :: p_prime_ij

    real(rk), dimension(2, 2) :: theta_ib_ie
    !< ((start,end), (arc1, arc2)); start/end angle for a given cell

    real(rk), dimension(2) :: second_vector, first_vector
    real(rk), dimension(2, 4) :: corner_vector_set !< ((x,y), vector 1:4)
    real(rk), dimension(2, 2) :: midpoint_vector_set !< ((x,y), vector 1:2)
    real(rk), dimension(2) :: p_prime, p_0
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

    !$omp parallel default(shared) &
    !$omp private(corner_vector_set, midpoint_vector_set, first_vector, second_vector) &
    !$omp private(i, j, arc, c, idx) &
    !$omp private(p_prime, p_0, p_prime_cross_vec_1, vec_2_cross_p_prime, n_arcs, theta_ib_ie)
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

      !$omp do
      do j = 1, self%nj
        do i = 1, self%ni
          corner_vector_set = self%edge_vectors(:, 1:4, i, j)
          p_prime = [self%p_prime_x(i, j), self%p_prime_y(i, j)]
          p_0 = [self%p0_x(i, j), self%p0_y(i, j)]

          do c = 1, self%n_neighbor_cells
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

            associate(x0=>self%p0_x(i, j), y0=>self%p0_y(i, j))
              ! Determine if P' is in this cell
              p_prime_cross_vec_1 = ((self%p_prime_x(i, j) - x0) * (first_vector(2) - y0)) - &
                                    ((self%p_prime_y(i, j) - y0) * (first_vector(1) - x0))

              vec_2_cross_p_prime = ((second_vector(1) - x0) * (self%p_prime_y(i, j) - y0)) - &
                                    ((second_vector(2) - y0) * (self%p_prime_x(i, j) - x0))
            end associate

            if(p_prime_cross_vec_1 <= 0.0_rk .and. vec_2_cross_p_prime <= 0.0_rk) then
              self%p_prime_in_cell(c, i, j) = .true.
              self%p_prime_ij(:, i, j) = [i, j]
            end if

            call find_arc_segments(origin=p_0, vector_1=first_vector, vector_2=second_vector, &
                                   origin_in_cell=self%p_prime_in_cell(c, i, j), &
                                   circle_xy=p_prime, &
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
      !$omp end do

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

      !$omp do
      do j = 1, self%nj
        do i = 1, self%ni
          midpoint_vector_set = self%edge_vectors(:, 1:2, i, j)
          p_prime = [self%p_prime_x(i, j), self%p_prime_y(i, j)]
          p_0 = [self%p0_x(i, j), self%p0_y(i, j)]

          do c = 1, 2
            select case(c)
            case(1) ! Cell 1
              first_vector = midpoint_vector_set(:, 2)
              second_vector = midpoint_vector_set(:, 1)
            case(2) ! Cell 2
              first_vector = midpoint_vector_set(:, 1)
              second_vector = midpoint_vector_set(:, 2)
            end select

            associate(x0=>self%p0_x(i, j), y0=>self%p0_y(i, j))
              ! Determine if P' is in this cell
              p_prime_cross_vec_1 = ((self%p_prime_x(i, j) - x0) * (first_vector(2) - y0)) - &
                                    ((self%p_prime_y(i, j) - y0) * (first_vector(1) - x0))

              vec_2_cross_p_prime = ((second_vector(1) - x0) * (self%p_prime_y(i, j) - y0)) - &
                                    ((second_vector(2) - y0) * (self%p_prime_x(i, j) - x0))
            end associate

            if(p_prime_cross_vec_1 <= 0.0_rk .and. vec_2_cross_p_prime <= 0.0_rk) then
              self%p_prime_in_cell(c, i, j) = .true.
              self%p_prime_ij(:, i, j) = [i, j]
            end if

            call find_arc_segments(origin=p_0, vector_1=first_vector, vector_2=second_vector, &
                                   origin_in_cell=self%p_prime_in_cell(c, i, j), &
                                   circle_xy=p_prime, &
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
      !$omp end do

    case default
      error stop 'Error in mach_cone_collection_t%find_arc_angles(), invalid cone location'
    end select
    !$omp end parallel

  end subroutine find_arc_angles

  subroutine sanity_checks(self)
    !< Do some sanity checks to make sure the mach cone is valid
    class(mach_cone_collection_t), intent(in) :: self
    real(rk), dimension(self%ni, self%nj) :: arc_sum
    integer(ik) :: i, j

    arc_sum = sum(self%dtheta, dim=1)

    do j = 1, self%nj
      do i = 1, self%ni
        if(count(self%n_arcs_per_cell(:, i, j) >= 2) > 1) then
          call self%print(i, j)
          error stop "Too many arcs in the mach cone (count(cone%n_arcs >= 2) > 1)"
        end if
      end do
    end do

    do j = 1, self%nj
      do i = 1, self%ni
        if(.not. equal(arc_sum(i, j), 2 * pi, 1e-6_rk)) then
          print *, 'arc_sum(i, j)', arc_sum(i, j), 2 * pi
          call self%print(i, j)
          error stop "Cone arcs do not add up to 2pi"
        end if
      end do
    end do

  end subroutine sanity_checks

  subroutine print(self, i, j)
    class(mach_cone_collection_t), intent(in) :: self
    integer(ik), intent(in) :: i, j
    integer(ik) :: c, arc, idx
    real(rk) :: vlen, dist

    write(*, '(a)') 'Mach Cone Details'
    write(*, '(a)') '================='

    write(*, '(a)') "location = '"//trim(self%cone_location)//"'"
    write(*, '(a, es16.8)') "tau = ", self%tau
    write(*, '(a, es16.8)') "radius = ", self%radius(i, j)
    write(*, '(a, es16.8, ",", es16.8, a)') 'cone["P (x,y)"] = [', self%p0_x(i, j), self%p0_y(i, j), "]"
    write(*, '(a, es16.8, ",", es16.8, a)') 'cone["P'//"'(x,y)"//'"] = [', self%p_prime_x(i, j), self%p_prime_y(i, j), "]"
    write(*, '(a, es16.8, ",", es16.8, a)') "dist =            [", self%p_prime_x(i,j) - self%p0_x(i,j), self%p_prime_y(i,j) - self%p0_y(i,j), "]"

    write(*, '(a, i7,",",i7, a)') 'cone["P'//"'(i,j)"//'"] = [', self%p_prime_ij(:, i, j), "]"

    print *
    write(*, *) 'Neighbor cell indices'
    do c = 1, self%n_neighbor_cells
      write(*, '(a,i0,a, 2i5)') "cell_", c, ' = ', self%cell_indices(:, c, i, j)
    end do

    print *
    write(*, '(a)') "# Edge Vectors: tail (x,y) -> head (x,y)"
    do c = 1, self%n_neighbor_cells
      write(*, '(a,i0,a, 2(es16.8, a))') "vector_", c, " = [", &
        self%edge_vectors(1, c, i, j), ",", self%edge_vectors(2, c, i, j), "]"
    end do

    do c = 1, self%n_neighbor_cells
      vlen = norm2(self%edge_vectors(:, c, i, j) - self%edge_vectors(:, 0, i, j))
      write(*, '(a, i0,a,es16.8)') "vector_", c, "_len = ", vlen
    end do

    print *
    write(*, '(a, es16.8)') "reference_density =     ", self%reference_density(i, j)
    write(*, '(a, es16.8)') "reference_x_velocity =  ", self%reference_u(i, j)
    write(*, '(a, es16.8)') "reference_y_velocity =  ", self%reference_v(i, j)
    write(*, '(a, es16.8)') "reference_sound_speed = ", self%reference_sound_speed(i, j)

    write(*, *) '# Valid arcs in each cell: [', self%n_arcs_per_cell(:, i, j), ']'
    print *

    write(*, '(a)') '# Delta Theta (for each arc) [deg]'
    do c = 1, self%n_neighbor_cells
      do arc = 1, 2
        idx = arc + (c - 1) * 2
        write(*, '(a, i0, a, f9.4, 2(a, i0))') &
          'dtheta[', c, '] = [', rad2deg(self%dtheta(idx, i, j)), '] arc ', arc, ', cell ', c
      end do
    end do
    print *

    do c = 1, self%n_neighbor_cells
      do arc = 1, 2
        idx = arc + (c - 1) * 2
        write(*, '(a, i0, a, es16.8, 2(a, i0))') &
          'sin_dtheta[', c, '] = [', self%sin_dtheta(idx, i, j), '] arc ', arc, ', cell ', c
      end do
    end do
    print *

    do c = 1, self%n_neighbor_cells
      do arc = 1, 2
        idx = arc + (c - 1) * 2
        write(*, '(a, i0, a, es16.8, 2(a, i0))') &
          'cos_dtheta[', c, '] = [', self%cos_dtheta(idx, i, j), '] arc ', arc, ', cell ', c
      end do
    end do
    print *

    do c = 1, self%n_neighbor_cells
      do arc = 1, 2
        idx = arc + (c - 1) * 2
        write(*, '(a, i0, a, es16.8, 2(a, i0))') &
          'sin_d2theta[', c, '] = [', self%sin_d2theta(idx, i, j), '] arc ', arc, ', cell ', c
      end do
    end do
    print *

    do c = 1, self%n_neighbor_cells
      do arc = 1, 2
        idx = arc + (c - 1) * 2
        write(*, '(a, i0, a, es16.8, 2(a, i0))') &
          'cos_d2theta[', c, '] = [', self%cos_d2theta(idx, i, j), '] arc ', arc, ', cell ', c
      end do
    end do
    print *

    write(*, '(a, es16.8)') 'arc_length_sum = ', rad2deg(sum(self%dtheta(:, i, j))), ' '

    write(*, '(a)') 'Neighbor state [rho,u,v,p]'
    do c = 1, self%n_neighbor_cells
      write(*, '(a, i0, a, 4(es16.8), a)') 'Cell: ', c, ' [ ', self%recon_rho(c, i, j), &
        self%recon_u(c, i, j), &
        self%recon_v(c, i, j), &
        self%recon_p(c, i, j), ']'
    end do
    print *

    write(*, '(a)') 'Cone state [u,v,p]'
    do c = 1, self%n_neighbor_cells
      do arc = 1, 2
        idx = arc + (c - 1) * 2
        write(*, '(a, 3(es16.8), 2(a, i0))') ' [', self%u(idx, i, j), &
          self%v(idx, i, j), &
          self%p(idx, i, j), &
          '] arc ', arc, ', cell ', c
      end do
    end do
    print *

    ! write(*, '(a)') "# P' in cell?"
    ! do i = 1, N_CELLS
    !   write(*, '(a, i0, a, l2, a)') 'Cell: ', i, ' [', self%p_prime_in_cell(i), ' ]'
    ! end do

    ! write(*, '(a)') new_line('a')

    ! write(*, '(a)') '# Neighbor cells contributing to the mach cone '
    ! write(*, '(a, i0, a)') 'neighbor_cells = ', self%n_neighbor_cells
    ! do i = 1, N_CELLS
    !   write(*, '(a, i0, " = [", 3(es16.8, ","), es16.8,"]", a)') &
    !     'cell_', i, self%recon_state(:, i)
    ! end do

    ! write(*, '(a)') new_line('a')
  end subroutine
end module mod_mach_cone_collection
