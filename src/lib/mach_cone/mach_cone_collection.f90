module mod_mach_cone_collection
#ifdef USE_OPENMP
  use omp_lib
#endif /* USE_OPENMP */
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use, intrinsic :: ieee_arithmetic
  use math_constants, only: pi, rad2deg
  use mod_globals, only: grid_is_orthogonal
  use mod_floating_point_utils, only: near_zero, equal, neumaier_sum
  use mod_abstract_reconstruction, only: abstract_reconstruction_t
  use mod_eos, only: eos
  use mod_vector, only: vector_t
  use mod_geometry, only: find_arc_segments, super_circle

  implicit none

  private
  public :: mach_cone_collection_t

  logical, parameter :: filter_small_dist = .true.
  real(rk), parameter :: TINY_DIST_RATIO = 1.0e-8_rk

  type :: mach_cone_collection_t
    integer(ik) :: ilo = 0 !< low i boundary based on cell indexing
    integer(ik) :: ihi = 0 !< high i boundary based on cell indexing
    integer(ik) :: jlo = 0 !< low j boundary based on cell indexing
    integer(ik) :: jhi = 0 !< high j boundary based on cell indexing

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
    logical :: is_initialized = .false.

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
    real(rk), dimension(:, :), allocatable :: reference_pressure      !< Reference pressure (e.g. neighbor averaged)
    real(rk), dimension(:, :), allocatable :: reference_sound_speed   !< Reference sound speed (e.g. neighbor averaged)
    logical, dimension(:, :), allocatable:: cone_is_transonic        !< Flag to enable special treatment for transonic cones
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
    procedure, private :: get_p_prime_quantities
    procedure, private :: get_transonic_cone_extents
    final :: finalize
  end type mach_cone_collection_t

contains

  subroutine initialize(self, tau, edge_vectors, &
                        reconstructed_rho, reconstructed_u, reconstructed_v, reconstructed_p, cell_indices, cone_location, &
                        reconstruction_operator, lbounds)
    !< Class constructor

    class(mach_cone_collection_t), intent(inout) :: self
    integer(ik), dimension(3) :: lbounds
    class(abstract_reconstruction_t), intent(in) :: reconstruction_operator !< Recon operator used to find P' values
    character(len=*), intent(in) :: cone_location !< corner, down/up midpoint, left/right midpoint

    real(rk), intent(in) :: tau

    real(rk), dimension(:, 0:, lbounds(2):, lbounds(3):), contiguous, intent(in) :: edge_vectors
    !< ((x,y), (origin, vector_1:N), i, j) ; set of vectors that define the corner

    integer(ik), dimension(:, :, lbounds(2):, lbounds(3):), contiguous, intent(in) :: cell_indices
    !< ((i,j), (cell_1:N), i, j); set of indices for the neighboring cells -> needed to find P' i,j index

    real(rk), dimension(:, lbounds(2):, lbounds(3):), contiguous, intent(in) :: reconstructed_rho
    real(rk), dimension(:, lbounds(2):, lbounds(3):), contiguous, intent(in) :: reconstructed_u
    real(rk), dimension(:, lbounds(2):, lbounds(3):), contiguous, intent(in) :: reconstructed_v
    real(rk), dimension(:, lbounds(2):, lbounds(3):), contiguous, intent(in) :: reconstructed_p
    !< ((cell_1:N), i, j); reconstructed state for point P.

    integer(ik) :: idx, arc, c, i, j, ni, nj, nc, n_p_prime_cells

    real(rk) :: recon_u, recon_v, recon_p
    real(rk), dimension(2) :: transonic_origin !< (x,y); origin for the transonic cone
    real(rk) :: transonic_radius               !< radius for the transonic cone

    self%ni = size(reconstructed_rho, dim=2)
    self%nj = size(reconstructed_rho, dim=3)
    if(self%ni <= 0 .or. self%nj <= 0) then
      error stop "Error in mach_cone_collection_t%initialize(), ni <= 1 .or. nj <= 1, e.g. not enough mach cones!"
    end if

    self%cone_location = trim(cone_location)
    self%tau = tau
    if(self%tau < tiny(1.0_rk)) error stop "Error in mach_cone_collection_t%initialize(), tau < tiny(1.0)"

    ni = self%ni
    nj = self%nj

    self%ilo = lbound(reconstructed_rho, dim=2)
    self%ihi = ubound(reconstructed_rho, dim=2)
    self%jlo = lbound(reconstructed_rho, dim=3)
    self%jhi = ubound(reconstructed_rho, dim=3)

    select case(trim(cone_location))
    case('corner')
      self%n_neighbor_cells = 4
    case('down/up midpoint', 'left/right midpoint')
      self%n_neighbor_cells = 2
    case default
      error stop "Error in mach_cone_collection_t%initialize(), unsupported cone location"
    end select

    nc = self%n_neighbor_cells
    !$omp parallel default(none) &
    !$omp shared(self, reconstructed_rho, reconstructed_u, reconstructed_v, reconstructed_p) &
    !$omp shared(cell_indices, edge_vectors) &
    !$omp firstprivate(ni, nj, nc) &
    !$omp private(i,j,c)
    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        do c = 1, nc
          if(reconstructed_rho(c, i, j) < 0.0_rk) then
            write(std_err, '(a, 4(es16.8,1x))') 'Reconstructed density states ['// &
              trim(self%cone_location)//'] (cell 1:N)', reconstructed_rho(:, i, j)
            write(std_err, '(a, i0, 1x, i0, a)') 'Cone index [i, j]: ', i, j
            write(std_err, '(a, 8(i0, 1x))') 'Cone neighbor cell indices [i, j]: ', cell_indices(:, :, i, j)
            error stop "Error in mach_cone_collection_t%initialize(), density in the reconstructed state is < 0"
          end if

          if(reconstructed_p(c, i, j) < 0.0_rk) then
            write(std_err, '(a, 4(es16.8,1x))') 'Reconstructed pressure states ['// &
              trim(self%cone_location)//'] (cell 1:N): ', reconstructed_p(:, i, j)
            write(std_err, '(a, i0, ", ", i0, a)') 'Cone index [i, j]: [', i, j, ']'
            write(std_err, '(a, 4("[",i0, ", ", i0, "] "))') 'Cone neighbor cell indices [i, j]: ', cell_indices(:, :, i, j)
            error stop "Error in mach_cone_collection_t%initialize(), pressure in the reconstructed state is < 0"
          end if
        end do
      end do
    end do
    !$omp end do

    !$omp single
    if(.not. self%is_initialized) then
      associate(nc=>self%n_neighbor_cells, &
                ilo=>self%ilo, ihi=>self%ihi, jlo=>self%jlo, jhi=>self%jhi)

        if(.not. allocated(self%edge_vectors)) then
          allocate(self%edge_vectors(2, 0:nc, self%ilo:self%ihi, self%jlo:self%jhi))
          self%edge_vectors = edge_vectors
        end if

        if(.not. allocated(self%p_prime_in_cell)) allocate(self%p_prime_in_cell(nc, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%recon_rho)) allocate(self%recon_rho(nc, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%recon_u)) allocate(self%recon_u(nc, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%recon_v)) allocate(self%recon_v(nc, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%recon_p)) allocate(self%recon_p(nc, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%dtheta)) allocate(self%dtheta(nc * 2, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%sin_dtheta)) allocate(self%sin_dtheta(nc * 2, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%cos_dtheta)) allocate(self%cos_dtheta(nc * 2, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%sin_d2theta)) allocate(self%sin_d2theta(nc * 2, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%cos_d2theta)) allocate(self%cos_d2theta(nc * 2, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%p_prime_x)) allocate(self%p_prime_x(ilo:ihi, jlo:jhi))
        if(.not. allocated(self%p_prime_y)) allocate(self%p_prime_y(ilo:ihi, jlo:jhi))
        if(.not. allocated(self%p0_x)) allocate(self%p0_x(ilo:ihi, jlo:jhi))
        if(.not. allocated(self%p0_y)) allocate(self%p0_y(ilo:ihi, jlo:jhi))
        if(.not. allocated(self%pressure_p_prime)) allocate(self%pressure_p_prime(ilo:ihi, jlo:jhi))
        if(.not. allocated(self%density_p_prime)) allocate(self%density_p_prime(ilo:ihi, jlo:jhi))
        if(.not. allocated(self%u)) allocate(self%u(nc * 2, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%v)) allocate(self%v(nc * 2, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%p)) allocate(self%p(nc * 2, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%p_prime_ij)) allocate(self%p_prime_ij(2, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%n_arcs_per_cell)) allocate(self%n_arcs_per_cell(nc, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%cell_indices)) allocate(self%cell_indices(2, nc, ilo:ihi, jlo:jhi))
        if(.not. allocated(self%radius)) allocate(self%radius(ilo:ihi, jlo:jhi))
        if(.not. allocated(self%reference_density)) allocate(self%reference_density(ilo:ihi, jlo:jhi))
        if(.not. allocated(self%reference_u)) allocate(self%reference_u(ilo:ihi, jlo:jhi))
        if(.not. allocated(self%reference_v)) allocate(self%reference_v(ilo:ihi, jlo:jhi))
        if(.not. allocated(self%reference_sound_speed)) allocate(self%reference_sound_speed(ilo:ihi, jlo:jhi))
        if(.not. allocated(self%reference_pressure)) allocate(self%reference_pressure(ilo:ihi, jlo:jhi))
        if(.not. allocated(self%cone_is_transonic)) allocate(self%cone_is_transonic(ilo:ihi, jlo:jhi))
      end associate

    end if
    !$omp end single

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        do c = 1, nc
          self%cell_indices(:, c, i, j) = cell_indices(:, c, i, j)
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        do c = 1, nc
          self%recon_rho(c, i, j) = reconstructed_rho(c, i, j)
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        do c = 1, nc
          self%recon_u(c, i, j) = reconstructed_u(c, i, j)
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        do c = 1, nc
          self%recon_v(c, i, j) = reconstructed_v(c, i, j)
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        do c = 1, nc
          self%recon_p(c, i, j) = reconstructed_p(c, i, j)
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

    call self%get_reference_state()

    !$omp parallel default(none) &
    !$omp private(transonic_origin, transonic_radius), &
    !$omp shared(self) &
    !$omp private(i,j)
    if(.not. self%is_initialized) then
      !$omp do
      do j = self%jlo, self%jhi
        do i = self%ilo, self%ihi
          self%p0_x(i, j) = self%edge_vectors(1, 0, i, j)
          self%p0_y(i, j) = self%edge_vectors(2, 0, i, j)
        end do
      end do
      !$omp end do
    end if

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi

        ! if(self%cone_is_transonic(i, j)) then
        !   call self%get_transonic_cone_extents(i, j, transonic_origin, transonic_radius)
        !   self%p_prime_x(i, j) = transonic_origin(1)
        !   self%p_prime_y(i, j) = transonic_origin(2)
        !   self%radius(i, j) = transonic_radius
        ! else
        self%radius(i, j) = self%tau * self%reference_sound_speed(i, j)

        if(abs(self%p0_x(i, j) - self%tau * self%reference_u(i, j)) < epsilon(1.0_rk)) then
          self%p_prime_x(i, j) = self%p0_x(i, j)
        else
          self%p_prime_x(i, j) = self%p0_x(i, j) - self%tau * self%reference_u(i, j)
        end if

        if(abs(self%p0_y(i, j) - self%tau * self%reference_v(i, j)) < epsilon(1.0_rk)) then
          self%p_prime_y(i, j) = self%p0_y(i, j)
        else
          self%p_prime_y(i, j) = self%p0_y(i, j) - self%tau * self%reference_v(i, j)
        end if

        ! end if
      end do
    end do
    !$omp end do
    !$omp end parallel

    call self%compute_trig_angles()

    ! Assign the reconstructed quantities to the primitive var values for each arc of the mach cone
    !$omp parallel default(none) &
    !$omp private(n_p_prime_cells) &
    !$omp shared(self) &
    !$omp private(i,j,c, arc, idx,recon_u, recon_v, recon_p)

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
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
    !$omp end do nowait

    !$omp end parallel
    call self%get_p_prime_quantities(recon_operator=reconstruction_operator)

    call self%sanity_checks()
    if(.not. self%is_initialized) self%is_initialized = .true.
  end subroutine initialize

  subroutine finalize(self)
    ! Cleanup all the allocated arrays upon destruction of the mach_cone_collection_t object
    type(mach_cone_collection_t), intent(inout) :: self

    if(allocated(self%recon_rho)) deallocate(self%recon_rho)
    if(allocated(self%recon_u)) deallocate(self%recon_u)
    if(allocated(self%recon_v)) deallocate(self%recon_v)
    if(allocated(self%recon_p)) deallocate(self%recon_p)
    if(allocated(self%dtheta)) deallocate(self%dtheta)
    if(allocated(self%sin_dtheta)) deallocate(self%sin_dtheta)
    if(allocated(self%cos_dtheta)) deallocate(self%cos_dtheta)
    if(allocated(self%sin_d2theta)) deallocate(self%sin_d2theta)
    if(allocated(self%cos_d2theta)) deallocate(self%cos_d2theta)
    if(allocated(self%p_prime_ij)) deallocate(self%p_prime_ij)
    if(allocated(self%p_prime_x)) deallocate(self%p_prime_x)
    if(allocated(self%p_prime_y)) deallocate(self%p_prime_y)
    if(allocated(self%p0_x)) deallocate(self%p0_x)
    if(allocated(self%p0_y)) deallocate(self%p0_y)
    if(allocated(self%n_arcs_per_cell)) deallocate(self%n_arcs_per_cell)
    if(allocated(self%p_prime_in_cell)) deallocate(self%p_prime_in_cell)
    if(allocated(self%pressure_p_prime)) deallocate(self%pressure_p_prime)
    if(allocated(self%density_p_prime)) deallocate(self%density_p_prime)
    if(allocated(self%u)) deallocate(self%u)
    if(allocated(self%v)) deallocate(self%v)
    if(allocated(self%p)) deallocate(self%p)
    if(allocated(self%edge_vectors)) deallocate(self%edge_vectors)
    if(allocated(self%cell_indices)) deallocate(self%cell_indices)
    if(allocated(self%radius)) deallocate(self%radius)
    if(allocated(self%reference_density)) deallocate(self%reference_density)
    if(allocated(self%reference_u)) deallocate(self%reference_u)
    if(allocated(self%reference_v)) deallocate(self%reference_v)
    if(allocated(self%reference_sound_speed)) deallocate(self%reference_sound_speed)
    if(allocated(self%reference_pressure)) deallocate(self%reference_pressure)
  end subroutine finalize

  subroutine get_reference_state(self)
    !< Find the reference state of the cone only based on the cells that are touched by the cone
    class(mach_cone_collection_t), intent(inout) :: self

    integer(ik) :: i, j, c
    logical, dimension(self%n_neighbor_cells) :: cell_is_supersonic  !< is the neighbor cell supersonic?
    real(rk) :: gamma        !< ideal gas gamma
    real(rk) :: n_cells_real !< # of neighbor cells
    real(rk) :: ref_u        !< Reference u-velocity
    real(rk) :: ref_v        !< Reference v-velocity
    real(rk) :: ref_rho      !< Reference density
    real(rk) :: ref_p        !< Reference pressure
    real(rk), dimension(self%n_neighbor_cells) :: vel          !< Total cell velocity
    real(rk), dimension(self%n_neighbor_cells) :: cs           !< Cell sound speed
    integer(ik) :: n_supersonic_cells !< how many neighbor cells are supersonic? (used to determine if a cell is transonic)

    gamma = eos%get_gamma()

    n_cells_real = real(self%n_neighbor_cells, rk)
    n_supersonic_cells = 0

    !$omp parallel default(none) &
    !$omp shared(self) &
    !$omp reduction(+:ref_rho, ref_u, ref_v, ref_p) &
    !$omp firstprivate(gamma, n_cells_real) &
    !$omp private(i,j,c) &
    !$omp private(cell_is_supersonic, vel, cs) &
    !$omp reduction(+:n_supersonic_cells)
    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        ref_rho = sum(self%recon_rho(:, i, j)) / n_cells_real
        self%reference_density(i, j) = ref_rho
      end do
    end do
    !$omp end do

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        ref_u = sum(self%recon_u(:, i, j)) / n_cells_real
        ref_v = sum(self%recon_v(:, i, j)) / n_cells_real
        self%reference_u(i, j) = ref_u
        self%reference_v(i, j) = ref_v
      end do
    end do
    !$omp end do

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        ref_p = (sum(self%recon_p(:, i, j)) / n_cells_real)
        self%reference_pressure(i, j) = ref_p
      end do
    end do
    !$omp end do

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        self%reference_sound_speed(i, j) = sqrt(gamma * self%reference_pressure(i, j) / self%reference_density(i, j))
      end do
    end do
    !$omp end do
    !$omp barrier

    ! Find the transonic cells for special treatment
    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi

        cell_is_supersonic = .false.

        vel = 0.0_rk
        cs = 0.0_rk
        do c = 1, self%n_neighbor_cells
          vel(c) = sqrt(self%recon_u(c, i, j)**2 + self%recon_v(c, i, j)**2)
          cs(c) = sqrt(gamma * self%recon_p(c, i, j) / self%recon_rho(c, i, j))
          if(vel(c) > cs(c) .or. abs(vel(c) - cs(c)) < epsilon(1.0_rk)) cell_is_supersonic(c) = .true.
        end do

        n_supersonic_cells = count(cell_is_supersonic)

        ! If some, but not all of the cells are supersonic, this is a "transonic" mach cone
        if(n_supersonic_cells > 0 .and. n_supersonic_cells < self%n_neighbor_cells) then
          self%cone_is_transonic(i, j) = .true.
        else
          self%cone_is_transonic(i, j) = .false.
        end if

        ! If the cone is transonic, use the subsonic linearization, e.g. averages from the subsonic cells
        ! if(self%cone_is_transonic(i, j)) then

        !   ! self%reference_sound_speed(i, j) = sum(cs, mask=cell_is_supersonic) / &
        !   !                                    real(n_supersonic_cells, rk)
        !   ! self%reference_density(i, j) = sum(self%recon_rho(:, i, j), mask=cell_is_supersonic) / &
        !   !                                real(n_supersonic_cells, rk)
        !   ! self%reference_u(i, j) = sum(self%recon_u(:, i, j), mask=cell_is_supersonic) / &
        !   !                                 real(n_supersonic_cells, rk)
        !   ! self%reference_v(i, j) = sum(self%recon_v(:, i, j), mask=cell_is_supersonic) / &
        !   !                                 real(n_supersonic_cells, rk)
        !   ! self%reference_pressure(i, j) = sum(self%recon_p(:, i, j), mask=cell_is_supersonic) / &
        !   !                                 real(n_supersonic_cells, rk)

        !   self%reference_sound_speed(i, j) = sum(cs, mask=.not.cell_is_supersonic) / &
        !                                      real(self%n_neighbor_cells - n_supersonic_cells, rk)
        !   self%reference_density(i, j) = sum(self%recon_rho(:, i, j), mask=.not.cell_is_supersonic) / &
        !                                  real(self%n_neighbor_cells - n_supersonic_cells, rk)
        !   ! self%reference_u(i, j) = sum(self%recon_u(:, i, j), mask=.not.cell_is_supersonic) / &
        !   !                                 real(self%n_neighbor_cells - n_supersonic_cells, rk)
        !   ! self%reference_v(i, j) = sum(self%recon_v(:, i, j), mask=.not.cell_is_supersonic) / &
        !   !                                 real(self%n_neighbor_cells - n_supersonic_cells, rk)

        !   ! self%reference_pressure(i, j) = sum(self%recon_p(:, i, j), mask=.not.cell_is_supersonic) / &
        !   !                                       real(self%n_neighbor_cells - n_supersonic_cells, rk)

        !   ! self%reference_sound_speed(i, j) = minval(cs, mask=.not. cell_is_supersonic)
        !   ! self%reference_density(i, j) = self%recon_rho(minloc(cs, dim=1), i, j)

        ! end if
      end do
    end do
    !$omp end do

    !$omp end parallel

  end subroutine get_reference_state

  subroutine get_transonic_cone_extents(self, i, j, origin, radius)
    !< If the set of neighbor cells are transonic (some supersonic, some subsonic),
    !< then the cone needs to be extended so that it encorporates both supersonic and subsonic cells.
    class(mach_cone_collection_t), intent(inout) :: self

    integer(ik), intent(in) :: i !< i index of the transonic cone
    integer(ik), intent(in) :: j !< j index of the transonic cone
    real(rk), dimension(2), intent(out) :: origin  !< origin of the new cone
    real(rk), intent(out) :: radius                !< radius of the new cone

    real(rk), dimension(2) :: origin_1v3 !< (x,y) origin of the circle that superscribes circles 1 and 3
    real(rk), dimension(2) :: origin_2v4 !< (x,y) origin of the circle that superscribes circles 2 and 4
    real(rk) :: radius_1v3               !< radius of the circle that superscribes circles 1 and 3
    real(rk) :: radius_2v4               !< radius of the circle that superscribes circles 1 and 3

    real(rk), dimension(2, 2) :: origins_to_compare !< ((x,y), cell); set of origins to compare between 2 circle
    real(rk), dimension(2) :: radii_to_compare      !< (cell); set of radii to compare between 2 circle

    real(rk), dimension(self%n_neighbor_cells) :: cone_radii      !< (cell); radius of each cone based on the neighbor cell state
    real(rk), dimension(2, self%n_neighbor_cells) :: cone_origins !< ((x,y), cell); origin of each cone based on the neighbor cell state

    integer(ik) :: c !< neighbor cell index
    real(rk) :: cs  !< single cell sound speed
    real(rk) :: vel !< single cell velocity
    real(rk) :: gamma !< EOS gamma

    gamma = eos%get_gamma()

    ! Get the cone size based on each neighbor cell's state
    do c = 1, self%n_neighbor_cells
      vel = sqrt(self%recon_u(c, i, j)**2 + self%recon_v(c, i, j)**2)
      cs = sqrt(gamma * self%recon_p(c, i, j) / self%recon_rho(c, i, j))

      cone_radii(c) = self%tau * cs
      cone_origins(1, c) = self%p0_x(i, j) - self%tau * self%recon_u(c, i, j)
      cone_origins(2, c) = self%p0_y(i, j) - self%tau * self%recon_v(c, i, j)
    end do

    select case(self%n_neighbor_cells)
    case(2)
      call super_circle(origins=cone_origins, radii=cone_radii, &
                        new_origin=origin, new_radius=radius)
    case(4)
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
    case default
      error stop "Error in mach_cone_collection_t%(get_transonic_cone_extents)"// &
        ", unsupported value of self%n_neighbor_cells"
    end select

  end subroutine get_transonic_cone_extents

  subroutine get_p_prime_quantities(self, recon_operator)
    !< Find the density and pressure at P' for the E0 operator to use
    class(mach_cone_collection_t), intent(inout) :: self
    class(abstract_reconstruction_t), intent(in) :: recon_operator !< Reconstruction operator used to find P' values

    integer(ik) :: i, j, c
    integer(ik) :: pp_i !< P' i
    integer(ik) :: pp_j !< P' j
    integer(ik) :: n_p_prime_cells !< how many neighbor cells claim that P' is contained in it?
    ! real(rk) :: pressure_p_prime
    ! real(rk) :: density_p_prime
    real(rk), dimension(self%n_neighbor_cells) :: pressure_p_prime !< temp variable for finding p(P') when multiple cells claim that P' is in it
    real(rk), dimension(self%n_neighbor_cells) :: density_p_prime  !< temp variable for finding rho(P') when multiple cells claim that P' is in it

    i = 0
    j = 0
    c = 0
    pp_i = 0 !< P' i
    pp_j = 0 !< P' j

    !$omp parallel default(none) &
    !$omp reduction(max:density_p_prime, pressure_p_prime) &
    !$omp private(n_p_prime_cells) &
    !$omp shared(self, recon_operator) &
    !$omp private(i,j,c, pp_i, pp_j)
    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi

        n_p_prime_cells = count(self%p_prime_in_cell(:, i, j))

        ! Find rho(P') based on the cells that report that they contain the
        ! point P'(x,y). Sometimes P' is on an edge, so more than one cell reports
        ! that P' is contained in it. If so, then take an average of the P' cells
        if(n_p_prime_cells == self%n_neighbor_cells) then

          ! If all say P', then P' is colocated with P0; P' is then just the reference value, which
          ! is just an average
          self%density_p_prime(i, j) = self%reference_density(i, j)
          self%pressure_p_prime(i, j) = self%reference_pressure(i, j)

        else if(n_p_prime_cells > 1 .and. n_p_prime_cells < self%n_neighbor_cells) then
          ! Take an average of all the P' cells
          density_p_prime = 0.0_rk
          pressure_p_prime = 0.0_rk
          do c = 1, self%n_neighbor_cells
            if(self%p_prime_in_cell(c, i, j)) then
              pp_i = self%cell_indices(1, c, i, j)
              pp_j = self%cell_indices(2, c, i, j)
              density_p_prime(c) = recon_operator%reconstruct_at_point(i=pp_i, j=pp_j, &
                                                                       x=self%p_prime_x(i, j), &
                                                                       y=self%p_prime_y(i, j), &
                                                                       var='rho')
              pressure_p_prime(c) = recon_operator%reconstruct_at_point(i=pp_i, j=pp_j, &
                                                                        x=self%p_prime_x(i, j), &
                                                                        y=self%p_prime_y(i, j), &
                                                                        var='p')
            end if
          end do

          ! Set P' to take the state of cell with the highest pressure
          self%pressure_p_prime(i, j) = sum(pressure_p_prime) / real(n_p_prime_cells, rk)
          self%density_p_prime(i, j) = sum(density_p_prime) / real(n_p_prime_cells, rk)

          ! write(*,'(a, 4(es16.6))') "  p(P'): ", pressure_p_prime
          ! write(*,'(a, 4(es16.6))') "rho(P'): ", density_p_prime
          ! write(*, '(a, 2(es16.6))') " actuals rho, p -> ", self%density_p_prime(i, j), self%pressure_p_prime(i, j)
          ! print*
        else ! only 1 cell reports that P' is in it, so just reconstruct at that single point
          pp_i = self%p_prime_ij(1, i, j)
          pp_j = self%p_prime_ij(2, i, j)
          self%density_p_prime(i, j) = recon_operator%reconstruct_at_point(i=pp_i, j=pp_j, &
                                                                           x=self%p_prime_x(i, j), &
                                                                           y=self%p_prime_y(i, j), &
                                                                           var='rho')
          self%pressure_p_prime(i, j) = recon_operator%reconstruct_at_point(i=pp_i, j=pp_j, &
                                                                            x=self%p_prime_x(i, j), &
                                                                            y=self%p_prime_y(i, j), &
                                                                            var='p')
        end if
      end do
    end do
    !$omp end do
    !$omp end parallel

    ! !$omp parallel default(none) &
    ! !$omp reduction(+:density_p_prime, pressure_p_prime) &
    ! !$omp private(n_p_prime_cells) &
    ! !$omp shared(self) &
    ! !$omp private(i,j,c)
    ! !$omp do
    ! do j = self%jlo, self%jhi
    !   do i = self%ilo, self%ihi
    !     n_p_prime_cells = count(self%p_prime_in_cell(:, i, j))
    !     if(n_p_prime_cells > 1) then
    !       pressure_p_prime = sum(self%recon_p(:, i, j), mask=self%p_prime_in_cell(:, i, j)) / real(n_p_prime_cells, rk)
    !       density_p_prime = sum(self%recon_rho(:, i, j), mask=self%p_prime_in_cell(:, i, j)) / real(n_p_prime_cells, rk)
    !       self%pressure_p_prime(i, j) = pressure_p_prime
    !       self%density_p_prime(i, j) = density_p_prime
    !     else
    !       do c = 1, self%n_neighbor_cells
    !         if(self%p_prime_in_cell(c, i, j)) then
    !           self%density_p_prime(i, j) = self%recon_rho(c, i, j)
    !           self%pressure_p_prime(i, j) = self%recon_p(c, i, j)
    !         end if
    !       end do
    !     end if

    !   end do
    ! end do
    ! !$omp end do
    ! !$omp end parallel
  end subroutine get_p_prime_quantities

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

    integer(ik) :: idx
    !< index conversion to go from ((arc1, arc2), (cell 1:N_CELLS), i, j) to ((cell_1_arc_1, cell_N_arc_1:cell_1_arc_N, cell_N_arc_N), i, j)

    integer(ik) :: i, j
    integer(ik) :: idx_max
    integer(ik), dimension(:, :, :), allocatable :: n_arcs_per_cell !< how many arc lengths in each neighbor

    idx_max = 2 * self%n_neighbor_cells

    associate(ilo=>self%ilo, ihi=>self%ihi, jlo=>self%jlo, jhi=>self%jhi)
      allocate(sin_theta_ib(idx_max, ilo:ihi, jlo:jhi))
      allocate(cos_theta_ib(idx_max, ilo:ihi, jlo:jhi))
      allocate(sin_theta_ie(idx_max, ilo:ihi, jlo:jhi))
      allocate(cos_theta_ie(idx_max, ilo:ihi, jlo:jhi))
      allocate(n_arcs_per_cell(self%n_neighbor_cells, ilo:ihi, jlo:jhi))
    end associate

    call self%find_arc_angles(theta_ib, theta_ie, n_arcs_per_cell)

    self%n_arcs_per_cell = n_arcs_per_cell

    !$omp parallel default(none) &
    !$omp shared(self, theta_ib, theta_ie) &
    !$omp shared(sin_theta_ib, sin_theta_ie) &
    !$omp shared(cos_theta_ib, cos_theta_ie) &
    !$omp firstprivate(idx_max) &
    !$omp private(i, j, idx)

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        do idx = 1, idx_max
          sin_theta_ib(idx, i, j) = sin(theta_ib(idx, i, j))
          cos_theta_ib(idx, i, j) = cos(theta_ib(idx, i, j))
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        do idx = 1, idx_max
          sin_theta_ie(idx, i, j) = sin(theta_ie(idx, i, j))
          cos_theta_ie(idx, i, j) = cos(theta_ie(idx, i, j))
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        do idx = 1, idx_max
          self%dtheta(idx, i, j) = abs(theta_ie(idx, i, j) - theta_ib(idx, i, j))
          ! diff = abs(theta_ie(idx, i, j) - theta_ib(idx, i, j))
          ! if (abs(diff) < epsilon(1.0_rk)) then
          !   self%dtheta(idx, i, j) = 0.0_rk
          ! else
          !   self%dtheta(idx, i, j) = diff
          ! end if
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        do idx = 1, idx_max
          if(abs(sin_theta_ib(idx, i, j)) < 1e-14_rk) sin_theta_ib(idx, i, j) = 0.0_rk
          if(abs(cos_theta_ib(idx, i, j)) < 1e-14_rk) cos_theta_ib(idx, i, j) = 0.0_rk
          if(abs(sin_theta_ie(idx, i, j)) < 1e-14_rk) sin_theta_ie(idx, i, j) = 0.0_rk
          if(abs(cos_theta_ie(idx, i, j)) < 1e-14_rk) cos_theta_ie(idx, i, j) = 0.0_rk
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        do idx = 1, idx_max
          self%sin_dtheta(idx, i, j) = sin_theta_ie(idx, i, j) - sin_theta_ib(idx, i, j)
          ! diff = sin_theta_ie(idx, i, j) - sin_theta_ib(idx, i, j)
          ! if (abs(diff) < epsilon(1.0_rk)) then
          !   self%sin_dtheta(idx, i, j) = 0.0_rk
          ! else
          !   self%sin_dtheta(idx, i, j) = diff
          ! end if
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        do idx = 1, idx_max
          self%cos_dtheta(idx, i, j) = cos_theta_ie(idx, i, j) - cos_theta_ib(idx, i, j)
          ! diff = cos_theta_ie(idx, i, j) - cos_theta_ib(idx, i, j)
          ! if (abs(diff) < epsilon(1.0_rk)) then
          !   self%cos_dtheta(idx, i, j) = 0.0_rk
          ! else
          !   self%cos_dtheta(idx, i, j) = diff
          ! end if
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        do idx = 1, idx_max
          self%sin_d2theta(idx, i, j) = (2.0_rk * sin_theta_ie(idx, i, j) * cos_theta_ie(idx, i, j)) - &
                                        (2.0_rk * sin_theta_ib(idx, i, j) * cos_theta_ib(idx, i, j))
          ! diff = (2.0_rk * sin_theta_ie(idx, i, j) * cos_theta_ie(idx, i, j)) - &
          !        (2.0_rk * sin_theta_ib(idx, i, j) * cos_theta_ib(idx, i, j))
          ! if (abs(diff) < epsilon(1.0_rk)) then
          !   self%sin_d2theta(idx, i, j) = 0.0_rk
          ! else
          !   self%sin_d2theta(idx, i, j) = diff
          ! end if
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        do idx = 1, idx_max
          self%cos_d2theta(idx, i, j) = (2.0_rk * cos_theta_ie(idx, i, j)**2 - 1.0_rk) - &
                                        (2.0_rk * cos_theta_ib(idx, i, j)**2 - 1.0_rk)
          ! diff =  (2.0_rk * cos_theta_ie(idx, i, j)**2 - 1.0_rk) - &
          !         (2.0_rk * cos_theta_ib(idx, i, j)**2 - 1.0_rk)
          ! if (diff < epsilon(1.0_rk)) then
          !   self%cos_d2theta(idx, i, j) = 0.0_rk
          ! else
          !   self%cos_d2theta(idx, i, j) = diff
          ! end if
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

    associate(ilo=>self%ilo, ihi=>self%ihi, jlo=>self%jlo, jhi=>self%jhi)
      allocate(theta_ib(2 * self%n_neighbor_cells, ilo:ihi, jlo:jhi))
      allocate(theta_ie(2 * self%n_neighbor_cells, ilo:ihi, jlo:jhi))
      allocate(n_arcs_per_cell(self%n_neighbor_cells, ilo:ihi, jlo:jhi))
    end associate

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

      !$omp parallel default(none) &
      !$omp shared(self, n_arcs_per_cell, theta_ib, theta_ie) &
      !$omp private(corner_vector_set, first_vector, second_vector) &
      !$omp private(i, j, arc, c, idx) &
      !$omp private(p_prime, p_0, p_prime_cross_vec_1, vec_2_cross_p_prime, n_arcs, theta_ib_ie)
      !$omp do
      do j = self%jlo, self%jhi
        do i = self%ilo, self%ihi
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
              ! self%p_prime_ij(:, i, j) = [i, j]
              self%p_prime_ij(:, i, j) = self%cell_indices(:, c, i, j)
            else
              self%p_prime_in_cell(c, i, j) = .false.
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
      !$omp end parallel

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

      !$omp parallel default(none) &
      !$omp shared(self, n_arcs_per_cell, theta_ib, theta_ie) &
      !$omp private(midpoint_vector_set, first_vector, second_vector) &
      !$omp private(i, j, arc, c, idx) &
      !$omp private(p_prime, p_0, p_prime_cross_vec_1, vec_2_cross_p_prime, n_arcs, theta_ib_ie)
      !$omp do
      do j = self%jlo, self%jhi
        do i = self%ilo, self%ihi
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
              ! self%p_prime_ij(:, i, j) = [i, j]
              self%p_prime_ij(:, i, j) = self%cell_indices(:, c, i, j)
            else
              self%p_prime_in_cell(c, i, j) = .false.
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
      !$omp end parallel

    case default
      error stop 'Error in mach_cone_collection_t%find_arc_angles(), invalid cone location'
    end select

  end subroutine find_arc_angles

  subroutine sanity_checks(self)
    !< Do some sanity checks to make sure the mach cone is valid
    class(mach_cone_collection_t), intent(in) :: self
    ! real(rk), dimension(self%ni, self%nj) :: arc_sum
    real(rk) :: arc_sum
    integer(ik) :: i, j

    ! idx_max = 2 * self%n_neighbor_cells

    !! $omp parallel workshare reduction(+:arc_sum)
    ! arc_sum = sum(self%dtheta, dim=1)
    !! $omp end parallel workshare

    !FIXME: this never seemed to work in openmp
    !! $omp parallel default(none) &
    !! $omp shared(self, arc_sum) &
    !! $omp private(i,j)
    !! $omp do
    do j = self%jlo, self%jhi
      do i = self%ilo, self%ihi
        arc_sum = sum(self%dtheta(:, i, j))
        if(abs(arc_sum - 2.0_rk * pi) > 1e-14_rk) then
          print *, 'abs(arc_sum - 2 * pi): ', abs(arc_sum - 2 * pi)
          print *, "Cone arcs do not add up to 2pi:", arc_sum
          call self%print(i, j)
          error stop "Cone arcs do not add up to 2pi"
        end if
      end do
    end do
    !! $omp end do
    !! $omp end parallel

  end subroutine sanity_checks

  subroutine print(self, i, j)
    class(mach_cone_collection_t), intent(in) :: self
    integer(ik), intent(in) :: i, j
    integer(ik) :: c, arc, idx
    real(rk) :: vlen

    write(*, '(a)') 'Mach Cone Details'
    write(*, '(a)') '================='

    write(*, '(a)') "location = '"//trim(self%cone_location)//"'"
    write(*, '(a, es16.8)') "tau = ", self%tau
    write(*, '(a, i0, 1x, i0)') "i,j = ", i, j
    write(*, '(a, es16.8)') "radius = ", self%radius(i, j)
    write(*, '(a, es16.8, ",", es16.8, a)') 'cone["P0(x,y)"] = [', &
      self%p0_x(i, j), self%p0_y(i, j), "]"
    write(*, '(a, es16.8, ",", es16.8, a)') 'cone["P'//"'(x,y)"//'"] = [', &
      self%p_prime_x(i, j), self%p_prime_y(i, j), "]"
    write(*, '(a, es16.8, ",", es16.8, a)') "dist =            [", &
      self%p_prime_x(i, j) - self%p0_x(i, j), self%p_prime_y(i, j) - self%p0_y(i, j), "]"

    write(*, '(a, i7,",",i7, a)') 'cone["P'//"'(i,j)"//'"] = [', self%p_prime_ij(:, i, j), "]"
    write(*, '(a, 4(l1, 1x))') "P' in cell: ", self%p_prime_in_cell(:, i, j)
    write(*, '(a, l1)') "Cone is transonic: ", self%cone_is_transonic(i, j)
    write(*, '(a, es16.6)') "rho(P'): ", self%density_p_prime(i, j)
    write(*, '(a, es16.6)') "p(P')  : ", self%pressure_p_prime(i, j)
    print *
    write(*, *) '# Neighbor cell indices'
    do c = 1, self%n_neighbor_cells
      write(*, '(a,i0,a, i0, ", ", i0, a)') "cell_", c, ' = [', self%cell_indices(:, c, i, j), "]"
    end do

    print *
    write(*, '(a)') "# Edge Vectors: head (x,y). Note, the tail is at P0(x,y)"
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
