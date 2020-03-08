module mod_mach_cone
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use math_constants, only: pi, rad2deg
  use mod_globals, only: TINY_DIST
  use mod_floating_point_utils, only: near_zero, equal
  use mod_eos, only: eos
  use mod_vector, only: vector_t, operator(.unitnorm.), operator(.dot.), operator(.cross.), angle_between
  use mod_geometry

  implicit none

  private
  public :: mach_cone_base_t, midpoint_mach_cone_t, quad_corner_mach_cone_t

  !< Type to encapsulate and calculate the geometry of the mach cone for the evolution operators
  type, abstract :: mach_cone_base_t

    real(rk), dimension(:, :), allocatable :: recon_state
    !< ((rho, u, v, p), neighbor cell 1 - 4); Reconstructed value of U at P for each neighbor cell

    logical, dimension(:), allocatable :: cell_is_supersonic
    !< Is each neighbor cell supersonic or not? This is needed for transonic mach cones in
    !< order to apply an entropy fix for transonic rarefaction regions.

    real(rk), dimension(:, :), allocatable :: theta_ie
    !< ((arc_1, arc_2), (cell_1:cell_4)); arc end angle

    real(rk), dimension(:, :), allocatable :: theta_ib
    !< ((arc_1, arc_2), (cell_1:cell_4)); arc begin angle

    integer(ik), dimension(:, :), allocatable :: p_prime_ij
    !< ((i,j), cell_1:cell_4) i,j location of P' or the apex of the Mach cone (global, not relative to P0)

    real(rk), dimension(:, :, :), allocatable :: edge_vectors

    integer(ik), dimension(:), allocatable :: n_arcs_per_cell
    !< (cell_1:cell_4); Number of arcs that each cell contributes (0, 1, or 2)

    logical, dimension(:), allocatable :: p_prime_in_cell
    !< (cell_1:cell_4); is P' inside the control volume?

    real(rk), dimension(:, :, :), allocatable :: arc_primitive_vars
    !< ((rho,u,v,p), (arc_1, arc_2), cell_1:cell_n)

    real(rk), dimension(2) :: p_prime_xy = 0.0_rk !< (x,y); Location of P'
    real(rk), dimension(2) :: p_xy = 0.0_rk       !< (x,y); Location of P
    integer(ik) :: n_neighbor_cells = 0           !< Number of neighbor cells
    real(rk) :: radius = -1.0_rk                  !< Radius of the cone
    real(rk) :: tau = 1.0e-10_rk                  !< Time evolution increment
    real(rk) :: reference_density = 0.0_rk        !< Reference density (e.g. neighbor averaged)
    real(rk) :: reference_u = 0.0_rk              !< Reference x velocity (e.g. neighbor averaged)
    real(rk) :: reference_v = 0.0_rk              !< Reference y velocity (e.g. neighbor averaged)
    real(rk) :: reference_sound_speed = 0.0_rk    !< Reference sound speed (e.g. neighbor averaged)
    real(rk) :: reference_mach_number = 0.0_rk    !< Reference Mach number (e.g. neighbor averaged)
    logical :: cone_is_transonic = .false.        !< Flag to enable special treatment for transonic cones
    character(len=32) :: cone_location = ''       !< Corner or midpoint cone?
  contains
    procedure, private :: init_arrays
    procedure, private :: get_reference_state
    procedure, private :: get_transonic_cone_extents
    procedure, private :: get_cone_extents
    procedure, private :: determine_p_prime_cell
    procedure, private :: sanity_checks
    procedure, pass :: write => write_cone
    generic, public :: write(formatted) => write
  end type mach_cone_base_t

  type, extends(mach_cone_base_t) :: midpoint_mach_cone_t
  contains
    procedure :: initialize => init_midpoint_cone
    final :: finalize_midpoint_cone
  end type midpoint_mach_cone_t

  type, extends(mach_cone_base_t) :: quad_corner_mach_cone_t
  contains
    procedure :: initialize => init_corner_cone
    final :: finalize_corner_cone
  end type quad_corner_mach_cone_t

contains

  subroutine init_midpoint_cone(self, tau, edge_vectors, reconstructed_state, cell_indices, cone_location)
    class(midpoint_mach_cone_t), intent(inout) :: self
    real(rk), intent(in) :: tau  !< time increment, tau -> 0 (very small number)
    character(len=*), intent(in) :: cone_location !< corner, down/up midpoint, left/right midpoint
    real(rk), dimension(2, 2, 2), intent(in) :: edge_vectors
    !< ((x,y), (tail,head), (vector_1:vector_n)); set of vectors that define the cell edges

    integer(ik), dimension(2, 2), intent(in) :: cell_indices
    !< ((i,j), cell_1:cell_n); set of indices for the neighboring cells -> needed to find P' i,j index

    real(rk), dimension(4, 2), intent(in) :: reconstructed_state  !< ((rho,u,v,p), cell); reconstructed state for point P.

    ! Locals
    integer(ik) :: l, i  !< index for looping through neighbor cells
    integer(ik) :: arc  !< index for looping through the # arcs in each neighbor cell
    integer(ik) :: n_arcs
    type(vector_t) :: p_prime_vector
    logical :: p_prime_in_cell
    integer(ik), dimension(2) :: cell_ij
    real(rk), dimension(2) :: vector_1_head, vector_2_head, origin
    real(rk), dimension(2) :: second_vector, first_vector
    real(rk), dimension(2, 2) :: theta_ib_ie

    theta_ib_ie = 0.0_rk
    cell_ij = 0

    call self%init_arrays(n_neighbor_cells=2)
    self%tau = tau
    self%p_xy = [edge_vectors(1, 1, 1), edge_vectors(2, 1, 1)]
    self%cone_location = trim(cone_location)
    self%edge_vectors = edge_vectors
    self%n_neighbor_cells = size(reconstructed_state, dim=2)

    do i = 1, 2
      do l = 1, 4
        self%recon_state(l, i) = reconstructed_state(l, i)
      end do
    end do

    call self%get_reference_state(reconstructed_state)

    if(any(reconstructed_state(1, :) < 0.0_rk)) then
      write(*, *) "Error in cone_t initialization, density in the reconstructed state is < 0"
      write(*, '(a, 4(es10.3,1x))') 'Reconstructed density states (cell_1:cell_n)', reconstructed_state(1, :)
      error stop
    end if

    if(any(reconstructed_state(4, :) < 0.0_rk)) then
      write(*, *) "Error in cone_t initialization, pressure in the reconstructed state is < 0"
      write(*, '(a, 4(es10.3,1x))') 'Reconstructed pressure states (cell_1:cell_n): ', reconstructed_state(4, :)
      error stop
    end if

    if(self%cone_is_transonic) then
      call self%get_transonic_cone_extents(origin=self%p_prime_xy, &
                                           radius=self%radius)
      self%reference_sound_speed = self%radius / self%tau
    else
      call self%get_cone_extents(x_vel=self%reference_u, &
                                 y_vel=self%reference_v, &
                                 sound_speed=self%reference_sound_speed, &
                                 origin=self%p_prime_xy, &
                                 radius=self%radius)
    end if

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

    vector_1_head = [edge_vectors(1, 2, 1), edge_vectors(2, 2, 1)] ! (x,y)
    vector_2_head = [edge_vectors(1, 2, 2), edge_vectors(2, 2, 2)] ! (x,y)
    origin = [edge_vectors(1, 1, 1), edge_vectors(2, 1, 1)] ! (x,y)

    p_prime_vector = vector_t(x=[edge_vectors(1, 1, 1), self%p_prime_xy(1)], &
                              y=[edge_vectors(2, 1, 1), self%p_prime_xy(2)])

    do i = 1, 2
      select case(i)
      case(1) ! Cell 1
        first_vector = vector_2_head
        second_vector = vector_1_head
      case(2) ! Cell 2
        first_vector = vector_1_head
        second_vector = vector_2_head
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
      self%n_arcs_per_cell(i) = n_arcs
      self%theta_ib(:, i) = theta_ib_ie(1, :)
      self%theta_ie(:, i) = theta_ib_ie(2, :)
    end do

    ! Assign the primitive state variables to each arc contained in each cell
    do i = 1, 2
      do arc = 1, self%n_arcs_per_cell(i)
        self%arc_primitive_vars(:, arc, i) = self%recon_state(:, i)
      end do
    end do

    call self%determine_p_prime_cell()
    call self%sanity_checks()
  end subroutine init_midpoint_cone

  subroutine init_corner_cone(self, tau, edge_vectors, reconstructed_state, cell_indices, cone_location)
    class(quad_corner_mach_cone_t), intent(inout) :: self
    real(rk), intent(in) :: tau  !< time increment, tau -> 0 (very small number)
    character(len=*), intent(in) :: cone_location !< corner, down/up midpoint, left/right midpoint
    real(rk), dimension(2, 2, 4), intent(in) :: edge_vectors
    !< ((x,y), (tail,head), (vector_1:vector_n)); set of vectors that define the cell edges

    integer(ik), dimension(2, 4), intent(in) :: cell_indices
    !< ((i,j), cell_1:cell_n); set of indices for the neighboring cells -> needed to find P' i,j index

    real(rk), dimension(4, 4), intent(in) :: reconstructed_state  !< ((rho,u,v,p), cell); reconstructed state for point P.

    ! Locals
    integer(ik) :: i  !< index for looping through neighbor cells
    integer(ik) :: arc  !< index for looping through the # arcs in each neighbor cell
    integer(ik) :: n_arcs
    type(vector_t) :: p_prime_vector
    logical :: p_prime_in_cell
    integer(ik), dimension(2) :: cell_ij
    real(rk), dimension(2) :: vector_1_head, vector_2_head, vector_3_head, vector_4_head, origin
    real(rk), dimension(2) :: second_vector, first_vector
    real(rk), dimension(2, 2) :: theta_ib_ie

    call self%init_arrays(n_neighbor_cells=4)
    self%tau = tau
    self%p_xy = [edge_vectors(1, 1, 1), edge_vectors(2, 1, 1)]
    self%cone_location = trim(cone_location)
    self%edge_vectors = edge_vectors
    call self%get_reference_state(reconstructed_state)

    if(any(reconstructed_state(1, :) < 0.0_rk)) then
      write(*, *) "Error in cone_t initialization, density in the reconstructed state is < 0"
      write(*, '(a, 4(es10.3,1x))') 'Reconstructed density states (cell_1:cell_n)', reconstructed_state(1, :)
      error stop
    end if

    if(any(reconstructed_state(4, :) < 0.0_rk)) then
      write(*, *) "Error in cone_t initialization, pressure in the reconstructed state is < 0"
      write(*, '(a, 4(es10.3,1x))') 'Reconstructed pressure states (cell_1:cell_n): ', reconstructed_state(4, :)
      error stop
    end if

    if(self%cone_is_transonic) then
      call self%get_transonic_cone_extents(origin=self%p_prime_xy, &
                                           radius=self%radius)
      self%reference_sound_speed = self%radius / self%tau
    else
      call self%get_cone_extents(x_vel=self%reference_u, &
                                 y_vel=self%reference_v, &
                                 sound_speed=self%reference_sound_speed, &
                                 origin=self%p_prime_xy, &
                                 radius=self%radius)
    end if

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

    do i = 1, 4
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
      self%n_arcs_per_cell(i) = n_arcs
      self%theta_ib(:, i) = theta_ib_ie(1, :)
      self%theta_ie(:, i) = theta_ib_ie(2, :)
    end do

    ! Assign the primitive state variables to each arc contained in each cell
    do i = 1, 4
      do arc = 1, self%n_arcs_per_cell(i)
        self%arc_primitive_vars(:, arc, i) = self%recon_state(:, i)
      end do
    end do

    call self%determine_p_prime_cell()
    call self%sanity_checks()

  end subroutine init_corner_cone

  subroutine finalize_corner_cone(self)
    type(quad_corner_mach_cone_t), intent(inout) :: self

    if(allocated(self%recon_state)) deallocate(self%recon_state)
    if(allocated(self%cell_is_supersonic)) deallocate(self%cell_is_supersonic)
    if(allocated(self%theta_ie)) deallocate(self%theta_ie)
    if(allocated(self%theta_ib)) deallocate(self%theta_ib)
    if(allocated(self%p_prime_ij)) deallocate(self%p_prime_ij)
    if(allocated(self%edge_vectors)) deallocate(self%edge_vectors)
    if(allocated(self%n_arcs_per_cell)) deallocate(self%n_arcs_per_cell)
    if(allocated(self%p_prime_in_cell)) deallocate(self%p_prime_in_cell)
    if(allocated(self%arc_primitive_vars)) deallocate(self%arc_primitive_vars)
  end subroutine finalize_corner_cone

  subroutine finalize_midpoint_cone(self)
    type(midpoint_mach_cone_t), intent(inout) :: self

    if(allocated(self%recon_state)) deallocate(self%recon_state)
    if(allocated(self%cell_is_supersonic)) deallocate(self%cell_is_supersonic)
    if(allocated(self%theta_ie)) deallocate(self%theta_ie)
    if(allocated(self%theta_ib)) deallocate(self%theta_ib)
    if(allocated(self%p_prime_ij)) deallocate(self%p_prime_ij)
    if(allocated(self%edge_vectors)) deallocate(self%edge_vectors)
    if(allocated(self%n_arcs_per_cell)) deallocate(self%n_arcs_per_cell)
    if(allocated(self%p_prime_in_cell)) deallocate(self%p_prime_in_cell)
    if(allocated(self%arc_primitive_vars)) deallocate(self%arc_primitive_vars)
  end subroutine finalize_midpoint_cone

  subroutine init_arrays(self, n_neighbor_cells)
    class(mach_cone_base_t), intent(inout) :: self
    integer(ik), intent(in) :: n_neighbor_cells

    self%n_neighbor_cells = n_neighbor_cells

    if(.not. allocated(self%recon_state)) allocate(self%recon_state(4, n_neighbor_cells))
    self%recon_state = 0.0_rk

    if(.not. allocated(self%cell_is_supersonic)) allocate(self%cell_is_supersonic(n_neighbor_cells))
    self%cell_is_supersonic = .false.

    if(.not. allocated(self%theta_ie)) allocate(self%theta_ie(2, n_neighbor_cells))
    self%theta_ie = 0.0_rk

    if(.not. allocated(self%theta_ib)) allocate(self%theta_ib(2, n_neighbor_cells))
    self%theta_ib = 0.0_rk

    if(.not. allocated(self%p_prime_ij)) allocate(self%p_prime_ij(2, n_neighbor_cells))
    self%p_prime_ij = 0

    if(.not. allocated(self%edge_vectors)) allocate(self%edge_vectors(2, 2, n_neighbor_cells))
    self%edge_vectors = 0

    if(.not. allocated(self%n_arcs_per_cell)) allocate(self%n_arcs_per_cell(n_neighbor_cells))
    self%n_arcs_per_cell = 0

    if(.not. allocated(self%p_prime_in_cell)) allocate(self%p_prime_in_cell(n_neighbor_cells))
    self%p_prime_in_cell = .false.

    if(.not. allocated(self%arc_primitive_vars)) allocate(self%arc_primitive_vars(4, 2, n_neighbor_cells))
    self%arc_primitive_vars = 0.0_rk

  end subroutine init_arrays

  subroutine write_cone(self, unit, iotype, v_list, iostat, iomsg)
    !< Implementation of `write(*,*) mach_cone_base_t`

    class(mach_cone_base_t), intent(in) :: self     !< cone class
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
    write(unit, '(a, es14.6, ",", es14.6, a)', iostat=iostat, iomsg=iomsg) 'cone["P'//"'(x,y)"//'"] = [', self%p_prime_xy, "]"//new_line('a')
    write(unit, '(a, es14.6, ",", es14.6, a)', iostat=iostat, iomsg=iomsg) &
      "dist = [", self%p_prime_xy - self%p_xy, "]"//new_line('a')

    do i = 1, 4
      write(unit, '(a, i7,",",i7, a, i0, a)', iostat=iostat, iomsg=iomsg) &
        'cone["P'//"'(i,j)"//'"] = [', self%p_prime_ij(:, i), "] # Cell: ", i, new_line('a')
    end do

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) "# Edge Vectors: tail (x,y) -> head (x,y)"//new_line('a')

    do i = 1, self%n_neighbor_cells
      write(unit, '(a, i0, a, 2(es10.3, ", ", es10.3, a))', iostat=iostat, iomsg=iomsg) "vector_", i, " = [[", &
        self%edge_vectors(:, 1, i), "], [", self%edge_vectors(:, 2, i), "] ]"//new_line('a')
    end do

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "reference_density =     ", self%reference_density, new_line('a')
    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "reference_x_velocity =  ", self%reference_u, new_line('a')
    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "reference_y_velocity =  ", self%reference_v, new_line('a')
    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "reference_sound_speed = ", self%reference_sound_speed, new_line('a')
    write(unit, '(a, es10.3, a)', iostat=iostat, iomsg=iomsg) "reference_mach_Number = ", self%reference_mach_number, new_line('a')

    write(unit, '(a)', iostat=iostat) new_line('a')
    write(unit, '(a, 4(i0, 1x),a)', iostat=iostat) '# Valid arcs in each cell: [', self%n_arcs_per_cell, ']'
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) '# Angles:       Theta_ib [deg]         Theta_ie [deg]'//new_line('a')
    do i = 1, 4
      write(unit, '(a, i0, 4(a,f9.5, ",",f9.5), a )', iostat=iostat, iomsg=iomsg) &
        'angles[:,', i, &
        '] = [[ ', rad2deg(self%theta_ib(:, i)), &
        '], [ ', rad2deg(self%theta_ie(:, i)), &
        ']]'//new_line('a')
    end do
    write(unit, '(a, f0.2, a)', iostat=iostat, iomsg=iomsg) 'arc_length_sum = ', &
      rad2deg(sum(self%theta_ie - self%theta_ib)), ' '//new_line('a')

    ! write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    ! write(unit, '(a)', iostat=iostat, iomsg=iomsg) 'Angles:       Theta_ib [rad]         Theta_ie [rad]'//new_line('a')
    ! do i = 1, 4
    !   write(unit, '(a, i0, 4(a,f9.4, ",",f9.4), a )', iostat=iostat, iomsg=iomsg) &
    !     'angles[:,', i, &
    !     '] = [[ ', self%theta_ib(:, i), &
    !     '], [ ', self%theta_ie(:, i), &
    !     ']]'//new_line('a')
    ! end do
    ! write(unit, '(a, f0.2, a)', iostat=iostat, iomsg=iomsg) 'arc_length_sum = ', &
    !   sum(self%theta_ie - self%theta_ib), ' '//new_line('a')

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) 'Cell Primitive Vars State [rho,u,v,p]'//new_line('a')
    do i = 1, 4
      write(unit, '(a, i0, a, 4(es10.3, 1x), a)', iostat=iostat, iomsg=iomsg) &
        'Cell: ', i, ' [ ', self%arc_primitive_vars(:, 1, i), '] (arc 1)'//new_line('a')
      write(unit, '(a, 4(es10.3, 1x), a)', iostat=iostat, iomsg=iomsg) &
        '        [ ', self%arc_primitive_vars(:, 2, i), '] (arc 2)'//new_line('a')
    end do

    ! write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    ! write(unit, '(a)', iostat=iostat, iomsg=iomsg) 'Reconstructed state [rho,u,v,p]'//new_line('a')
    ! do i = 1, 4
    !   write(unit, '(a, i0, a, 4(f6.2, 1x), a)', iostat=iostat, iomsg=iomsg) &
    !     'Cell: ', i, ' [ ', self%reconstructed_state(:, i), ']'//new_line('a')
    ! end do

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
    write(unit, '(a)', iostat=iostat, iomsg=iomsg) "# P' in cell?"//new_line('a')
    do i = 1, 4
      write(unit, '(a, i0, a, l2, a)', iostat=iostat, iomsg=iomsg) 'Cell: ', i, ' [', self%p_prime_in_cell(i), ' ]'//new_line('a')
    end do

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) '# Neighbor cells contributing to the mach cone '//new_line('a')
    write(unit, '(a, i0, a)', iostat=iostat, iomsg=iomsg) 'neighbor_cells = ', self%n_neighbor_cells, new_line('a')
    do i = 1, ubound(self%recon_state, dim=2)
      write(unit, '(a, i0, " = [", 3(es10.3, ","), es10.3,"]", a)', iostat=iostat, iomsg=iomsg) 'cell_', i, self%recon_state(:, i), new_line('a')
    end do

    write(unit, '(a)', iostat=iostat, iomsg=iomsg) new_line('a')
  end subroutine write_cone

  pure subroutine get_cone_extents(self, x_vel, y_vel, sound_speed, origin, radius)
    !< Given a velocity and sound speed, determine the extents of the mach cone

    class(mach_cone_base_t), intent(in) :: self
    real(rk), intent(in) :: x_vel        !< x velocity
    real(rk), intent(in) :: y_vel        !< y velocity
    real(rk), intent(in) :: sound_speed  !< sound speed
    real(rk), dimension(2), intent(out) :: origin      !< origin of the cone/circle
    real(rk), intent(out) :: radius      !< radius of the cone/circle

    ! Locals
    real(rk) :: velocity !< Magnitude of the velocity vector

    velocity = sqrt(x_vel**2 + y_vel**2)

    associate(x=>self%p_xy(1), y=>self%p_xy(2), &
              u=>x_vel, v=>y_vel, tau=>self%tau)
      origin = [x - tau * u, y - tau * v]
      radius = sound_speed * tau

      ! For very small distances, just make it coincide with P
      if(abs(origin(1) - x) < TINY_DIST) origin(1) = x
      if(abs(origin(2) - y) < TINY_DIST) origin(2) = y
    end associate
  end subroutine get_cone_extents

  subroutine get_transonic_cone_extents(self, origin, radius)
    !< If the set of neighbor cells are transonic (some supersonic, some subsonic),
    !< then the cone needs to be extended so that it encorporates both supersonic and subsonic cells.

    class(mach_cone_base_t), intent(inout) :: self
    real(rk), dimension(2), intent(out) :: origin  !< origin of the new cone
    real(rk), intent(out) :: radius                !< radius of the new cone

    real(rk), dimension(2) :: origin_1v3, origin_2v4
    real(rk) :: radius_1v3, radius_2v4
    real(rk), dimension(2, 2) :: origins_to_compare
    real(rk), dimension(2) :: radii_to_compare

    integer(ik) :: i
    real(rk) :: sound_speed
    logical :: is_inside
    integer(ik) :: largest_cone_idx
    !< index of the cone with the largets radius

    real(rk), dimension(:), allocatable :: cone_radii
    !< (cell); radius of each cone based on the neighbor cell state

    real(rk), dimension(:, :), allocatable :: cone_origins
    !< ((x,y), cell); origin of each cone based on the neighbor cell state

    associate(rho=>self%recon_state(1, :), &
              u=>self%recon_state(2, :), &
              v=>self%recon_state(3, :), &
              p=>self%recon_state(4, :), &
              n=>self%n_neighbor_cells)

      allocate(cone_radii(n))
      allocate(cone_origins(2, n))

      ! print*, 'Getting transonic cone extents'
      ! Get the cone size based on each neighbor cell's state
      do i = 1, n
        sound_speed = eos%sound_speed(pressure=p(i), density=rho(i))
        call self%get_cone_extents(x_vel=u(i), y_vel=v(i), sound_speed=sound_speed, &
                                   origin=cone_origins(:, i), radius=cone_radii(i))
        ! print*, "Cell: ", i, "Cs:", sound_speed, "V: ", sqrt(u(i)**2 + v(i)**2), "Supersonic", self%cell_is_supersonic(i)
      end do
    end associate
    select case(self%n_neighbor_cells)
    case(2)  ! 2 circles to compare

      ! print*, 'cone_origins', cone_origins
      ! print*, 'cone_radii', cone_radii
      call super_circle(origins=cone_origins, radii=cone_radii, &
                        new_origin=origin, new_radius=radius)
      ! print*, 'new cone_origins', origin
      ! print*, 'new cone_radii', radius
      ! error stop
    case(4)  ! 4 circles to compare

      ! Since the circle comparison is only 2 at a time, this
      ! does the diagonal cells first (1v3 and 2v4) then then
      ! compares the results from those for the final circle

      ! Compare cells 1 and 3
      ! print*, 'cone_origins', cone_origins
      ! print*, 'cone_radii', cone_radii
      origins_to_compare(:, 1) = cone_origins(:, 1)
      origins_to_compare(:, 2) = cone_origins(:, 3)
      radii_to_compare = [cone_radii(1), cone_radii(3)]
      call super_circle(origins=origins_to_compare, radii=radii_to_compare, &
                        new_origin=origin_1v3, new_radius=radius_1v3)
      ! print*, 'origin_1v3', origin_1v3
      ! print*, 'radius_1v3', radius_1v3

      ! Compare cells 2 and 4
      origins_to_compare(:, 1) = cone_origins(:, 2)
      origins_to_compare(:, 2) = cone_origins(:, 4)
      radii_to_compare = [cone_radii(2), cone_radii(4)]
      call super_circle(origins=origins_to_compare, radii=radii_to_compare, &
                        new_origin=origin_2v4, new_radius=radius_2v4)
      ! print*, 'origin_2v4', origin_2v4
      ! print*, 'radius_2v4', radius_2v4

      ! Compare result from 1 vs 3 and 2 vs 4
      origins_to_compare(:, 1) = origin_1v3
      origins_to_compare(:, 2) = origin_2v4
      radii_to_compare = [radius_1v3, radius_2v4]
      call super_circle(origins=origins_to_compare, radii=radii_to_compare, &
                        new_origin=origin, new_radius=radius)

      ! print*, 'origins', origin
      ! print*, 'radius', radius
    end select

    deallocate(cone_origins)
    deallocate(cone_radii)
  end subroutine get_transonic_cone_extents

  subroutine get_reference_state(self, reconstructed_state)
    !< Find the reference state of the cone only based on the cells that are touched by the cone.
    !< Sometimes when a cone is near a large density jump, a skewed reference state can cause
    !< negative densities and pressures. This aims to alleviate that problem...

    class(mach_cone_base_t), intent(inout) :: self
    real(rk), dimension(:, :), intent(in) :: reconstructed_state !< [rho, u, v, a]

    integer(ik) :: ref_state_cell
    real(rk) :: ave_p
    real(rk), dimension(:), allocatable, contiguous :: mach_number, sound_speed
    integer(ik) :: i
    real(rk) :: gamma

    allocate(mach_number(self%n_neighbor_cells))
    allocate(sound_speed(self%n_neighbor_cells))

    ! sound_speed = eos%sound_speed(pressure=self%recon_state(4, :), density=self%recon_state(1, :))
    ! mach_number = sqrt(self%recon_state(2, :)**2 + self%recon_state(3, :)**2) / sound_speed
    gamma = eos%get_gamma()
    do i = 1, self%n_neighbor_cells
      sound_speed(i) = sqrt(gamma * self%recon_state(4, i) / self%recon_state(1, i))
    end do
    do i = 1, self%n_neighbor_cells
      mach_number(i) = sqrt(self%recon_state(2, i)**2 + self%recon_state(3, i)**2) / sound_speed(i)
    end do

    associate(rho=>self%recon_state(1, :), &
              u=>self%recon_state(2, :), &
              v=>self%recon_state(3, :), &
              p=>self%recon_state(4, :), &
              n=>real(self%n_neighbor_cells, rk))

      where(mach_number > 1.0_rk) self%cell_is_supersonic = .true.

      ! The cone is transonic if there is a combo of super/subsonic cells
      if(count(self%cell_is_supersonic(1:self%n_neighbor_cells)) == n .or. &
         count(self%cell_is_supersonic(1:self%n_neighbor_cells)) == 0) then
        self%cone_is_transonic = .false.
      else
        self%cone_is_transonic = .true.
      end if

      self%reference_density = sum(rho) / n
      self%reference_u = sum(u) / n
      self%reference_v = sum(v) / n
      ave_p = sum(p) / n
      self%reference_sound_speed = eos%sound_speed(pressure=ave_p, &
                                                   density=self%reference_density)
    end associate

    ! ! Tiny velocity fix
    ! if(abs(self%reference_u) < 1e-10_rk) self%reference_u = 0.0_rk
    ! if(abs(self%reference_v) < 1e-10_rk) self%reference_v = 0.0_rk

    self%reference_mach_number = sqrt(self%reference_u**2 + self%reference_v**2) / &
                                 self%reference_sound_speed

    deallocate(mach_number)
    deallocate(sound_speed)
  end subroutine get_reference_state

  pure logical function determine_if_p_prime_is_in_cell(origin, vec_1_head, vec_2_head, p_prime_vector) result(in_cell)
    !< Implementation of whether the P' point is inside the current cell/control volume. This
    !< uses the cross product of 2 vectors in 2d, which gives a scalar
    type(vector_t), intent(in) :: p_prime_vector
    real(rk), dimension(2), intent(in) :: origin     !< (x,y) origin of both vectors
    real(rk), dimension(2), intent(in) :: vec_1_head !< (x,y) head of the 1st vector
    real(rk), dimension(2), intent(in) :: vec_2_head !< (x,y) head of the 2nd vector

    type(vector_t) :: edge_vector_1, edge_vector_2

    in_cell = .false.
    edge_vector_1 = vector_t(x=[origin(1), vec_1_head(1)], &
                             y=[origin(2), vec_1_head(2)])

    edge_vector_2 = vector_t(x=[origin(1), vec_2_head(1)], &
                             y=[origin(2), vec_2_head(2)])

    ! In the text is has >= instead of <= for some reason, but the following works
    ! like it's supposed to.
    if((p_prime_vector.cross.edge_vector_1) <= 0.0_rk .and. &
       (edge_vector_2.cross.p_prime_vector) <= 0.0_rk) then
      in_cell = .true.
    else
      in_cell = .false.
    end if
  end function determine_if_p_prime_is_in_cell

  subroutine determine_p_prime_cell(self)
    !< Sometimes when u and v tilde (reference velocities) are 0, P' and P are collocated. When this
    !< is the case, P' needs to be chosen from one of the neighbor cells. This subroutine scans all
    !< the neighbor cells and picks the one with the highest pressure (if any). If not, then it juse
    !< choses the first cell from the list of cells that "contain" P'

    class(mach_cone_base_t), intent(inout) :: self
    integer(ik) :: p_prime_cell, i
    real(rk) :: max_p
    integer(ik) :: max_p_idx
    integer(ik), dimension(2) :: p_prime_ij

    max_p = 0.0_rk
    p_prime_ij = [0, 0]

    if(count(self%p_prime_in_cell) > 1) then
      do i = 1, self%n_neighbor_cells
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
    class(mach_cone_base_t), intent(in) :: self
    real(rk) :: arc_sum

    arc_sum = sum(self%theta_ie - self%theta_ib)

    if(count(self%n_arcs_per_cell >= 2) > 1) then
      error stop "Too many arcs in the mach cone (count(cone%n_arcs >= 2) > 1)"
    end if

    if(.not. equal(arc_sum, 2 * pi, 1e-6_rk)) then
      print *, "Cone arcs do not add up to 2pi: ", arc_sum
      print *, self
      error stop "Cone arcs do not add up to 2pi"
    end if

    if(self%radius < 0.0_rk) then
      write(*, *) "Error: cone radius < 0"
      error stop "Error: cone radius < 0"
    end if

  end subroutine sanity_checks

end module mod_mach_cone
