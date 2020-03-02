module mod_cone
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use math_constants, only: pi, rad2deg
  use mod_globals, only: TINY_DIST
  use mod_floating_point_utils, only: near_zero, equal
  use mod_eos, only: eos
  use mod_vector, only: vector_t, operator(.unitnorm.), operator(.dot.), operator(.cross.), angle_between
  use mod_geometry

  implicit none

  type :: cone_t
    !< Type to encapsulate and calculate the geometry of the mach cone for the evolution operators
    real(rk), dimension(2, 4) :: theta_ie = 0.0_rk
    !< ((arc_1, arc_2), (cell_1:cell_4)); arc end angle

    real(rk), dimension(2, 4) :: theta_ib = 0.0_rk
    !< ((arc_1, arc_2), (cell_1:cell_4)); arc begin angle

    real(rk), dimension(2) :: p_prime_xy = 0.0_rk
    !< (x,y); Location of P'

    real(rk), dimension(2) :: p_xy = 0.0_rk
    !< (x,y); Location of P

    integer(ik), dimension(2, 4) :: p_prime_ij = 0
    !< ((i,j), cell_1:cell_4) i,j location of P' or the apex of the Mach cone (global, not relative to P0)

    integer(ik) :: n_neighbor_cells = 0
    !< How many cells influence this cone (up to 4 for now)

    real(rk) :: radius = -1.0_rk
    !< Radius of the cone

    real(rk), dimension(:, :, :), allocatable :: edge_vectors

    integer(ik), dimension(4) :: n_arcs = 0
    !< (cell_1:cell_4); Number of arcs that each cell contributes (0, 1, or 2)

    logical, dimension(4) :: p_prime_in_cell = .false.
    !< (cell_1:cell_4); is P' inside the control volume?

    real(rk), dimension(4, 2, 4) :: arc_primitive_vars = 0.0_rk
    !< ((rho,u,v,p), (arc_1, arc_2), cell_1:cell_n)

    real(rk) :: reference_density = 0.0_rk
    real(rk) :: reference_u = 0.0_rk
    real(rk) :: reference_v = 0.0_rk
    real(rk) :: reference_sound_speed = 0.0_rk
    real(rk) :: reference_mach_number = 0.0_rk

    real(rk) :: tau = 1.0e-10_rk
    !< Time evolution increment

    logical, dimension(4) :: cell_is_supersonic
    !< Is each neighbor cell supersonic or not? This is needed for transonic mach cones in
    !< order to apply an entropy fix for transonic rarefaction regions.

    logical :: cone_is_transonic = .false.

    real(rk), dimension(:, :), allocatable :: recon_state
    !< ((rho, u, v, p), neighbor cell 1 - 4); Reconstructed value of U at P for each neighbor cell

    character(len=32) :: cone_location
  contains
    procedure, nopass, private :: determine_if_p_prime_is_in_cell
    procedure, private :: determine_p_prime_cell
    procedure, private :: sanity_checks
    procedure, private :: get_transonic_cone_extents
    procedure, private :: get_cone_extents
    procedure, private :: get_reference_state
    procedure, pass :: write => write_cone
    generic, public :: write(formatted) => write
  end type

contains

  type(cone_t) function new_cone(tau, edge_vectors, reconstructed_state, cell_indices, cone_location)
    !< Constructor for the Mach cone type
    real(rk), intent(in) :: tau
    !< time increment, tau -> 0 (very small number)

    character(len=*), intent(in) :: cone_location !< corner, down/up midpoint, left/right midpoint

    real(rk), dimension(:, :, :), intent(in) :: edge_vectors
    !< ((x,y), (tail,head), (vector_1:vector_n)); set of vectors that define the cell edges

    integer(ik), dimension(:, :), intent(in) :: cell_indices
    !< ((i,j), cell_1:cell_n); set of indices for the neighboring cells -> needed to find P' i,j index

    real(rk), dimension(:, :), intent(in) :: reconstructed_state
    !< ((rho,u,v,p), cell); reconstructed state for point P.
    !< P has a different reconstruction for each cell

    integer(ik) :: n_total_vectors  !< number of edge vectors (should be 2 or 4)
    integer(ik) :: neighbor_cell  !< index for looping through neighbor cells
    integer(ik), dimension(2, 4) :: edge_vector_ordering  !< index order for each set of vectors
    real(rk), dimension(2, 2, 2) :: single_cell_edge_vectors

    type(vector_t) :: p_prime_vector
    logical :: p_prime_in_cell
    integer(ik) :: n_arcs, arc
    integer(ik), dimension(2) :: cell_ij
    real(rk), dimension(2, 2) :: theta_ib_ie

    integer(ik) :: tilde_i

    p_prime_in_cell = .false.
    n_arcs = 0
    cell_ij = 0
    theta_ib_ie = 0.0_rk
    n_total_vectors = 0
    neighbor_cell = 0
    edge_vector_ordering = 0
    single_cell_edge_vectors = 0.0_rk
    n_total_vectors = size(edge_vectors, dim=3)
    new_cone%tau = tau

    new_cone%p_xy = [edge_vectors(1, 1, 1), edge_vectors(2, 1, 1)]

    allocate(new_cone%edge_vectors, mold=edge_vectors)
    new_cone%edge_vectors = edge_vectors
    new_cone%cone_location = trim(cone_location)

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

    call new_cone%get_reference_state(reconstructed_state)

    if(new_cone%cone_is_transonic) then
      call new_cone%get_transonic_cone_extents(origin=new_cone%p_prime_xy, &
                                               radius=new_cone%radius)
    else
      call new_cone%get_cone_extents(x_vel=new_cone%reference_u, &
                                     y_vel=new_cone%reference_v, &
                                     sound_speed=new_cone%reference_sound_speed, &
                                     origin=new_cone%p_prime_xy, &
                                     radius=new_cone%radius)
    end if

    ! the order of the vectors matters, mainly for the cross product to determine
    ! if P' is in the neighboring cell or not
    select case(n_total_vectors)
    case(2)
      edge_vector_ordering(:, 1:2) = reshape([[2, 1],[1, 2]], shape=[2, 2])
      new_cone%n_neighbor_cells = 2
    case(4)
      edge_vector_ordering = reshape([[4, 1],[1, 2],[2, 3],[3, 4]], shape=[2, 4])
      new_cone%n_neighbor_cells = 4
    case default
      error stop 'Code not set up yet to handle mach cones with other than 2 or 4 edge vectors'
    end select

    p_prime_vector = vector_t(x=[edge_vectors(1, 1, 1), new_cone%p_prime_xy(1)], &
                              y=[edge_vectors(2, 1, 1), new_cone%p_prime_xy(2)])

    ! Loop through each neighbor cell and determine intersections and angles
    do neighbor_cell = 1, new_cone%n_neighbor_cells
      p_prime_in_cell = .false.

      cell_ij = cell_indices(:, neighbor_cell)  ! cell_indices is indexed via ((i,j), cell_1:cell_n)

      ! indexing for edge_vectors is ((x,y), (tail,head), (vector_1:vector_n))
      single_cell_edge_vectors(:, :, 1) = edge_vectors(:, :, edge_vector_ordering(1, neighbor_cell))
      single_cell_edge_vectors(:, :, 2) = edge_vectors(:, :, edge_vector_ordering(2, neighbor_cell))

      p_prime_in_cell = determine_if_p_prime_is_in_cell(single_cell_edge_vectors, p_prime_vector)
      new_cone%p_prime_in_cell(neighbor_cell) = p_prime_in_cell
      if(p_prime_in_cell) new_cone%p_prime_ij(:, neighbor_cell) = cell_ij

      call get_arc_segments(lines=single_cell_edge_vectors, origin_in_cell=p_prime_in_cell, &
                            circle_xy=new_cone%p_prime_xy, circle_radius=new_cone%radius, &
                            arc_segments=theta_ib_ie, n_arcs=n_arcs)

      if(n_arcs > 0) then
        do arc = 1, n_arcs
          new_cone%arc_primitive_vars(:, arc, neighbor_cell) = reconstructed_state(:, neighbor_cell)

          ! Small number fix
          if(new_cone%arc_primitive_vars(2, arc, neighbor_cell) < 5e-16_rk) then
            new_cone%arc_primitive_vars(2, arc, neighbor_cell) = 0.0_rk
          end if

          if(new_cone%arc_primitive_vars(3, arc, neighbor_cell) < 5e-16_rk) then
            new_cone%arc_primitive_vars(3, arc, neighbor_cell) = 0.0_rk
          end if
        end do
      end if

      ! If all the arcs already add up to 2pi, then skip the remaining neighbor cells
      if(equal(sum(new_cone%theta_ie - new_cone%theta_ib), 2.0_rk * pi, epsilon=1e-10_rk)) exit
      ! //TODO: test the above statement in unit testing
      new_cone%n_arcs(neighbor_cell) = n_arcs
      new_cone%theta_ib(:, neighbor_cell) = theta_ib_ie(1, :)
      new_cone%theta_ie(:, neighbor_cell) = theta_ib_ie(2, :)

    end do

    call new_cone%determine_p_prime_cell()
    call new_cone%sanity_checks()

    ! if (new_cone%cone_is_transonic) print*, new_cone
  end function new_cone

  subroutine determine_p_prime_cell(self)
    !< Sometimes when u and v tilde (reference velocities) are 0, P' and P are collocated. When this
    !< is the case, P' needs to be chosen from one of the neighbor cells. This subroutine scans all
    !< the neighbor cells and picks the one with the highest pressure (if any). If not, then it juse
    !< choses the first cell from the list of cells that "contain" P'

    class(cone_t), intent(inout) :: self
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
    class(cone_t), intent(in) :: self

    if(count(self%n_arcs >= 2) > 1) then
      error stop "Too many arcs in the mach cone (count(cone%n_arcs >= 2) > 1)"
    end if

    if(.not. equal(sum(self%theta_ie - self%theta_ib), 2.0_rk * pi, epsilon=1e-10_rk)) then
      print *, "Cone arcs do not add up to 2pi: ", sum(self%theta_ie - self%theta_ib)
      print *, self
      error stop "Cone arcs do not add up to 2pi"
    end if

    if(self%radius < 0.0_rk) then
      write(*, *) "Error: cone radius < 0"
      error stop "Error: cone radius < 0"
    end if

  end subroutine sanity_checks

  pure subroutine get_cone_extents(self, x_vel, y_vel, sound_speed, origin, radius)
    !< Given a velocity and sound speed, determine the extents of the mach cone

    class(cone_t), intent(in) :: self
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

  subroutine get_reference_state(self, reconstructed_state)
    !< Find the reference state of the cone only based on the cells that are touched by the cone.
    !< Sometimes when a cone is near a large density jump, a skewed reference state can cause
    !< negative densities and pressures. This aims to alleviate that problem...

    class(cone_t), intent(inout) :: self
    real(rk), dimension(:, :), intent(in) :: reconstructed_state !< [rho, u, v, a]

    integer(ik) :: ref_state_cell
    real(rk) :: ave_p
    real(rk) :: mach_number, sound_speed
    integer(ik) :: i

    self%n_neighbor_cells = size(reconstructed_state, dim=2)

    allocate(self%recon_state, mold=reconstructed_state)
    self%recon_state = reconstructed_state

    associate(rho=>reconstructed_state(1, :), &
              u=>reconstructed_state(2, :), &
              v=>reconstructed_state(3, :), &
              p=>reconstructed_state(4, :), &
              n=>self%n_neighbor_cells)

      do i = 1, n
        sound_speed = eos%sound_speed(pressure=p(i), density=rho(i))
        mach_number = sqrt(u(i)**2 + v(i)**2) / sound_speed

        if(mach_number >= 1.0_rk) then
          self%cell_is_supersonic(i) = .true.
        else
          self%cell_is_supersonic(i) = .false.
        end if
        ! print*, i, 'self%cell_is_supersonic(i)', self%cell_is_supersonic(i)

      end do

      ! The cone is transonic if there is a combo of super/subsonic cells
      if(count(self%cell_is_supersonic(1:n)) == n .or. &
         count(self%cell_is_supersonic(1:n)) == 0) then
        self%cone_is_transonic = .false.
      else
        self%cone_is_transonic = .true.
      end if

      ! if (self%cone_is_transonic) then
      !   self%reference_density = maxval(rho)
      !   self%reference_u = maxval(u)
      !   self%reference_v = maxval(v)
      !   ave_p = maxval(p)
      !   self%reference_sound_speed = eos%sound_speed(pressure=ave_p, &
      !                                            density=self%reference_density)
      ! else
      self%reference_density = sum(rho) / n
      self%reference_u = sum(u) / n
      self%reference_v = sum(v) / n
      ave_p = sum(p) / n
      self%reference_sound_speed = eos%sound_speed(pressure=ave_p, &
                                                   density=self%reference_density)
      ! end if
    end associate

    ! Tiny velocity fix
    if(abs(self%reference_u) < TINY_DIST) self%reference_u = 0.0_rk
    if(abs(self%reference_v) < TINY_DIST) self%reference_v = 0.0_rk

    self%reference_mach_number = sqrt(self%reference_u**2 + self%reference_v**2) / &
                                 self%reference_sound_speed
  end subroutine get_reference_state

  subroutine get_transonic_cone_extents(self, origin, radius)
    !< If the set of neighbor cells are transonic (some supersonic, some subsonic),
    !< then the cone needs to be extended so that it encorporates both supersonic and subsonic cells.

    class(cone_t), intent(inout) :: self
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

  subroutine find_p_prime_cell(self)
    class(cone_t), intent(inout) :: self
    real(rk) :: u_tilde, v_tilde

    u_tilde = self%reference_u
    v_tilde = self%reference_v

    select case(trim(self%cone_location))
    case('left/right midpoint')
    case('down/up midpoint')
    case('corner')
    case default
      error stop "Unknown cone_location in cone_t%find_p_prime_cell()"
    end select

  end subroutine find_p_prime_cell

  subroutine write_cone(self, unit, iotype, v_list, iostat, iomsg)
    !< Implementation of `write(*,*) vector_t`

    class(cone_t), intent(in) :: self     !< cone class
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
    write(unit, '(a, 4(i0, 1x),a)', iostat=iostat) '# Valid arcs in each cell: [', self%n_arcs, ']'
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

  pure logical function determine_if_p_prime_is_in_cell(edge_vectors, p_prime_vector) result(in_cell)
    !< Implementation of whether the P' point is inside the current cell/control volume. This
    !< uses the cross product of 2 vectors in 2d, which gives a scalar
    type(vector_t), intent(in) :: p_prime_vector
    real(rk), dimension(2, 2, 2), intent(in) :: edge_vectors

    type(vector_t) :: edge_vector_1, edge_vector_2

    in_cell = .false.
    edge_vector_1 = vector_t(x=[edge_vectors(1, 1, 1), edge_vectors(1, 2, 1)], &
                             y=[edge_vectors(2, 1, 1), edge_vectors(2, 2, 1)])

    edge_vector_2 = vector_t(x=[edge_vectors(1, 1, 2), edge_vectors(1, 2, 2)], &
                             y=[edge_vectors(2, 1, 2), edge_vectors(2, 2, 2)])

    ! In the text, these inequalities are swapped, but I define the P' position vector
    ! as starting at P0 and going to P'
    ! write(*,*) 'p_prime_vector: ', p_prime_vector
    ! write(*,*) 'edge_vector_1:  ', edge_vector_1
    ! write(*,*) 'edge_vector_2:  ', edge_vector_2

    ! In the text is has >= instead of <= for some reason, but the following works
    ! like it's supposed to.
    if((p_prime_vector.cross.edge_vector_1) <= 0.0_rk .and. &
       (edge_vector_2.cross.p_prime_vector) <= 0.0_rk) then
      in_cell = .true.
    else
      in_cell = .false.
    end if
  end function determine_if_p_prime_is_in_cell

end module mod_cone
