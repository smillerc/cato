! MIT License
! Copyright (c) 2019 Sam Miller
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

module mod_regular_2d_grid
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_out => output_unit, std_error => error_unit
  use mod_grid, only: grid_t, C1, M1, C2, M2, C3, M3, C4, M4
  use mod_units
  use mod_quad_cell, only: quad_cell_t
  use mod_input, only: input_t
  use hdf5_interface, only: hdf5_file
  use mod_globals, only: debug_print
  use mod_floating_point_utils, only: equal
  use mod_nondimensionalization, only: set_length_scale, l_0
  use mod_functional, only: arange

  implicit none

  private
  public :: regular_2d_grid_t, new_regular_2d_grid

  type, extends(grid_t) :: regular_2d_grid_t
    !< Summary: The regular_2d_grid_t type holds all of the geometry info relevant to the grid.

    ! real(rk), dimension(:, :, :, :, :), allocatable :: cell_edge_vectors
    !< (x:y, loc1:loc2, face1:face4, i, j)
    !< This describes the point locations of the set of edge vectors at each location the mach cone is evaluated.
    !< For instance, for edge E1, a mach cone is constructed at C1, M1, and C2.
    !< At the corner location 1 at C1, e.g. index `loc1`, or k=1,1 (in paper lingo), the 3 points that define the 2
    !< edge vectors are C2, C1, and C4 or (C2,C1) and (C4,C1). To take advantage of memory access patterns, since
    !< mach cones are made 12x at each cell for each timestep, these points are lumped into a big array. The index
    !< loc1:loc2 spans the 1st corner, midpoint, and 2nd corner.

    ! k=1,1 -> [C2, C1, C4]
    ! k=1,c -> [C2, M1, C1]
    ! k=1,2 -> [C3, C2, C1]

    ! k=2,1 -> [C3, C2, C1]
    ! k=2,c -> [C3, M2, C2]
    ! k=2,2 -> [C4, C3, C2]

    ! k=3,1 -> [C4, C3, C2]
    ! k=3,c -> [C4, M3, C3]
    ! k=3,2 -> [C1, C4, C3]

    ! k=4,1 -> [C1, C4, C3]
    ! k=4,c -> [C1, M4, C4]
    ! k=4,2 -> [C2, C1, C4]
  contains
    procedure, public :: initialize
    procedure, private :: initialize_from_hdf5
    procedure, private :: initialize_from_ini
    procedure, private :: populate_element_specifications
    procedure, private :: scale_and_nondimensionalize
    procedure, public :: get_x
    procedure, public :: get_y
    procedure, public :: get_cell_volumes
    procedure, public :: get_cell_edge_lengths
    procedure, public :: get_cell_edge_norm_vectors
    procedure, public :: get_midpoint_persistent_vectors
    procedure, public :: get_corner_persistent_vectors
    procedure, public :: finalize
    procedure, public :: copy
    procedure, public :: print_grid_stats
    final :: force_finalization

  end type regular_2d_grid_t

  interface new_regular_2d_grid
    module procedure :: constructor
  end interface
contains

  function constructor(input) result(grid)
    type(regular_2d_grid_t), pointer :: grid
    class(input_t), intent(in) :: input
    allocate(grid)
    call grid%initialize(input)
  end function constructor

  subroutine initialize(self, input)
    !< Implementation of the grid initialization process for a regular 2d grid
    class(regular_2d_grid_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    integer(ik) :: i, j, alloc_status

    ! Low node/cell indices (always starts at 1)
    self%ilo_node = 1
    self%jlo_node = 1
    self%ilo_cell = 1
    self%jlo_cell = 1

    if(input%read_init_cond_from_file .or. input%restart_from_file) then
      write(*, '(a)') 'Initializing the grid via .hdf5'
      call self%initialize_from_hdf5(input)
    else
      write(*, '(a)') 'Initializing the grid via .ini'
      call self%initialize_from_ini(input)
    end if

    ! Allocate cell based arrays
    associate(imin=>self%ilo_bc_cell, imax=>self%ihi_bc_cell, &
              jmin=>self%jlo_bc_cell, jmax=>self%jhi_bc_cell)

      allocate(self%cell_volume(imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%cell_volume"
      self%cell_volume = 0.0_rk

      allocate(self%cell_dx(imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%cell_dx"
      self%cell_dx = 0.0_rk

      allocate(self%cell_dy(imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%cell_dy"
      self%cell_dy = 0.0_rk

      allocate(self%cell_centroid_x(imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%cell_centroid_x"
      self%cell_centroid_x = 0.0_rk

      allocate(self%cell_centroid_y(imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%cell_centroid_y"
      self%cell_centroid_y = 0.0_rk

      allocate(self%cell_edge_lengths(4, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%cell_edge_lengths"
      self%cell_edge_lengths = 0.0_rk

      allocate(self%cell_edge_norm_vectors(2, 4, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%cell_edge_norm_vectors"
      self%cell_edge_norm_vectors = 0.0_rk

      allocate(self%cell_node_x(1:8, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%cell_node_x"
      self%cell_node_x = 0.0_rk

      allocate(self%cell_node_y(1:8, imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%cell_node_y"
      self%cell_node_y = 0.0_rk
    end associate

    ! Allocate the edge vector arrays
    associate(imin_node=>self%ilo_node, imax_node=>self%ihi_node, &
              jmin_node=>self%jlo_node, jmax_node=>self%jhi_node, &
              imin_cell=>self%ilo_cell, imax_cell=>self%ihi_cell, &
              jmin_cell=>self%jlo_cell, jmax_cell=>self%jhi_cell)

      allocate(self%corner_edge_vectors(2, 0:4, imin_node:imax_node, jmin_node:jmax_node), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%corner_edge_vectors"

      allocate(self%corner_edge_vectors_scale(imin_node:imax_node, jmin_node:jmax_node), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%corner_edge_vectors_scale"

      allocate(self%downup_midpoint_edge_vectors(2, 0:2, imin_node:imax_node, jmin_cell:jmax_cell), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%downup_midpoint_edge_vectors"

      allocate(self%downup_midpoint_edge_vectors_scale(imin_node:imax_node, jmin_cell:jmax_cell), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%downup_midpoint_edge_vectors_scale"

      allocate(self%leftright_midpoint_edge_vectors(2, 0:2, imin_cell:imax_cell, jmin_node:jmax_node), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%leftright_midpoint_edge_vectors"

      allocate(self%leftright_midpoint_edge_vectors_scale(imin_cell:imax_cell, jmin_node:jmax_node), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%leftright_midpoint_edge_vectors_scale"
    end associate

    call self%populate_element_specifications()
    call self%scale_and_nondimensionalize()

    call self%get_corner_persistent_vectors(scale=.false., shift=.false.)
    call self%get_midpoint_persistent_vectors(edge='left', scale=.false., shift=.false.)
    call self%get_midpoint_persistent_vectors(edge='bottom', scale=.false., shift=.false.)

    call self%print_grid_stats()
  end subroutine initialize

  subroutine copy(out_grid, in_grid)
    class(grid_t), intent(in) :: in_grid
    class(regular_2d_grid_t), intent(inout) :: out_grid

    call debug_print('Running regular_2d_grid_t%copy()', __FILE__, __LINE__)

    out_grid%ilo_bc_node = in_grid%ilo_bc_node
    out_grid%jlo_bc_node = in_grid%jlo_bc_node
    out_grid%ihi_bc_node = in_grid%ihi_bc_node
    out_grid%jhi_bc_node = in_grid%jhi_bc_node
    out_grid%ilo_node = in_grid%ilo_node
    out_grid%jlo_node = in_grid%jlo_node
    out_grid%ihi_node = in_grid%ihi_node
    out_grid%jhi_node = in_grid%jhi_node
    out_grid%ni_node = in_grid%ni_node
    out_grid%nj_node = in_grid%nj_node
    out_grid%ilo_bc_cell = in_grid%ilo_bc_cell
    out_grid%jlo_bc_cell = in_grid%jlo_bc_cell
    out_grid%ihi_bc_cell = in_grid%ihi_bc_cell
    out_grid%jhi_bc_cell = in_grid%jhi_bc_cell
    out_grid%ilo_cell = in_grid%ilo_cell
    out_grid%jlo_cell = in_grid%jlo_cell
    out_grid%ihi_cell = in_grid%ihi_cell
    out_grid%jhi_cell = in_grid%jhi_cell
    out_grid%ni_cell = in_grid%ni_cell
    out_grid%nj_cell = in_grid%nj_cell
    out_grid%xmin = in_grid%xmin
    out_grid%xmax = in_grid%xmax
    out_grid%min_dx = in_grid%min_dx
    out_grid%max_dx = in_grid%max_dx
    out_grid%ymin = in_grid%ymin
    out_grid%ymax = in_grid%ymax
    out_grid%min_dy = in_grid%min_dy
    out_grid%max_dy = in_grid%max_dy
    out_grid%x_length = in_grid%x_length
    out_grid%y_length = in_grid%y_length

    ! allocate(out_grid%node_x, source=in_grid%node_x)
    out_grid%node_x = in_grid%node_x

    ! if(allocated(out_grid%node_y)) deallocate(out_grid%node_y)
    ! allocate(out_grid%node_y, source=in_grid%node_y)
    out_grid%node_y = in_grid%node_y

    ! if(allocated(out_grid%cell_volume)) deallocate(out_grid%cell_volume)
    ! allocate(out_grid%cell_volume, source=in_grid%cell_volume)
    out_grid%cell_volume = in_grid%cell_volume

    ! if(allocated(out_grid%cell_centroid_xy)) deallocate(out_grid%cell_centroid_xy)
    ! allocate(out_grid%cell_centroid_xy, source=in_grid%cell_centroid_xy)
    out_grid%cell_centroid_x = in_grid%cell_centroid_x
    out_grid%cell_centroid_y = in_grid%cell_centroid_y

    ! if(allocated(out_grid%cell_edge_lengths)) deallocate(out_grid%cell_edge_lengths)
    ! allocate(out_grid%cell_edge_lengths, source=in_grid%cell_edge_lengths)
    out_grid%cell_edge_lengths = in_grid%cell_edge_lengths

    ! if(allocated(out_grid%cell_node_xy)) deallocate(out_grid%cell_node_xy)
    ! allocate(out_grid%cell_node_xy, source=in_grid%cell_node_xy)
    out_grid%cell_node_x = in_grid%cell_node_x
    out_grid%cell_node_y = in_grid%cell_node_y

    ! if(allocated(out_grid%cell_edge_norm_vectors)) deallocate(out_grid%cell_edge_norm_vectors)
    ! allocate(out_grid%cell_edge_norm_vectors, source=in_grid%cell_edge_norm_vectors)
    out_grid%cell_edge_norm_vectors = in_grid%cell_edge_norm_vectors
  end subroutine copy

  subroutine print_grid_stats(self)
    class(regular_2d_grid_t), intent(inout) :: self

    print *
    write(*, '(a)') "Grid stats:"
    write(*, '(a)') "========================="
    write(*, '(a, i6)') "n_ghost_layers: ", self%n_ghost_layers
    write(*, '(a, i6, a, i6)') "i nodes: ", self%ilo_node, ' -> ', self%ihi_node
    write(*, '(a, i6, a, i6)') "j nodes: ", self%jlo_node, ' -> ', self%jhi_node
    write(*, '(a, i6, a, i6)') "i cells: ", self%ilo_cell, ' -> ', self%ihi_cell
    write(*, '(a, i6, a, i6)') "j cells: ", self%jlo_cell, ' -> ', self%jhi_cell
    write(*, *)
    write(*, '(a)') "Ghost Regions"
    write(*, '(a)') "-------------"
    write(*, '(a, i6, a, i6)') "i nodes: ", self%ilo_bc_node, ' & ', self%ihi_bc_node
    write(*, '(a, i6, a, i6)') "j nodes: ", self%jlo_bc_node, ' & ', self%jhi_bc_node
    write(*, '(a, i6, a, i6)') "i cells: ", self%ilo_bc_cell, ' & ', self%ihi_bc_cell
    write(*, '(a, i6, a, i6)') "j cells: ", self%jlo_bc_cell, ' & ', self%jhi_bc_cell
    write(*, *)
    write(*, '(a)') "Totals"
    write(*, '(a)') "------"
    write(*, '(a, i0)') "ni_nodes: ", self%ni_node
    write(*, '(a, i0)') "nj_nodes: ", self%nj_node
    write(*, '(a, i0)') "ni_cells: ", self%ni_cell
    write(*, '(a, i0)') "nj_cells: ", self%nj_cell
    write(*, '(a, i0)') "total cells: ", self%nj_cell * self%ni_node

    write(*, *)
    write(*, '(a)') "Extents"
    write(*, '(a)') "-------"
    write(*, '(2(a, es10.3))') "x range [non-dim] (w/o ghost)", &
      self%node_x(self%ilo_node, self%jlo_node), '  -> ', self%node_x(self%ihi_node, self%jhi_node)
    write(*, '(2(a, es10.3))') "y range [non-dim] (w/o ghost)", &
      self%node_x(self%ilo_node, self%jlo_node), '  -> ', self%node_y(self%ihi_node, self%jhi_node)
    write(*, '(2(a, es10.3))') "x range [dim]     (w/o ghost)", &
      self%node_x(self%ilo_node, self%jlo_node) * l_0, '  -> ', self%node_x(self%ihi_node, self%jhi_node) * l_0
    write(*, '(2(a, es10.3))') "y range [dim]     (w/o ghost)", &
      self%node_x(self%ilo_node, self%jlo_node) * l_0, '  -> ', self%node_y(self%ihi_node, self%jhi_node) * l_0
    write(*, '(a)') "========================="
    write(*, *)

  end subroutine print_grid_stats

  subroutine initialize_from_hdf5(self, input)
    !< Initialize the grid from an .hdf5 file. This allows flexibility on the shape
    !< and contents of the grid (structured, but flexible quadrilateral cell shapes)
    !< Note: Initial grid files MUST include the ghost layer
    class(regular_2d_grid_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    type(hdf5_file) :: h5
    integer(ik) :: alloc_status
    logical :: file_exists
    character(:), allocatable :: filename
    character(32) :: str_buff = ''

    real(rk), dimension(:, :), allocatable :: x
    real(rk), dimension(:, :), allocatable :: y

    if(input%restart_from_file) then
      filename = trim(input%restart_file)
    else
      filename = trim(input%initial_condition_file)
    end if

    file_exists = .false.
    inquire(file=filename, exist=file_exists)

    if(.not. file_exists) then
      write(*, '(a)') 'Error in regular_2d_grid_t%initialize_from_hdf5(); file not found: "'//filename//'"'
      error stop 'Error in regular_2d_grid_t%initialize_from_hdf5(); file not found, exiting...'
    end if

    call h5%initialize(filename=filename, status='old', action='r')

    ! Low i/j boundary condition indices
    call h5%get('/n_ghost_layers', self%n_ghost_layers)

    if(self%n_ghost_layers /= input%n_ghost_layers) then
      write(std_error, '(2(a, i0))') "regular_2d_grid_t%n_ghost_layers: ", &
        self%n_ghost_layers, ", input%n_ghost_layers: ", input%n_ghost_layers
      error stop "The number of ghost layers in the .hdf5 file does not match the"// &
        " input requirement set by the edge interpolation scheme"
    end if

    call h5%get('/x', x)
    call h5%get('/y', y)

    if(input%restart_from_file) then
      call h5%readattr('/x', 'units', str_buff)
      select case(trim(str_buff))
      case('um')
        x = x * um_to_cm
        y = y * um_to_cm
      case('cm')
        ! Do nothing, since cm is what the code works in
      case default
        error stop "Unknown x,y units in grid from .h5 file. Acceptable units are 'um' or 'cm'."
      end select
    end if

    call h5%finalize()

    ! Node
    self%ilo_node = 1
    self%jlo_node = 1
    self%ilo_bc_node = 1 - self%n_ghost_layers
    self%jlo_bc_node = 1 - self%n_ghost_layers
    self%ihi_bc_node = ubound(x, dim=1) - self%n_ghost_layers
    self%jhi_bc_node = ubound(x, dim=2) - self%n_ghost_layers
    self%ihi_node = self%ihi_bc_node - self%n_ghost_layers
    self%jhi_node = self%jhi_bc_node - self%n_ghost_layers

    ! Cell
    self%ilo_cell = 1
    self%jlo_cell = 1
    self%ilo_bc_cell = self%ilo_bc_node
    self%jlo_bc_cell = self%jlo_bc_node

    self%ihi_cell = self%ihi_node - 1
    self%jhi_cell = self%jhi_node - 1
    self%ihi_bc_cell = self%ihi_bc_node - 1
    self%jhi_bc_cell = self%jhi_bc_node - 1

    write(*, '(a, 2(i5))') "cell jlo: ", self%jlo_cell, self%jlo_bc_cell
    write(*, '(a, 2(i5))') "cell jhi: ", self%jhi_cell, self%jhi_bc_cell
    write(*, '(a, 2(i5))') "node jlo: ", self%jlo_node, self%jlo_bc_node
    write(*, '(a, 2(i5))') "node jhi: ", self%jhi_node, self%jhi_bc_node

    self%ni_node = size(x, dim=1) - (2 * self%n_ghost_layers)
    self%nj_node = size(y, dim=2) - (2 * self%n_ghost_layers)
    self%ni_cell = self%ni_node - 1
    self%nj_cell = self%nj_node - 1

    ! Allocate node based arrays
    associate(imin=>self%ilo_bc_node, imax=>self%ihi_bc_node, &
              jmin=>self%jlo_bc_node, jmax=>self%jhi_bc_node)

      ! error stop
      allocate(self%node_x(imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%node_x from .ini input file"

      allocate(self%node_y(imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%node_y from .ini input file"
      self%node_x(imin:imax, jmin:jmax) = x!(1:imax,1:jmax)
      self%node_y(imin:imax, jmin:jmax) = y!(1:imax,1:jmax)
    end associate

    deallocate(x)
    deallocate(y)
  end subroutine initialize_from_hdf5

  subroutine initialize_from_ini(self, input)
    !< Initialize the nodes from the .ini input file -> requires equal spacing
    class(regular_2d_grid_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    integer(ik) :: alloc_status
    integer(ik) :: i, j

    self%n_ghost_layers = input%n_ghost_layers

    self%ilo_bc_node = 1 - self%n_ghost_layers
    self%jlo_bc_node = 1 - self%n_ghost_layers
    self%ilo_bc_cell = 1 - self%n_ghost_layers
    self%jlo_bc_cell = 1 - self%n_ghost_layers

    ! High node/cell indices
    self%ihi_node = input%ni_nodes
    self%jhi_node = input%nj_nodes
    self%ihi_cell = self%ihi_node - 1
    self%jhi_cell = self%jhi_node - 1

    ! High i/j boundary condition indices
    self%ihi_bc_node = self%ihi_node + self%n_ghost_layers
    self%jhi_bc_node = self%jhi_node + self%n_ghost_layers
    self%ihi_bc_cell = self%ihi_cell + self%n_ghost_layers
    self%jhi_bc_cell = self%jhi_cell + self%n_ghost_layers

    self%ni_node = input%ni_nodes
    self%nj_node = input%nj_nodes
    self%ni_cell = self%ni_node - 1
    self%nj_cell = self%nj_node - 1

    ! Allocate node based arrays
    associate(imin=>self%ilo_bc_node, imax=>self%ihi_bc_node, &
              jmin=>self%jlo_bc_node, jmax=>self%jhi_bc_node)

      allocate(self%node_x(imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%node_x from .ini input file"

      allocate(self%node_y(imin:imax, jmin:jmax), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%node_y from .ini input file"
    end associate

    self%xmin = input%xmin / l_0
    self%xmax = input%xmax / l_0
    self%ymin = input%ymin / l_0
    self%ymax = input%ymax / l_0

    self%x_length = abs(self%xmax - self%xmin)
    if(self%x_length <= 0) error stop "grid%x_length <= 0"

    self%y_length = abs(self%ymax - self%ymin)
    if(self%y_length <= 0) error stop "grid%x_length <= 0"

    self%min_dx = self%x_length / real(self%ni_cell, rk)
    self%max_dx = self%min_dx ! placeholder for now
    if(self%min_dx <= 0) error stop "grid%dx <= 0"

    self%min_dy = self%y_length / real(self%nj_cell, rk)
    self%max_dy = self%min_dy ! placeholder for now
    if(self%min_dy <= 0) error stop "grid%dy <= 0"

    do j = self%jlo_bc_node, self%jhi_bc_node
      self%node_x(:, j) = arange(start=self%xmin - self%n_ghost_layers * self%min_dx, &
                                 end=self%xmax + self%n_ghost_layers * self%min_dx, increment=self%min_dx)
    end do

    do i = self%ilo_bc_node, self%ihi_bc_node
      self%node_y(i, :) = arange(start=self%ymin - self%n_ghost_layers * self%min_dy, &
                                 end=self%ymax + self%n_ghost_layers * self%min_dy, increment=self%min_dy)
    end do

  end subroutine initialize_from_ini

  subroutine populate_element_specifications(self)
    !< Summary: Fill the element arrays up with the geometric information
    !< This seemed to be better for memory access patterns elsewhere in the code. Fortran prefers
    !< and structure of arrays rather than an array of structures

    class(regular_2d_grid_t), intent(inout) :: self
    type(quad_cell_t) :: quad
    real(rk), dimension(4) :: x_coords
    real(rk), dimension(4) :: y_coords

    real(rk), dimension(8) :: p_x !< x coords of the cell corners and midpoints (c1,m1,c2,m2,c3,m3,c4,m4)
    real(rk), dimension(8) :: p_y !< x coords of the cell corners and midpoints (c1,m1,c2,m2,c3,m3,c4,m4)

    real(rk) :: min_vol, max_vol
    integer(ik) :: i, j

    x_coords = 0.0_rk
    y_coords = 0.0_rk

    do j = self%jlo_bc_cell, self%jhi_bc_cell
      do i = self%ilo_bc_cell, self%ihi_bc_cell
        associate(x=>self%node_x, y=>self%node_y)
          x_coords = [x(i, j), x(i + 1, j), x(i + 1, j + 1), x(i, j + 1)]
          y_coords = [y(i, j), y(i + 1, j), y(i + 1, j + 1), y(i, j + 1)]

          call quad%initialize(x_coords, y_coords)

        end associate

        self%cell_volume(i, j) = quad%volume
        self%cell_centroid_x(i, j) = quad%centroid(1)
        self%cell_centroid_y(i, j) = quad%centroid(2)
        self%cell_edge_lengths(:, i, j) = quad%edge_lengths

        call quad%get_cell_point_coords(p_x, p_y)
        self%cell_node_x(:, i, j) = p_x
        self%cell_node_y(:, i, j) = p_y

        self%cell_edge_norm_vectors(:, :, i, j) = quad%edge_norm_vectors
        self%cell_dx(i, j) = quad%min_dx
        self%cell_dy(i, j) = quad%min_dy
      end do
    end do

  end subroutine populate_element_specifications

  subroutine scale_and_nondimensionalize(self)
    !< Scale the grid so that the cells are of size close to 1. If the grid is uniform,
    !< then everything (edge length and volume) are all 1. If not uniform, then the smallest
    !< edge legnth is 1. The scaling is done via the smallest edge length. This also sets
    !< the length scale for the non-dimensionalization module

    class(regular_2d_grid_t), intent(inout) :: self

    real(rk) :: diff, min_edge_length, max_edge_length

    min_edge_length = minval(self%cell_edge_lengths)
    max_edge_length = maxval(self%cell_edge_lengths)

    if(min_edge_length < tiny(1.0_rk)) error stop "Error in grid initialization, the cell min_edge_length = 0"

    diff = max_edge_length - min_edge_length
    if(diff < 2.0_rk * epsilon(1.0_rk)) then
      self%grid_is_uniform = .true.
    end if

    ! Set l_0 for the entire code
    call set_length_scale(length_scale=min_edge_length)

    ! Scale so that the minimum edge length is 1
    self%node_x = self%node_x / min_edge_length
    self%node_y = self%node_y / min_edge_length
    self%cell_edge_lengths = self%cell_edge_lengths / min_edge_length
    self%cell_centroid_x = self%cell_centroid_x / min_edge_length
    self%cell_centroid_y = self%cell_centroid_y / min_edge_length
    self%cell_node_x = self%cell_node_x / min_edge_length
    self%cell_node_y = self%cell_node_y / min_edge_length
    self%cell_dx = self%cell_dx / min_edge_length
    self%cell_dy = self%cell_dy / min_edge_length

    ! If the grid is uniform, then we can make it all difinitively 1
    if(self%grid_is_uniform) then
      write(*, '(a)') "The grid is uniform, setting volume and edge lengths to 1, now that everything is scaled"
      self%cell_volume = 1.0_rk
      self%cell_edge_lengths = 1.0_rk
      self%cell_dx = 1.0_rk
      self%cell_dy = 1.0_rk
      self%min_dx = 1.0_rk
      self%max_dx = 1.0_rk
    else
      self%cell_volume = self%cell_volume / min_edge_length**2
      self%min_dx = minval(self%node_x(lbound(self%node_x, 1) + 1:ubound(self%node_x, 1), :) - &
                           self%node_x(lbound(self%node_x, 1):ubound(self%node_x, 1) - 1, :))
      self%max_dx = maxval(self%node_x(lbound(self%node_x, 1) + 1:ubound(self%node_x, 1), :) - &
                           self%node_x(lbound(self%node_x, 1):ubound(self%node_x, 1) - 1, :))

      self%min_dy = minval(self%node_y(:, lbound(self%node_y, 2) + 1:ubound(self%node_y, 2)) - &
                           self%node_y(:, lbound(self%node_y, 2):ubound(self%node_y, 2) - 1))
      self%max_dy = maxval(self%node_y(:, lbound(self%node_y, 2) + 1:ubound(self%node_y, 2)) - &
                           self%node_y(:, lbound(self%node_y, 2):ubound(self%node_y, 2) - 1))
    end if

    self%xmin = minval(self%node_x)
    self%xmax = maxval(self%node_x)
    self%ymin = minval(self%node_y)
    self%ymax = maxval(self%node_y)

    self%x_length = abs(self%xmax - self%xmin)
    if(self%x_length <= 0) error stop "grid%x_length <= 0"

    self%y_length = abs(self%ymax - self%ymin)
    if(self%y_length <= 0) error stop "grid%x_length <= 0"

  end subroutine scale_and_nondimensionalize

  subroutine force_finalization(self)
    type(regular_2d_grid_t), intent(inout) :: self
    call self%finalize()
  end subroutine force_finalization

  subroutine finalize(self)
    class(regular_2d_grid_t), intent(inout) :: self

    call debug_print('Running regular_2d_grid_t%finalize()', __FILE__, __LINE__)

    if(allocated(self%cell_volume)) deallocate(self%cell_volume)
    if(allocated(self%node_x)) deallocate(self%node_x)
    if(allocated(self%node_y)) deallocate(self%node_y)
    if(allocated(self%cell_centroid_x)) deallocate(self%cell_centroid_x)
    if(allocated(self%cell_centroid_y)) deallocate(self%cell_centroid_y)
    if(allocated(self%cell_edge_lengths)) deallocate(self%cell_edge_lengths)
    if(allocated(self%cell_node_x)) deallocate(self%cell_node_x)
    if(allocated(self%cell_node_y)) deallocate(self%cell_node_y)
    if(allocated(self%cell_edge_norm_vectors)) deallocate(self%cell_edge_norm_vectors)
    if(allocated(self%corner_edge_vectors)) deallocate(self%corner_edge_vectors)
    if(allocated(self%downup_midpoint_edge_vectors)) deallocate(self%downup_midpoint_edge_vectors)
    if(allocated(self%leftright_midpoint_edge_vectors)) deallocate(self%leftright_midpoint_edge_vectors)
    if(allocated(self%cell_dx)) deallocate(self%cell_dx)
    if(allocated(self%cell_dy)) deallocate(self%cell_dy)

  end subroutine finalize

  pure function get_x(self, i, j) result(x)
    !< Public interface to get x
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik), intent(in) :: i, j
    real(rk) :: x
    x = self%node_x(i, j)
  end function get_x

  pure function get_y(self, i, j) result(y)
    !< Public interface to get y
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik), intent(in) :: i, j
    real(rk) :: y
    y = self%node_y(i, j)
  end function get_y

  pure function get_cell_volumes(self, i, j) result(cell_volume)
    !< Public interface to get cell_volume
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik), intent(in) :: i, j
    real(rk) :: cell_volume
    cell_volume = self%cell_volume(i, j)
  end function get_cell_volumes

  pure function get_cell_edge_lengths(self, i, j, f) result(cell_edge_lengths)
    !< Public interface to get cell_edge_lengths
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik), intent(in) :: i, j, f
    real(rk) :: cell_edge_lengths
    cell_edge_lengths = self%cell_edge_lengths(f, i, j)
  end function get_cell_edge_lengths

  pure function get_cell_edge_norm_vectors(self, i, j, f, xy) result(cell_edge_norm_vectors)
    !< Public interface to get cell_edge_norm_vectors
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik), intent(in) :: i, j, f, xy
    real(rk) :: cell_edge_norm_vectors
    cell_edge_norm_vectors = self%cell_edge_norm_vectors(xy, f, i, j)
  end function get_cell_edge_norm_vectors

  subroutine get_midpoint_persistent_vectors(self, edge, scale, shift)
    !< Public interface to get_midpoint_vectors
    class(regular_2d_grid_t), intent(inout) :: self
    character(len=*), intent(in) :: edge ! 'bottom', 'top', 'left', 'right'
    logical, intent(in) :: scale
    logical, intent(in) :: shift
    ! scale_and_shift will enable scaling the lengths of all the vectors to near 1. It will
    ! use the smallest length of the set and divide all the vectors by it. Shifting will move
    ! the entire set to the origin at (0,0).

    integer(ik) :: i, j, v

    real(rk), dimension(2) :: vector_length
    real(rk), dimension(2) :: origin

    ! Corner/midpoint index convention         Cell Indexing convention
    ! --------------------------------         ------------------------
    !
    !   C----M----C----M----C
    !   |         |         |                             E3
    !   O    x    O    x    O                      C4-----M3----C3
    !   |         |         |                      |            |
    !   C----M----C----M----C                  E4  M4     X     M2  E2
    !   |         |         |                      |            |
    !   O    x    O    x    O                      C1----M1----C2
    !   |         |         |                            E1
    !   C----M----C----M----C
    !
    ! For left/right midpoints, the edge vectors go left then right.
    ! The neighboring cells are above (i,j) and below (i,j-1)
    ! For quad cells, C - corner, M - midpoint, E - edge

    ! self%cell_node_xy indexing convention
    !< ((x,y), (point_1:point_n), (node=1,midpoint=2), i, j); The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints

    select case(trim(edge))
    case('left')
      ! Left edge of cell (i,j)
      !  |         C4         |
      !  |         |          |
      !  |         M4         |
      !  | (i-1,j) |  (i,j)   |
      !  |        C1          |
      do j = self%jlo_cell, self%jhi_cell
        do i = self%ilo_node, self%ihi_node

          origin(1) = self%cell_node_x(M4, i, j)
          origin(2) = self%cell_node_y(M4, i, j)
          self%downup_midpoint_edge_vectors(:, 0, i, j) = origin

          self%downup_midpoint_edge_vectors(1, 1, i, j) = self%cell_node_x(C1, i, j) ! vector 1 (C1)
          self%downup_midpoint_edge_vectors(2, 1, i, j) = self%cell_node_y(C1, i, j) ! vector 1 (C1)

          self%downup_midpoint_edge_vectors(1, 2, i, j) = self%cell_node_x(C4, i, j) ! vector 2 (C4)
          self%downup_midpoint_edge_vectors(2, 2, i, j) = self%cell_node_y(C4, i, j) ! vector 2 (C4)

          ! ! shift to origin
          ! if(shift) then
          !   do v = 1, 2
          !     self%downup_midpoint_edge_vectors(:, v, i, j) = self%downup_midpoint_edge_vectors(:, v, i, j) - origin
          !   end do
          ! end if

          ! ! get the length of each vector
          ! if(scale) then
          !   do v = 1, 2
          !     vector_length(v) = sqrt(self%downup_midpoint_edge_vectors(1, v, i, j)**2 + &
          !                             self%downup_midpoint_edge_vectors(2, v, i, j)**2)
          !   end do

          !   ! scale by the smallest length
          !   self%downup_midpoint_edge_vectors_scale(i, j) = minval(vector_length)
          !   self%downup_midpoint_edge_vectors(:, :, i, j) = self%downup_midpoint_edge_vectors(:, :, i, j) &
          !                                                   / minval(vector_length)
          ! end if

        end do
      end do

      if(.not. scale) then
        self%downup_midpoint_edge_vectors_scale = 1.0_rk
      end if

    case('bottom')
      ! Bottom edge of cell (i,j)
      !       |          |
      !       |  (i,j)   |
      !       C1---M1---C2
      !       | (i,j-1)  |
      !       |          |
      do j = self%jlo_node, self%jhi_node
        do i = self%ilo_cell, self%ihi_cell
          origin(1) = self%cell_node_x(M1, i, j) ! origin (M1)
          origin(2) = self%cell_node_y(M1, i, j) ! origin (M1)
          self%leftright_midpoint_edge_vectors(:, 0, i, j) = origin

          self%leftright_midpoint_edge_vectors(1, 1, i, j) = self%cell_node_x(C1, i, j) ! vector 1 (C1)
          self%leftright_midpoint_edge_vectors(2, 1, i, j) = self%cell_node_y(C1, i, j) ! vector 1 (C1)

          self%leftright_midpoint_edge_vectors(1, 2, i, j) = self%cell_node_x(C2, i, j) ! vector 2 (C2)
          self%leftright_midpoint_edge_vectors(2, 2, i, j) = self%cell_node_y(C2, i, j) ! vector 2 (C2)

          ! shift to origin
          if(shift) then
            do v = 1, 2
              self%leftright_midpoint_edge_vectors(:, v, i, j) = self%leftright_midpoint_edge_vectors(:, v, i, j) - origin
            end do
          end if

          ! get the length of each vector
          if(scale) then
            do v = 1, 2
              vector_length(v) = sqrt(self%leftright_midpoint_edge_vectors(1, v, i, j)**2 + &
                                      self%leftright_midpoint_edge_vectors(2, v, i, j)**2)
            end do

            ! scale by the smallest length
            self%leftright_midpoint_edge_vectors_scale(i, j) = minval(vector_length)
            self%leftright_midpoint_edge_vectors(:, :, i, j) = self%leftright_midpoint_edge_vectors(:, :, i, j) &
                                                               / minval(vector_length)
          end if
        end do
      end do

      if(.not. scale) then
        self%leftright_midpoint_edge_vectors_scale = 1.0_rk
      end if

    case default
      error stop "Invalid location for midpoint edge vector request"
    end select

  end subroutine get_midpoint_persistent_vectors

  subroutine get_corner_persistent_vectors(self, scale, shift)
    !< Public interface to get_corner_vectors
    class(regular_2d_grid_t), intent(inout) :: self

    ! scale_and_shift will enable scaling the lengths of all the vectors to near 1. It will
    ! use the smallest length of the set and divide all the vectors by it. Shifting will move
    ! the entire set to the origin at (0,0).

    logical, intent(in) :: scale
    logical, intent(in) :: shift
    real(rk), dimension(2) :: origin
    real(rk), dimension(4) :: vector_length

    integer(ik) :: i, j, v, N

    ! Corner/midpoint index convention         Cell Indexing convention
    ! --------------------------------         ------------------------
    !
    !   C----M----C----M----C
    !   |         |         |                             E3
    !   O    x    O    x    O                      N4-----M3----N3
    !   |         |         |                      |            |
    !   C----M----C----M----C                  E4  M4     C     M2  E2
    !   |         |         |                      |            |
    !   O    x    O    x    O                      N1----M1----N2
    !   |         |         |                            E1
    !   C----M----C----M----C
    !
    ! For left/right midpoints, the edge vectors go left then right.
    ! The neighboring cells are above (i,j) and below (i,j-1)
    ! For quad cells, C - corner, M - midpoint, E - edge

    ! Corner vector set
    !       cell 4                   cell 3
    !      (i-1,j)                   (i,j)
    !  C4----M3----C3    P3   C4----M3----C3
    !  |            |    |    |            |
    !  M4    X4    M2    |    M4    X3    M2
    !  |            |    |    |            |
    !  C1----M1----C2    |    C1----M1----C2
    !                    |
    !  P4----------------O-----------------P2
    !                    |
    !  C4----M3----C3    |    C4----M3----C3
    !  |            |    |    |            |
    !  M4    X1    M2    |    M4    X2    M2
    !  |            |    |    |            |
    !  C1----M1----C2    P1   C1----M1----C2
    !      cell 1                  cell 2
    !     (i-1,j-1)               (i,j-1)

    ! self%cell_node_xy indexing convention
    !< ((x,y), (point_1:point_n), (node=1,midpoint=2), i, j); The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints
    self%corner_edge_vectors = 0.0_rk
    do j = self%jlo_node, self%jhi_node
      do i = self%ilo_node, self%ihi_node
        origin(1) = self%cell_node_x(C1, i, j)
        origin(2) = self%cell_node_y(C1, i, j)
        self%corner_edge_vectors(:, 0, i, j) = origin

        self%corner_edge_vectors(1, 1, i, j) = self%cell_node_x(C2, i - 1, j - 1)  ! vector 1
        self%corner_edge_vectors(2, 1, i, j) = self%cell_node_y(C2, i - 1, j - 1)  ! vector 1

        self%corner_edge_vectors(1, 2, i, j) = self%cell_node_x(C2, i, j)          ! vector 2
        self%corner_edge_vectors(2, 2, i, j) = self%cell_node_y(C2, i, j)          ! vector 2

        self%corner_edge_vectors(1, 3, i, j) = self%cell_node_x(C4, i, j)          ! vector 3
        self%corner_edge_vectors(2, 3, i, j) = self%cell_node_y(C4, i, j)          ! vector 3

        self%corner_edge_vectors(1, 4, i, j) = self%cell_node_x(C4, i - 1, j - 1)  ! vector 4
        self%corner_edge_vectors(2, 4, i, j) = self%cell_node_y(C4, i - 1, j - 1)  ! vector 4

        ! shift to origin
        if(shift) then
          do v = 1, 4
            self%corner_edge_vectors(:, v, i, j) = self%corner_edge_vectors(:, v, i, j) - origin
          end do
        end if

        ! get the length of each vector
        if(scale) then
          do v = 1, 4
            vector_length(v) = sqrt(self%corner_edge_vectors(1, v, i, j)**2 + &
                                    self%corner_edge_vectors(2, v, i, j)**2)
          end do

          ! scale by the smallest length
          self%corner_edge_vectors_scale(i, j) = minval(vector_length)
          self%corner_edge_vectors(:, :, i, j) = self%corner_edge_vectors(:, :, i, j) / minval(vector_length)
        end if
      end do
    end do

    if(.not. scale) then
      self%corner_edge_vectors_scale = 1.0_rk
    end if
  end subroutine get_corner_persistent_vectors

end module mod_regular_2d_grid
