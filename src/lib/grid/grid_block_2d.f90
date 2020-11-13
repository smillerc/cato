module mod_grid_block_2d

  use iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use mod_grid_block, only: grid_block_t
  use mod_input, only: input_t
  use mod_parallel, only: tile_indices
  use h5fortran, only: hdf5_file, hsize_t
  use mod_quad_cell, only: quad_cell_t
  use mod_error, only: error_msg
  use mod_nondimensionalization, only: set_length_scale, l_0
  use mod_units, only: um_to_cm

  implicit none(type, external)

  private
  public :: grid_block_2d_t, new_2d_grid_block

  type, extends(grid_block_t) :: grid_block_2d_t
    real(rk), dimension(:, :), allocatable :: volume !< (i, j); volume of each cell
    real(rk), dimension(:, :), allocatable :: dx     !< (i, j); dx spacing of each cell
    real(rk), dimension(:, :), allocatable :: dy     !< (i, j); dy spacing of each cell
    real(rk), dimension(:, :), allocatable :: node_x      !< (i, j); x location of each node
    real(rk), dimension(:, :), allocatable :: node_y      !< (i, j); y location of each node
    real(rk), dimension(:, :), allocatable :: centroid_x  !< (i, j); x location of the cell centroid
    real(rk), dimension(:, :), allocatable :: centroid_y  !< (i, j); y location of the cell centroid
    real(rk), dimension(:, :, :), allocatable :: edge_lengths         !< ((edge_1:edge_n), i, j); length of each edge
    real(rk), dimension(:, :, :, :), allocatable :: edge_norm_vectors !< ((x,y), edge, i, j); normal direction vector of each face
  contains
    ! Private methods
    private
    procedure :: populate_element_specifications
    procedure :: scale_and_nondimensionalize

    ! Public methods
    procedure, public :: initialize => init_2d_block
    procedure, public :: gather
    procedure, public :: read_from_h5
    procedure, public :: write_to_h5

    ! Finalize
    final :: finalize_2d_block
  end type grid_block_2d_t

contains
  function new_2d_grid_block(input) result(grid)
    type(grid_block_2d_t), pointer :: grid
    class(input_t), intent(in) :: input
    allocate(grid)
    call grid%initialize(input)
  end function new_2d_grid_block

  subroutine init_2d_block(self, input)
    class(grid_block_2d_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    type(hdf5_file) :: h5
    integer(hsize_t), allocatable :: cell_shape(:)
    integer(ik) :: global_cell_dims(2) = 0
    integer(ik) :: global_node_dims(2) = 0
    integer(ik), dimension(2) :: lbounds = 0
    integer(ik), dimension(2) :: ubounds = 0
    integer(ik) :: indices(4)
    integer(ik) :: alloc_status

    call h5%initialize(trim(input%initial_condition_file), status='old', action='r')
    call h5%shape('/density', cell_shape)
    call h5%read('/n_ghost_layers', self%n_halo_cells)
    call h5%finalize()

    ! The grid includes the ghost layer information. We want to tile based on
    ! the real domain b/c it makes the book-keeping a bit easier.
    ! There are two sets of ghost/halo/boundary cells for each direction

    ! Get the total # of cells in each direction
    global_cell_dims = int(cell_shape, ik)

    ! Remove the halo cell count for now, b/c this makes partitioning easier
    global_cell_dims = global_cell_dims - (2 * self%n_halo_cells)

    self%global_dims(1:2) = global_cell_dims
    self%lbounds_global = [1, 1, 0]
    self%ubounds_global = [self%global_dims(1), self%global_dims(2), 0]

    ! Now partition via the tile_indices function (this splits based on current image number)
    indices = tile_indices(global_cell_dims)
    lbounds = indices([1, 3])
    ubounds = indices([2, 4])
    self%lbounds(1:2) = lbounds
    self%ubounds(1:2) = ubounds

    self%lbounds_halo(1:2) = self%lbounds(1:2) - self%n_halo_cells
    self%ubounds_halo(1:2) = self%ubounds(1:2) + self%n_halo_cells

    self%host_image_id = this_image()
    self%block_dims(1:2) = self%ubounds(1:2) - self%lbounds(1:2) + 1
    self%total_cells = size(self%block_dims)

    ! Determine if this grid block has an edge on one of the global boundaries
    if(self%lbounds(1) == 1) self%on_ilo_bc = .true.
    if(self%lbounds(2) == 1) self%on_jlo_bc = .true.
    if(self%ubounds(1) == self%ubounds_global(1)) self%on_ihi_bc = .true.
    if(self%ubounds(2) == self%ubounds_global(2)) self%on_jhi_bc = .true.

    ! Allocate the node-based arrays (thus the + 1 in the ilo/ihi)
    associate(ilo => self%lbounds_halo(1), ihi => self%ubounds_halo(1) + 1, &
              jlo => self%lbounds_halo(2), jhi => self%ubounds_halo(2) + 1)
      allocate(self%node_x(ilo:ihi, jlo:jhi))  !< (i, j); x location of each node
      allocate(self%node_y(ilo:ihi, jlo:jhi))  !< (i, j); y location of each node
    end associate

    ! Allocate all of the cell-based arrays
    associate(ilo => self%lbounds_halo(1), ihi => self%ubounds_halo(1), &
              jlo => self%lbounds_halo(2), jhi => self%ubounds_halo(2))
      allocate(self%volume(ilo:ihi, jlo:jhi))             !< (i, j); volume of each cell
      allocate(self%dx(ilo:ihi, jlo:jhi))                 !< (i, j); dx spacing of each cell
      allocate(self%dy(ilo:ihi, jlo:jhi))                 !< (i, j); dy spacing of each cell
      allocate(self%centroid_x(ilo:ihi, jlo:jhi))              !< (i, j); x location of the cell centroid
      allocate(self%centroid_y(ilo:ihi, jlo:jhi))              !< (i, j); y location of the cell centroid
      allocate(self%edge_lengths(4, ilo:ihi, jlo:jhi))         !< ((edge_1:edge_n), i, j); length of each edge
      allocate(self%edge_norm_vectors(2, 4, ilo:ihi, jlo:jhi)) !< ((x,y), edge, i, j);
    end associate

    self%volume = 0.0_rk
    self%dx = 0.0_rk
    self%dy = 0.0_rk
    self%node_x = 0.0_rk
    self%node_y = 0.0_rk
    self%centroid_x = 0.0_rk
    self%centroid_y = 0.0_rk
    self%edge_lengths = 0.0_rk
    self%edge_norm_vectors = 0.0_rk

    call self%read_from_h5(input)
    call self%populate_element_specifications()
    call self%scale_and_nondimensionalize()
  end subroutine init_2d_block

  subroutine print_grid_stats(self)
    class(grid_block_2d_t), intent(in) :: self
    integer :: ni, nj
    if(this_image() == 1) then
      write(*, '(a)') "Grid stats:"
      write(*, '(a)') "========================="
      write(*, '(a)') "Blocks"
      write(*, '(a)') "-------------"
      write(*, '(a, i0)') "Number of blocks          : ", num_images()
      write(*, '(2(a,i0))') "Average block size (cells): ", self%block_dims(1), " x ", self%block_dims(2)
      write(*, '(a, i6)') "Number of halo cells      : ", self%n_halo_cells
      write(*, '(a)') "Global"
      write(*, '(a)') "-------------"
      write(*, '(a, i0)') "Number of i cells: ", self%global_dims(1)
      write(*, '(a, i0)') "Number of j cells: ", self%global_dims(2)
      write(*, '(a, i0)') "Total cells      : ", size(self%global_dims)
    end if
  end subroutine print_grid_stats

  subroutine finalize_2d_block(self)
    type(grid_block_2d_t), intent(inout) :: self

    if(allocated(self%node_x)) deallocate(self%node_x)
    if(allocated(self%node_y)) deallocate(self%node_y)
    if(allocated(self%volume)) deallocate(self%volume)
    if(allocated(self%dx)) deallocate(self%dx)
    if(allocated(self%dy)) deallocate(self%dy)
    if(allocated(self%centroid_x)) deallocate(self%centroid_x)
    if(allocated(self%centroid_y)) deallocate(self%centroid_y)
    if(allocated(self%edge_lengths)) deallocate(self%edge_lengths)
    if(allocated(self%edge_norm_vectors)) deallocate(self%edge_norm_vectors)
  end subroutine

  ! --------------------------------------------------------------------
  ! I/O for HDF5
  ! --------------------------------------------------------------------
  subroutine read_from_h5(self, input)
    !< Read in the data from an hdf5 file. This will read from the
    !< file and only grab the per-image data, e.g. each grid block will
    !< only read from the indices it has been assigned with respect
    !< to the global domain.
    class(grid_block_2d_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    type(hdf5_file) :: h5
    integer(ik) :: alloc_status, ilo, ihi, jlo, jhi
    logical :: file_exists
    character(:), allocatable :: filename
    character(32) :: str_buff = ''
    character(300) :: msg = ''
    real(rk), allocatable, dimension(:, :) :: x

    integer(hsize_t), allocatable :: dims(:)

    if(input%restart_from_file) then
      filename = trim(input%restart_file)
    else
      filename = trim(input%initial_condition_file)
    end if

    file_exists = .false.
    inquire(file=filename, exist=file_exists)

    if(.not. file_exists) then
      call error_msg(module_name='mod_grid_block_2d', class_name='grid_block_2d_t', &
                     procedure_name='read_from_h5', &
                     message='Error in regular_2d_grid_t%initialize_from_hdf5(); file not found: "'//filename//'"', &
                     file_name=__FILE__, line_number=__LINE__)
    end if

    call h5%initialize(filename=trim(filename), status='old', action='r')

    call h5%read('/n_ghost_layers', self%n_halo_cells)
    if(self%n_halo_cells /= input%n_ghost_layers) then
      write(msg, '(2(a,i0),a)') "The number of ghost layers in the .hdf5 file (", self%n_halo_cells, ") does not match the"// &
        " input requirement set by the edge interpolation scheme (", input%n_ghost_layers, ")"

      call error_msg(module_name='mod_grid_block_2d', class_name='grid_block_2d_t', &
                     procedure_name='read_from_h5', &
                     message=msg, &
                     file_name=__FILE__, line_number=__LINE__)
    end if

    ! I couldn't get the slice read to work, so we read the whole array
    ! in and extract the slice later
    call h5%shape('/x', dims)
    allocate(x(dims(1), dims(2)))

    call h5%read(dname='/x', value=x)
    self%node_x = x(self%lbounds_halo(1) + self%n_halo_cells:self%ubounds_halo(1) + self%n_halo_cells + 1, &
                    self%lbounds_halo(2) + self%n_halo_cells:self%ubounds_halo(2) + self%n_halo_cells + 1)

    call h5%read(dname='/y', value=x)
    self%node_y = x(self%lbounds_halo(1) + self%n_halo_cells:self%ubounds_halo(1) + self%n_halo_cells + 1, &
                    self%lbounds_halo(2) + self%n_halo_cells:self%ubounds_halo(2) + self%n_halo_cells + 1)

    if(input%restart_from_file) then
      call h5%readattr('/x', 'units', str_buff)
      select case(trim(str_buff))
      case('um')
        self%node_x = self%node_x * um_to_cm
        self%node_y = self%node_y * um_to_cm
      case('cm')
        ! Do nothing, since cm is what the code works in
      case default
        error stop "Unknown x,y units in grid from .h5 file. Acceptable units are 'um' or 'cm'."
      end select
    end if

    call h5%finalize()
  end subroutine read_from_h5

  subroutine write_to_h5(self, filename, dataset)
    !< Write the global data to an hdf5 file. This gathers all to image 1 and
    !< writes to a single file. Technically, this could be done in parallel
    !< in the future, but this is the simple case
    class(grid_block_2d_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dataset
    type(hdf5_file) :: h5

    ! Gather all to the master image and write to file
    call h5%initialize(filename='grid.h5', status='new', action='w', comp_lvl=6)
    call h5%write(dname='/x', value=self%gather(var='x', image=1))
    call h5%write(dname='/y', value=self%gather(var='y', image=1))
    call h5%write(dname='/volume', value=self%gather(var='volume', image=1))

    call h5%finalize()
  end subroutine write_to_h5

  function gather(self, var, image)
    !< Performs a gather of field data to image.
    class(grid_block_2d_t), intent(in) :: self
    integer(ik), intent(in) :: image
    character(len=*), intent(in) :: var
    real(rk), allocatable, dimension(:, :) :: gather
    real(rk), allocatable :: gather_coarray(:, :)[:]
    integer(ik) :: ni, nj

    ! This will have halo regions write over themselves, but this shouldn't be a problem b/c
    ! they are the same
    select case(var)
    case('x', 'y')

      ni = self%global_dims(1) + 1 ! the +1 is b/c of nodes vs cells
      nj = self%global_dims(2) + 1 ! the +1 is b/c of nodes vs cells
      allocate(gather_coarray(ni, nj)[*])
      ! Allocate the node-based arrays (thus the + 1 in the ilo/ihi)
      associate(ilo => self%lbounds(1), ihi => self%ubounds(1) + 1, &
                jlo => self%lbounds(2), jhi => self%ubounds(2) + 1)

        select case(var)
        case('x')
          gather_coarray(ilo:ihi, jlo:jhi)[image] = self%node_x(ilo:ihi, jlo:jhi)
        case('y')
          gather_coarray(ilo:ihi, jlo:jhi)[image] = self%node_y(ilo:ihi, jlo:jhi)
        end select

        sync all
        if(this_image() == image) gather = gather_coarray
      end associate

    case('volume')

      ni = self%global_dims(1)
      nj = self%global_dims(2)
      allocate(gather_coarray(ni, nj)[*])
      associate(ilo => self%lbounds(1), ihi => self%ubounds(1), &
                jlo => self%lbounds(2), jhi => self%ubounds(2))
        gather_coarray(ilo:ihi, jlo:jhi)[image] = self%volume(ilo:ihi, jlo:jhi)
        sync all
        if(this_image() == image) gather = gather_coarray
      end associate
    end select

    if (allocated(gather_coarray)) deallocate(gather_coarray)
  end function gather

  subroutine populate_element_specifications(self)
    !< Summary: Fill the element arrays up with the geometric information
    !< This seemed to be better for memory access patterns elsewhere in the code. Fortran prefers
    !< and structure of arrays rather than an array of structures

    class(grid_block_2d_t), intent(inout) :: self
    type(quad_cell_t) :: quad
    real(rk), dimension(4) :: x_coords
    real(rk), dimension(4) :: y_coords

    real(rk), dimension(8) :: p_x !< x coords of the cell corners and midpoints (c1,m1,c2,m2,c3,m3,c4,m4)
    real(rk), dimension(8) :: p_y !< x coords of the cell corners and midpoints (c1,m1,c2,m2,c3,m3,c4,m4)

    integer(ik) :: i, j

    x_coords = 0.0_rk
    y_coords = 0.0_rk

    do j = self%lbounds_halo(2), self%ubounds_halo(2)
      do i = self%lbounds_halo(1), self%ubounds_halo(1)
        associate(x => self%node_x, y => self%node_y)
          x_coords = [x(i, j), x(i + 1, j), x(i + 1, j + 1), x(i, j + 1)]
          y_coords = [y(i, j), y(i + 1, j), y(i + 1, j + 1), y(i, j + 1)]

          call quad%initialize(x_coords, y_coords)

        end associate

        self%volume(i, j) = quad%volume
        self%centroid_x(i, j) = quad%centroid(1)
        self%centroid_y(i, j) = quad%centroid(2)
        self%edge_lengths(:, i, j) = quad%edge_lengths

        ! call quad%get_cell_point_coords(p_x, p_y)
        ! self%node_x(:, i, j) = p_x
        ! self%node_y(:, i, j) = p_y

        self%edge_norm_vectors(:, :, i, j) = quad%edge_norm_vectors
        self%dx(i, j) = quad%min_dx
        self%dy(i, j) = quad%min_dy
      end do
    end do

  end subroutine populate_element_specifications

  subroutine scale_and_nondimensionalize(self)
    !< Scale the grid so that the cells are of size close to 1. If the grid is uniform,
    !< then everything (edge length and volume) are all 1. If not uniform, then the smallest
    !< edge legnth is 1. The scaling is done via the smallest edge length. This also sets
    !< the length scale for the non-dimensionalization module

    class(grid_block_2d_t), intent(inout) :: self

    real(rk) :: diff
    real(rk), save :: min_edge_length[*]
    real(rk), save :: max_edge_length[*]

    min_edge_length = minval(self%edge_lengths)
    max_edge_length = maxval(self%edge_lengths)

    ! Now broadcast the global max/min to all the images
    call co_max(max_edge_length)
    call co_min(min_edge_length)

    if(min_edge_length < tiny(1.0_rk)) error stop "Error in grid initialization, the cell min_edge_length = 0"

    diff = max_edge_length - min_edge_length
    if(diff < 2.0_rk * epsilon(1.0_rk)) then
      self%is_uniform = .true.
    end if

    ! Set l_0 for the entire code
    call set_length_scale(length_scale=min_edge_length)

    ! Scale so that the minimum edge length is 1
    self%node_x = self%node_x / min_edge_length
    self%node_y = self%node_y / min_edge_length
    self%edge_lengths = self%edge_lengths / min_edge_length
    self%centroid_x = self%centroid_x / min_edge_length
    self%centroid_y = self%centroid_y / min_edge_length
    self%node_x = self%node_x / min_edge_length
    self%node_y = self%node_y / min_edge_length
    self%dx = self%dx / min_edge_length
    self%dy = self%dy / min_edge_length

    ! If the grid is uniform, then we can make it all difinitively 1
    if(self%is_uniform) then
      if(this_image() == 1) then
        write(*, '(a)') "The grid is uniform, setting volume and edge lengths to 1, now that everything is scaled"
      endif
      self%volume = 1.0_rk
      ! self%cell_edge_lengths = 1.0_rk
      ! self%dx = 1.0_rk
      ! self%dy = 1.0_rk
      ! self%min_dx = 1.0_rk
      ! self%max_dx = 1.0_rk
    else
      self%volume = self%volume / min_edge_length**2
      ! self%min_dx = minval(self%node_x(lbound(self%node_x, 1) + 1:ubound(self%node_x, 1), :) - &
      !                     self%node_x(lbound(self%node_x, 1):ubound(self%node_x, 1) - 1, :))
      ! self%max_dx = maxval(self%node_x(lbound(self%node_x, 1) + 1:ubound(self%node_x, 1), :) - &
      !                     self%node_x(lbound(self%node_x, 1):ubound(self%node_x, 1) - 1, :))

      ! self%min_dy = minval(self%node_y(:, lbound(self%node_y, 2) + 1:ubound(self%node_y, 2)) - &
      !                     self%node_y(:, lbound(self%node_y, 2):ubound(self%node_y, 2) - 1))
      ! self%max_dy = maxval(self%node_y(:, lbound(self%node_y, 2) + 1:ubound(self%node_y, 2)) - &
      !                     self%node_y(:, lbound(self%node_y, 2):ubound(self%node_y, 2) - 1))
    end if

    ! self%xmin = minval(self%node_x)
    ! self%xmax = maxval(self%node_x)
    ! self%ymin = minval(self%node_y)
    ! self%ymax = maxval(self%node_y)

    ! self%x_length = abs(self%xmax - self%xmin)
    ! if(self%x_length <= 0) error stop "grid%x_length <= 0"

    ! self%y_length = abs(self%ymax - self%ymin)
    ! if(self%y_length <= 0) error stop "grid%x_length <= 0"

  end subroutine scale_and_nondimensionalize
end module
