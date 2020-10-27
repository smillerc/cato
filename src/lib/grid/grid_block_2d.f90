module mod_grid_block_2d

  use iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use mod_grid_block, only: grid_block_t
  use mod_input, only: input_t
  use mod_parallel, only: tile_indices
  use h5fortran, only: hdf5_file, hsize_t

  implicit none(type, external)

  private
  public :: grid_block_2d_t, new_2d_grid_block

  type, extends(grid_block_t) :: grid_block_2d_t
    real(rk), dimension(:, :), allocatable :: cell_volume !< (i, j); volume of each cell
    real(rk), dimension(:, :), allocatable :: cell_dx     !< (i, j); dx spacing of each cell
    real(rk), dimension(:, :), allocatable :: cell_dy     !< (i, j); dy spacing of each cell
    real(rk), dimension(:, :), allocatable :: node_x      !< (i, j); x location of each node
    real(rk), dimension(:, :), allocatable :: node_y      !< (i, j); y location of each node
    real(rk), dimension(:, :), allocatable :: centroid_x  !< (i, j); x location of the cell centroid
    real(rk), dimension(:, :), allocatable :: centroid_y  !< (i, j); y location of the cell centroid
    real(rk), dimension(:, :, :), allocatable :: edge_lengths         !< ((edge_1:edge_n), i, j); length of each edge
    real(rk), dimension(:, :, :, :), allocatable :: edge_norm_vectors !< ((x,y), edge, i, j); normal direction vector of each face
  contains
    procedure :: initialize => init_2d_block
    procedure :: gather
    procedure :: read_from_h5
    procedure :: write_to_h5
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
    integer(ik) :: global_node_dims(2)
    integer(hsize_t), allocatable :: node_shape(:)
    integer(ik) :: indices(4)
    integer(ik) :: alloc_status

    call h5%initialize(trim(input%initial_condition_file), status='old', action='r')
    call h5%shape('/node_x', node_shape)
    call h5%finalize()

    allocate(self%global_node_dims(2)); self%global_node_dims = 0
    allocate(self%global_cell_dims(2)); self%global_cell_dims = 0
    allocate(self%domain_node_shape(2)); self%domain_node_shape = 0
    allocate(self%domain_cell_shape(2)); self%domain_cell_shape = 0
    allocate(self%node_lbounds(2)); self%node_lbounds = 0
    allocate(self%node_ubounds(2)); self%node_ubounds = 0
    allocate(self%node_lbounds_halo(2)); self%node_lbounds_halo = 0
    allocate(self%node_ubounds_halo(2)); self%node_ubounds_halo = 0
    allocate(self%cell_lbounds(2)); self%cell_lbounds = 0
    allocate(self%cell_ubounds(2)); self%cell_ubounds = 0
    allocate(self%cell_lbounds_halo(2)); self%cell_lbounds_halo = 0
    allocate(self%cell_ubounds_halo(2)); self%cell_ubounds_halo = 0
    allocate(self%on_bc(2,2)); self%on_bc = .false.

    ! The shape from the hdf5 file must use hsize_t integers, this just converts so that
    ! the tile_indicies function works (which expects int32 integers)
    global_node_dims(1) = int(node_shape(1), ik)
    global_node_dims(2) = int(node_shape(2), ik)

    ! Coarray domain-decomposition. The indices are tiled based on the number of available images
    indices = tile_indices(global_node_dims)

    self%global_node_dims = global_node_dims
    self%global_cell_dims = global_node_dims - 1

    self%node_lbounds = indices([1, 3])
    self%node_ubounds = indices([2, 4])

    self%cell_lbounds = self%node_lbounds
    self%cell_ubounds = self%node_ubounds - 1

    ! The domain shape is exent in each direction
    self%domain_node_shape(1) = self%node_ubounds(1) - self%node_lbounds(1) + 1 ! # of nodes in the i direction
    self%domain_node_shape(2) = self%node_ubounds(2) - self%node_lbounds(2) + 1 ! # of nodes in the j direction

    self%domain_cell_shape(1) = self%cell_ubounds(1) - self%cell_lbounds(1) + 1 ! # of cells in the i direction
    self%domain_cell_shape(2) = self%cell_ubounds(2) - self%cell_lbounds(2) + 1 ! # of cells in the j direction
    self%host_image_id = this_image()

    ! Determine if this field (or subdomain) has an edge on one of the boundaries
    ! When data is synced between images / subdomains, those on the boundary must
    ! have the procedure done via the boundary condition classes. Those on the interior
    ! of the domain just exchange info with their neighbors
    if(self%node_lbounds(1) == 1) self%on_bc(1, 1) = .true.
    if(self%node_lbounds(2) == 1) self%on_bc(1, 2) = .true.
    if(self%node_ubounds(1) == self%global_node_dims(1)) self%on_bc(2, 1) = .true.
    if(self%node_ubounds(2) == self%global_node_dims(2)) self%on_bc(2, 2) = .true.

    self%cell_lbounds_halo = self%cell_lbounds - self%n_halo_cells
    self%cell_ubounds_halo = self%cell_ubounds + self%n_halo_cells

    self%node_lbounds_halo = self%node_lbounds - self%n_halo_cells
    self%node_ubounds_halo = self%node_ubounds + self%n_halo_cells

    ! Allocate the node-based arrays
    associate(ilo => self%node_lbounds(1), ihi => self%node_ubounds(1), &
              jlo => self%node_lbounds(2), jhi => self%node_ubounds(2))
      allocate(self%node_x(ilo:ihi, jlo:jhi))  !< (i, j); x location of each node
      allocate(self%node_y(ilo:ihi, jlo:jhi))  !< (i, j); y location of each node
    end associate

    ! Allocate all of the cell-based arrays
    associate(ilo => self%cell_lbounds(1), ihi => self%cell_ubounds(1), &
              jlo => self%cell_lbounds(2), jhi => self%cell_ubounds(2))
      allocate(self%cell_volume(ilo:ihi, jlo:jhi))             !< (i, j); volume of each cell
      allocate(self%cell_dx(ilo:ihi, jlo:jhi))                 !< (i, j); dx spacing of each cell
      allocate(self%cell_dy(ilo:ihi, jlo:jhi))                 !< (i, j); dy spacing of each cell
      allocate(self%centroid_x(ilo:ihi, jlo:jhi))              !< (i, j); x location of the cell centroid
      allocate(self%centroid_y(ilo:ihi, jlo:jhi))              !< (i, j); y location of the cell centroid
      allocate(self%edge_lengths(4, ilo:ihi, jlo:jhi))         !< ((edge_1:edge_n), i, j); length of each edge
      allocate(self%edge_norm_vectors(2, 4, ilo:ihi, jlo:jhi)) !< ((x,y), edge, i, j);
    end associate

    self%cell_volume = 0.0_rk
    self%cell_dx = 0.0_rk
    self%cell_dy = 0.0_rk
    self%node_x = 0.0_rk
    self%node_y = 0.0_rk
    self%centroid_x = 0.0_rk
    self%centroid_y = 0.0_rk
    self%edge_lengths = 0.0_rk
    self%edge_norm_vectors = 0.0_rk
  end subroutine init_2d_block

  subroutine finalize_2d_block(self)
    type(grid_block_2d_t), intent(inout) :: self

    if(allocated(self%node_x)) deallocate(self%node_x)
    if(allocated(self%node_y)) deallocate(self%node_y)
    if(allocated(self%cell_volume)) deallocate(self%cell_volume)
    if(allocated(self%cell_dx)) deallocate(self%cell_dx)
    if(allocated(self%cell_dy)) deallocate(self%cell_dy)
    if(allocated(self%centroid_x)) deallocate(self%centroid_x)
    if(allocated(self%centroid_y)) deallocate(self%centroid_y)
    if(allocated(self%edge_lengths)) deallocate(self%edge_lengths)
    if(allocated(self%edge_norm_vectors)) deallocate(self%edge_norm_vectors)
    if(allocated(self%global_node_dims )) deallocate(self%global_node_dims )
    if(allocated(self%global_cell_dims )) deallocate(self%global_cell_dims )
    if(allocated(self%domain_node_shape)) deallocate(self%domain_node_shape)
    if(allocated(self%domain_cell_shape)) deallocate(self%domain_cell_shape)
    if(allocated(self%node_lbounds     )) deallocate(self%node_lbounds     )
    if(allocated(self%node_ubounds     )) deallocate(self%node_ubounds     )
    if(allocated(self%node_lbounds_halo)) deallocate(self%node_lbounds_halo)
    if(allocated(self%node_ubounds_halo)) deallocate(self%node_ubounds_halo)
    if(allocated(self%cell_lbounds     )) deallocate(self%cell_lbounds     )
    if(allocated(self%cell_ubounds     )) deallocate(self%cell_ubounds     )
    if(allocated(self%cell_lbounds_halo)) deallocate(self%cell_lbounds_halo)
    if(allocated(self%cell_ubounds_halo)) deallocate(self%cell_ubounds_halo)
    if(allocated(self%on_bc)) deallocate(self%on_bc)
    
  end subroutine

  ! --------------------------------------------------------------------
  ! I/O for HDF5
  ! --------------------------------------------------------------------
  subroutine read_from_h5(self, filename)
    !< Read in the data from an hdf5 file. This will read from the
    !< file and only grab the per-image data, e.g. each grid block will
    !< only read from the indices it has been assigned with respect
    !< to the global domain.
    class(grid_block_2d_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    type(hdf5_file) :: h5

    call h5%initialize(filename=trim(filename), status='old', action='r')

    ! TODO: what about reading from a 0-based hdf5 file?
    ! Read on a per-image basis. The -1 is b/c hdf5 files are C-based
    ! and therefore start at 0 instead of 1
    associate(ilo => self%node_lbounds(1) - 1, ihi => self%node_ubounds(1) - 1, &
              jlo => self%node_lbounds(2) - 1, jhi => self%node_ubounds(2) - 1)

      call h5%read(dname='x', &
                   value=self%node_x(ilo:ihi, jlo:jhi), &
                   istart=[ilo, jhi], iend=[jlo, jhi])
      call h5%read(dname='y', &
                   value=self%node_y(ilo:ihi, jlo:jhi), &
                   istart=[ilo, jhi], iend=[jlo, jhi])
    end associate
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

    select case(var)
    case('x', 'y')
      allocate(gather_coarray(self%global_node_dims(1), self%global_node_dims(2))[*])
      associate(is => self%node_lbounds(1), ie => self%node_ubounds(1), &
                js => self%node_lbounds(2), je => self%node_ubounds(2))

        select case(var)
        case('x')
          gather_coarray(is:ie, js:je)[image] = self%node_x(is:ie, js:je)
        case('y')
          gather_coarray(is:ie, js:je)[image] = self%node_y(is:ie, js:je)
        end select

        sync all
        if(this_image() == image) gather = gather_coarray
      end associate

    case('volume')
      allocate(gather_coarray(self%global_cell_dims(1), self%global_cell_dims(2))[*])
      associate(is => self%cell_lbounds(1), ie => self%cell_ubounds(1), &
                js => self%cell_lbounds(2), je => self%cell_ubounds(2))
        gather_coarray(is:ie, js:je)[image] = self%cell_volume(is:ie, js:je)
        sync all
        if(this_image() == image) gather = gather_coarray
      end associate
    end select

    deallocate(gather_coarray)
  end function gather

end module
