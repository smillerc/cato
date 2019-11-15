module mod_regular_2d_grid
  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_grid, only: grid_t
  use mod_quad_cell, only: quad_cell_t
  use mod_input, only: input_t

  implicit none

  private
  public :: regular_2d_grid_t

  type, extends(grid_t) :: regular_2d_grid_t
    !< Summary: The regular_2d_grid_t type holds all of the geometry info relevant to the grid.
    private
    integer(ik) :: ihi !< i starting index
    integer(ik) :: ilo !< i ending index
    integer(ik) :: ni
    integer(ik) :: jlo !< j starting index
    integer(ik) :: jhi !< j ending index
    integer(ik) :: nj

    real(rk) :: xmin !< Min x location
    real(rk) :: xmax !< Max x location
    real(rk) :: dx
    real(rk) :: ymin !< Min y location
    real(rk) :: ymax !< Max y location
    real(rk) :: dy

    real(rk) :: x_length !< Length of the domain in x
    real(rk) :: y_length !< Length of the domain in y

    ! Rather than keep an array of element types, make arrays that hold the
    ! element information. Use the element type to populate this
    ! once at the beginning of the simulation or on demand

    real(rk), dimension(:, :), allocatable :: cell_volumes !< (i,j); volume of each cell

    !  Numbering convention for a 2D quadrilateral cell
    !
    !                     E3
    !                     M3
    !              N4-----o-----N3
    !              |            |
    !      E4   M4 o      C     o M2   E2
    !              |            |
    !              N1-----o-----N2
    !                     M1
    !                     E1
    ! N: node or vertex
    ! E: face or edge or interface
    ! M: midpoint of the edge (o)
    ! C: cell or control volume (in finite-volume lingo)

    real(rk), dimension(:, :), allocatable :: node_x !< (i, j, N1:N4); x location of each node
    real(rk), dimension(:, :), allocatable :: node_y !< (i, j, N1:N4); y location of each node

    ! real(rk), dimension(:, :, :, :), allocatable :: cell_node_xy
    !           !< (i, j, N1:M4, x:y), so it stores each node (N) and midpoint (M) in order,
    !           !< e.g. the node dimension is [N1, M1, N2, ... , N4, M4]

    real(rk), dimension(:, :, :), allocatable :: cell_centroid_xy
    !< (i, j, x:y); x,y of the cell centroid

    real(rk), dimension(:, :, :), allocatable :: cell_edge_lengths
    !< (i, j, face1:face4); length of each edge

    real(rk), dimension(:, :, :, :), allocatable :: cell_edge_norm_vectors
    !< (i, j, face1:face4, x:y); normal direction vector of each face

    real(rk), dimension(:, :, :, :, :), allocatable :: cell_edge_vectors
    !< (i, j, face1:face4, loc1:loc2, x:y)
    !< This describes the point locations of the set of edge vectors at each location the mach cone is evaluated.
    !< For instance, for edge E1, a mach cone is constructed at N1, M1, and N2.
    !< At the corner location 1 at N1, e.g. index `loc1`, or k=1,1 (in paper lingo), the 3 points that define the 2
    !< edge vectors are N2, N1, and N4 or (N2,N1) and (N4,N1). To take advantage of memory access patterns, since
    !< mach cones are made 12x at each cell for each timestep, these points are lumped into a big array. The index
    !< loc1:loc2 spans the 1st corner, midpoint, and 2nd corner.

    ! k=1,1 -> [N2, N1, N4]
    ! k=1,c -> [N2, M1, N1]
    ! k=1,2 -> [N3, N2, N1]

    ! k=2,1 -> [N3, N2, N1]
    ! k=2,c -> [N3, M2, N2]
    ! k=2,2 -> [N4, N3, N2]

    ! k=3,1 -> [N4, N3, N2]
    ! k=3,c -> [N4, M3, N3]
    ! k=3,2 -> [N1, N4, N3]

    ! k=4,1 -> [N1, N4, N3]
    ! k=4,c -> [N1, M4, N4]
    ! k=4,2 -> [N2, N1, N4]
  contains
    procedure, public :: initialize
    procedure, private :: populate_element_specifications
    procedure, public :: get_ihi
    procedure, public :: get_ilo
    procedure, public :: get_jlo
    procedure, public :: get_jhi
    procedure, public :: get_ni
    procedure, public :: get_nj
    procedure, public :: get_xmin
    procedure, public :: get_xmax
    procedure, public :: get_ymin
    procedure, public :: get_ymax
    procedure, public :: get_x_length
    procedure, public :: get_y_length
    procedure, public :: get_x
    procedure, public :: get_y
    procedure, public :: get_dx
    procedure, public :: get_dy
    procedure, public :: get_cell_volumes
    procedure, public :: get_cell_centroids
    procedure, public :: get_cell_edge_lengths
    procedure, public :: get_cell_node_xy
    procedure, public :: get_cell_edge_norm_vectors
    procedure, public :: finalize
    final :: force_finalization

  end type regular_2d_grid_t

contains

  subroutine initialize(self, input)

    class(regular_2d_grid_t), intent(inout) :: self
    class(input_t), intent(in) :: input

    integer(ik) :: alloc_status, i, j

    self%ilo = 1
    self%ihi = self%ilo + input%ni - 1
    self%ni = input%ni

    self%jlo = 1
    self%jhi = self%jlo + input%nj - 1
    self%nj = input%nj

    self%xmin = input%xmin
    self%xmax = input%xmax
    self%ymin = input%ymin
    self%ymax = input%ymax
    self%x_length = self%xmax - self%xmin
    if(self%x_length <= 0) error stop "grid%x_length <= 0"

    self%y_length = self%ymax - self%ymin
    if(self%y_length <= 0) error stop "grid%x_length <= 0"

    self%dx = self%x_length / (self%ni - 1)
    if(self%dx <= 0) error stop "grid%dx <= 0"

    self%dy = self%y_length / (self%nj - 1)
    if(self%dy <= 0) error stop "grid%dy <= 0"

    allocate(self%node_x(self%ilo:self%ihi, self%jlo:self%jhi), &
             stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%x"

    do i = self%ilo, self%ihi
      self%node_x(i, :) = self%xmin + (i - 1) * self%dx
    end do

    allocate(self%node_y(self%ilo:self%ihi, self%jlo:self%jhi), &
             stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%y"

    do j = self%jlo, self%jhi
      self%node_y(:, j) = self%ymin + (j - 1) * self%dy
    end do

    allocate(self%cell_volumes(self%ilo:self%ihi - 1, self%jlo:self%jhi - 1), &
             stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%cell_volumes"
    self%cell_volumes = 0.0_rk

    allocate(self%cell_centroid_xy(self%ilo:self%ihi - 1, self%jlo:self%jhi - 1, 2), &
             stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%cell_centroids"
    self%cell_centroid_xy = 0.0_rk

    allocate(self%cell_edge_lengths(self%ilo:self%ihi - 1, self%jlo:self%jhi - 1, 4), &
             stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%cell_edge_lengths"
    self%cell_edge_lengths = 0.0_rk

    allocate(self%cell_node_xy(self%ilo:self%ihi - 1, self%jlo:self%jhi - 1, 8, 2), &
             stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%cell_edge_midpoints"
    self%cell_node_xy = 0.0_rk

    allocate(self%cell_edge_norm_vectors(self%ilo:self%ihi - 1, self%jlo:self%jhi - 1, 4, 2), &
             stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate regular_2d_grid_t%cell_edge_norm_vectors"
    self%cell_edge_norm_vectors = 0.0_rk

    call self%populate_element_specifications()
  end subroutine

  subroutine populate_element_specifications(self)
    !< Summary: Fill the element arrays up with the geometric information
    !< This seemed to be better for memory access patterns elsewhere in the code. Fortran prefers
    !< and structure of arrays rather than an array of structures

    class(regular_2d_grid_t), intent(inout) :: self
    class(quad_cell_t), allocatable :: quad

    integer(ik) :: i, j

    do j = self%jlo, self%jhi - 1
      do i = self%ilo, self%ihi - 1

        allocate(quad_cell_t :: quad)

        associate(x=>self%node_x, y=>self%node_y)
          call quad%initialize(x_coords=[x(i, j), x(i + 1, j), x(i + 1, j + 1), x(i, j + 1)], &
                               y_coords=[y(i, j), y(i + 1, j), y(i + 1, j + 1), y(i, j + 1)])

        end associate

        self%cell_volumes(i, j) = quad%volume
        ! self%cell_centroid_xy(i, j, :) = quad%centroid
        ! self%cell_edge_lengths(i, j, :) = quad%edge_lengths
        ! self%cell_edge_midpoints(i, j, :, :) = quad%edge_midpoints
        ! self%cell_edge_norm_vectors(i, j, :, :) = quad%edge_norm_vectors

        deallocate(quad)

      end do
    end do

  end subroutine

  subroutine force_finalization(self)
    type(regular_2d_grid_t), intent(inout) :: self
    call self%finalize
  end subroutine

  subroutine finalize(self)
    class(regular_2d_grid_t), intent(inout) :: self
    integer(ik) :: alloc_status

    print *, 'Finalizing regular_2d_grid_t'
    if(allocated(self%node_x)) deallocate(self%node_x, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate regular_2d_grid_t%x"

    if(allocated(self%node_y)) deallocate(self%node_y, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate regular_2d_grid_t%y"

    if(allocated(self%cell_volumes)) deallocate(self%cell_volumes, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate regular_2d_grid_t%cell_volumes"

    if(allocated(self%cell_centroid_xy)) deallocate(self%cell_centroid_xy, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate regular_2d_grid_t%cell_centroids"

    if(allocated(self%cell_edge_lengths)) deallocate(self%cell_edge_lengths, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate regular_2d_grid_t%cell_edge_lengths"

    if(allocated(self%cell_node_xy)) deallocate(self%cell_node_xy, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate regular_2d_grid_t%cell_node_xy"

    if(allocated(self%cell_edge_norm_vectors)) deallocate(self%cell_edge_norm_vectors, stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to deallocate regular_2d_grid_t%cell_edge_norm_vectors"

    print *, 'Done'
  end subroutine finalize

  pure function get_ihi(self) result(ihi)
    !< Public interface to get ihi
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik) :: ihi
    ihi = self%ihi
  end function

  pure function get_ilo(self) result(ilo)
    !< Public interface to get ilo
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik) :: ilo
    ilo = self%ilo
  end function

  pure function get_ni(self) result(ni)
    !< Public interface to get ni
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik) :: ni
    ni = self%ni
  end function

  pure function get_nj(self) result(nj)
    !< Public interface to get nj
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik) :: nj
    nj = self%nj
  end function

  pure function get_jlo(self) result(jlo)
    !< Public interface to get jlo
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik) :: jlo
    jlo = self%jlo
  end function

  pure function get_jhi(self) result(jhi)
    !< Public interface to get jhi
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik) :: jhi
    jhi = self%jhi
  end function

  pure function get_xmin(self) result(xmin)
    !< Public interface to get xmin
    class(regular_2d_grid_t), intent(in) :: self
    real(rk) :: xmin
    xmin = self%xmin
  end function

  pure function get_xmax(self) result(xmax)
    !< Public interface to get xmax
    class(regular_2d_grid_t), intent(in) :: self
    real(rk) :: xmax
    xmax = self%xmax
  end function

  pure function get_ymin(self) result(ymin)
    !< Public interface to get ymin
    class(regular_2d_grid_t), intent(in) :: self
    real(rk) :: ymin
    ymin = self%ymin
  end function

  pure function get_ymax(self) result(ymax)
    !< Public interface to get ymax
    class(regular_2d_grid_t), intent(in) :: self
    real(rk) :: ymax
    ymax = self%ymax
  end function

  pure function get_dx(self) result(dx)
    !< Public interface to get ymax
    class(regular_2d_grid_t), intent(in) :: self
    real(rk) :: dx
    dx = self%dx
  end function

  pure function get_dy(self) result(dy)
    !< Public interface to get ymax
    class(regular_2d_grid_t), intent(in) :: self
    real(rk) :: dy
    dy = self%dy
  end function

  pure function get_x_length(self) result(x_length)
    !< Public interface to get x_length
    class(regular_2d_grid_t), intent(in) :: self
    real(rk) :: x_length
    x_length = self%x_length
  end function

  pure function get_y_length(self) result(y_length)
    !< Public interface to get y_length
    class(regular_2d_grid_t), intent(in) :: self
    real(rk) :: y_length
    y_length = self%y_length
  end function

  pure function get_x(self, i, j) result(x)
    !< Public interface to get x
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik), intent(in) :: i, j
    real(rk) :: x
    x = self%node_x(i, j)
  end function

  pure function get_y(self, i, j) result(y)
    !< Public interface to get y
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik), intent(in) :: i, j
    real(rk) :: y
    y = self%node_y(i, j)
  end function

  pure function get_cell_volumes(self, i, j) result(cell_volumes)
    !< Public interface to get cell_volumes
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik), intent(in) :: i, j
    real(rk) :: cell_volumes
    cell_volumes = self%cell_volumes(i, j)
  end function

  pure function get_cell_centroids(self, i, j, xy) result(cell_centroids)
    !< Public interface to get cell_centroids
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik), intent(in) :: i, j, xy
    real(rk) :: cell_centroids
    cell_centroids = self%cell_centroid_xy(i, j, xy)
  end function

  pure function get_cell_edge_lengths(self, i, j, f) result(cell_edge_lengths)
    !< Public interface to get cell_edge_lengths
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik), intent(in) :: i, j, f
    real(rk) :: cell_edge_lengths
    cell_edge_lengths = self%cell_edge_lengths(i, j, f)
  end function

  pure function get_cell_node_xy(self, i, j, n, xy) result(cell_node_xy)
    !< Public interface to get cell_edge_midpoints
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik), intent(in) :: i, j, n, xy
    real(rk) :: cell_node_xy
    cell_node_xy = self%cell_node_xy(i, j, n, xy)
  end function

  pure function get_node_xy_pair(self, i, j, n_id) result(xy_pair)
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik), intent(in) :: i, j, n_id
    real(rk), dimension(2) :: xy_pair

    xy_pair = self
  end function

  pure function get_cell_edge_norm_vectors(self, i, j, f, xy) result(cell_edge_norm_vectors)
    !< Public interface to get cell_edge_norm_vectors
    class(regular_2d_grid_t), intent(in) :: self
    integer(ik), intent(in) :: i, j, f, xy
    real(rk) :: cell_edge_norm_vectors
    cell_edge_norm_vectors = self%cell_edge_norm_vectors(i, j, f, xy)
  end function

end module mod_regular_2d_grid
