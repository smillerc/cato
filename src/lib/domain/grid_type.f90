module mod_grid

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t

  implicit none

  private
  public :: grid_t

  type, abstract :: grid_t
    ! private
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

    real(rk), dimension(:, :), allocatable :: cell_volume !< (i,j); volume of each cell

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

    real(rk), dimension(:, :, :), allocatable :: cell_centroid_xy
    !< (xy, i, j); (x,y) location of the cell centroid

    real(rk), dimension(:, :, :), allocatable :: cell_edge_lengths
    !< (i, j, face1:face4); length of each edge

    real(rk), dimension(:, :, :, :, :), allocatable :: cell_node_xy
    !< (xy, point, node:midpoint, i, j); The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints

    real(rk), dimension(:, :, :, :), allocatable :: cell_edge_norm_vectors
    !< (xy, face, i, j); normal direction vector of each face

  contains
    procedure(initialize), deferred :: initialize
    procedure(get_2d_data), deferred, public :: get_x
    procedure(get_2d_data), deferred, public :: get_y
    procedure(get_2d_data), deferred, public :: get_cell_volumes
    procedure(get_cell_centroid_xy), deferred, public :: get_cell_centroid_xy
    procedure(get_cell_edge_lengths), deferred, public :: get_cell_edge_lengths
    procedure(get_4d_data), deferred, public :: get_cell_edge_norm_vectors
    ! procedure(get_cell_node_xy), deferred, public :: get_cell_node_xy
  end type grid_t

  abstract interface
    subroutine initialize(self, input)
      import :: grid_t
      import :: input_t
      class(grid_t), intent(inout) :: self
      class(input_t), intent(in) :: input
    end subroutine

    pure function get_2d_data(self, i, j) result(x)
      !< Public interface to get get_2d_data
      import :: grid_t
      import :: ik, rk
      class(grid_t), intent(in) :: self
      integer(ik), intent(in) :: i, j
      real(rk) :: x
    end function

    pure function get_cell_centroid_xy(self, i, j) result(cell_centroid_xy)
      !< Public interface to get get_cell_centroid_xy
      import :: grid_t
      import :: ik, rk
      class(grid_t), intent(in) :: self
      integer(ik), intent(in) :: i, j
      real(rk), dimension(2) :: cell_centroid_xy
    end function

    pure function get_cell_edge_lengths(self, i, j, f) result(cell_centroids)
      !< Public interface to get get_cell_edge_lengths
      import :: grid_t
      import :: ik, rk
      class(grid_t), intent(in) :: self
      integer(ik), intent(in) :: i, j, f
      real(rk) :: cell_centroids
    end function

    pure function get_4d_data(self, i, j, f, xy) result(cell_edge_midpoints)
      !< Public interface to get get_4d_data
      import :: grid_t
      import :: ik, rk
      class(grid_t), intent(in) :: self
      integer(ik), intent(in) :: i, j, f, xy
      real(rk) :: cell_edge_midpoints
    end function

    pure function get_cell_node_xy(self, i, j, n, xy) result(cell_edge_norm_vectors)
      !< Public interface to get get_5d_data
      import :: grid_t
      import :: ik, rk
      class(grid_t), intent(in) :: self
      integer(ik), intent(in) :: i, j, n, xy
      real(rk) :: cell_edge_norm_vectors
    end function

  end interface

end module mod_grid
