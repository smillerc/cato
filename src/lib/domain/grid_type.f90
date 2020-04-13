module mod_grid

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t

  implicit none

  private
  public :: grid_t

  type, abstract :: grid_t

    ! Node indicies
    integer(ik) :: ilo_bc_node = 0 !< low i boundary condition node index
    integer(ik) :: jlo_bc_node = 0 !< low j boundary condition node index
    integer(ik) :: ihi_bc_node = 0 !< high i boundary condition node index
    integer(ik) :: jhi_bc_node = 0 !< high j boundary condition node index

    integer(ik) :: ilo_node = 0 !< low i node index (not including boundary)
    integer(ik) :: jlo_node = 0 !< low j node index (not including boundary)
    integer(ik) :: ihi_node = 0 !< high i node index (not including boundary)
    integer(ik) :: jhi_node = 0 !< high j node index (not including boundary)
    integer(ik) :: ni_node = 0 !< Number of i nodes (not including boundary nodes)
    integer(ik) :: nj_node = 0 !< Number of j nodes (not including boundary nodes)

    ! Cell indices (this may seem a bit redundant, but the main idea is to improved code readibility)
    integer(ik) :: ilo_bc_cell = 0 !< low i boundary condition cell index
    integer(ik) :: jlo_bc_cell = 0 !< low j boundary condition cell index
    integer(ik) :: ihi_bc_cell = 0 !< high i boundary condition cell index
    integer(ik) :: jhi_bc_cell = 0 !< high j boundary condition cell index

    integer(ik) :: ilo_cell = 0 !< low i cell index (not including boundary)
    integer(ik) :: jlo_cell = 0 !< low j cell index (not including boundary)
    integer(ik) :: ihi_cell = 0 !< high i cell index (not including boundary)
    integer(ik) :: jhi_cell = 0 !< high j cell index (not including boundary)
    integer(ik) :: ni_cell = 0 !< Number of i cells (not including boundary)
    integer(ik) :: nj_cell = 0 !< Number of j cells (not including boundary)

    real(rk) :: xmin = 0.0_rk   !< Min x location
    real(rk) :: xmax = 0.0_rk   !< Max x location
    real(rk) :: min_dx = 0.0_rk !< Minimum spacing in x
    real(rk) :: max_dx = 0.0_rk !< Maximum spacing in x
    real(rk) :: ymin = 0.0_rk   !< Min y location
    real(rk) :: ymax = 0.0_rk   !< Max y location
    real(rk) :: min_dy = 0.0_rk !< Minimum spacing in y
    real(rk) :: max_dy = 0.0_rk !< Maximum spacing in y

    real(rk) :: x_length !< Length of the domain in x
    real(rk) :: y_length !< Length of the domain in y

    ! Rather than keep an array of element types, make arrays that hold the
    ! element information. Use the element type to populate this
    ! once at the beginning of the simulation or on demand

    real(rk), dimension(:, :), allocatable :: cell_volume !< (i,j); volume of each cell
    real(rk), dimension(:, :, :), allocatable :: cell_size !< ((dx,dy), i, j); dx, dy spacing of each cell

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

    real(rk), dimension(:, :), allocatable :: node_x !< (i, j); x location of each node
    real(rk), dimension(:, :), allocatable :: node_y !< (i, j); y location of each node

    real(rk), dimension(:, :, :), allocatable :: cell_centroid_xy
    !< ((x,y), i, j); (x,y) location of the cell centroid

    real(rk), dimension(:, :, :), allocatable :: cell_edge_lengths
    !< ((edge_1:edge_n), i, j); length of each edge

    real(rk), dimension(:, :, :, :, :), allocatable :: cell_node_xy
    !< ((x,y), (point_1:point_n), (node=1,midpoint=2), i, j); The node/midpoint dimension just selects which set of points,
    !< e.g. 1 - all corners, 2 - all midpoints

    real(rk), dimension(:, :, :, :), allocatable :: cell_edge_norm_vectors
    !< ((x,y), edge, i, j); normal direction vector of each face

    real(rk), dimension(:, :, :, :), allocatable :: corner_edge_vectors
    !< ((x,y), (origin, vector_1:4), i, j); edge vector set for each corner
    real(rk), dimension(:, :), allocatable :: corner_edge_vectors_scale

    real(rk), dimension(:, :, :, :), allocatable :: downup_midpoint_edge_vectors
    !< ((x,y), (origin, vector_1:2), i, j); edge vector set for each downup midpoint
    real(rk), dimension(:, :), allocatable :: downup_midpoint_edge_vectors_scale

    real(rk), dimension(:, :, :, :), allocatable :: leftright_midpoint_edge_vectors
    !< ((x,y), (origin, vector_1:2), i, j); edge vector set for each leftright midpoint
    real(rk), dimension(:, :), allocatable :: leftright_midpoint_edge_vectors_scale

    logical :: scale_edge_vectors = .true.

  contains
    procedure(initialize), deferred :: initialize
    procedure(get_2d_data), deferred, public :: get_x
    procedure(get_2d_data), deferred, public :: get_y
    procedure(get_2d_data), deferred, public :: get_cell_volumes
    procedure(get_cell_centroid_xy), deferred, public :: get_cell_centroid_xy
    procedure(get_cell_edge_lengths), deferred, public :: get_cell_edge_lengths
    procedure(get_4d_data), deferred, public :: get_cell_edge_norm_vectors
    procedure(get_midpoint_vectors), deferred, public :: get_midpoint_vectors
    procedure(get_corner_vectors), deferred, public :: get_corner_vectors
    procedure(get_midpoint_vectors_scaled_and_shifted), deferred, public :: get_midpoint_vectors_scaled_and_shifted
    procedure(get_corner_vectors_scaled_and_shifted), deferred, public :: get_corner_vectors_scaled_and_shifted
    procedure(copy_grid), public, deferred :: copy
    generic :: assignment(=) => copy
  end type grid_t

  abstract interface
    subroutine initialize(self, input)
      import :: grid_t
      import :: input_t
      class(grid_t), intent(inout) :: self
      class(input_t), intent(in) :: input
    end subroutine

    subroutine copy_grid(out_grid, in_grid)
      import :: grid_t
      class(grid_t), intent(in) :: in_grid
      class(grid_t), intent(inout) :: out_grid
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

    pure function get_midpoint_vectors(self, cell_ij, edge) result(vectors)
      !< Public interface to get_midpoint_vectors
      import :: grid_t
      import :: ik, rk
      class(grid_t), intent(in) :: self
      integer(ik), dimension(2), intent(in) :: cell_ij
      character(len=*), intent(in) :: edge ! 'bottom', or 'top'
      real(rk), dimension(2, 0:2) :: vectors !< ((x,y), (origin, vector_1, vector_2))
    end function

    pure function get_corner_vectors(self, cell_ij, corner) result(vectors)
      !< Public interface to get_corner_vectors
      import :: grid_t
      import :: ik, rk
      class(grid_t), intent(in) :: self
      integer(ik), dimension(2), intent(in) :: cell_ij
      character(len=*), intent(in) :: corner ! 'lowerleft', 'lowerright', 'upperright', 'upperleft'
      real(rk), dimension(2, 0:4) :: vectors !< ((x,y), (origin, vector_1, vector_2, vector_3, vector_4)
    end function

    pure subroutine get_corner_vectors_scaled_and_shifted(self)
      !< Public interface to get_corner_vectors_scaled_and_shifted
      import :: grid_t
      class(grid_t), intent(inout) :: self
    end subroutine

    pure subroutine get_midpoint_vectors_scaled_and_shifted(self, edge)
      !< Public interface to get_midpoint_vectors_scaled_and_shifted
      import :: grid_t
      class(grid_t), intent(inout) :: self
      character(len=*), intent(in) :: edge ! 'bottom', 'top', 'left', 'right'
    end subroutine

  end interface

end module mod_grid
