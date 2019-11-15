module mod_grid

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_input, only: input_t

  implicit none

  private
  public :: grid_t

  type, abstract :: grid_t
  contains
    procedure(initialize), deferred :: initialize
    procedure(get_2d_data), deferred, public :: get_x
    procedure(get_2d_data), deferred, public :: get_y
    procedure(get_2d_data), deferred, public :: get_cell_volumes
    procedure(get_cell_centroids), deferred, public :: get_cell_centroids
    procedure(get_cell_edge_lengths), deferred, public :: get_cell_edge_lengths
    procedure(get_4d_data), deferred, public :: get_cell_edge_norm_vectors
    procedure(get_cell_node_xy), deferred, public :: get_cell_node_xy
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

    pure function get_cell_centroids(self, i, j, xy) result(cell_centroids)
      !< Public interface to get get_cell_centroids
      import :: grid_t
      import :: ik, rk
      class(grid_t), intent(in) :: self
      integer(ik), intent(in) :: i, j, xy
      real(rk) :: cell_centroids
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
