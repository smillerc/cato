module mod_quad_cell
  !< Summary: Define the quadrilateral cell/control volume object

  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_vector_2d, only: vector_2d_t, operator(.unitnorm.)

  implicit none

  private
  public :: quad_cell_t, quad_cell_edge_vector_mapping_t
  !<  Numbering convention for a 2D quadrilateral cell
  !<
  !<                     F3
  !<                     M3
  !<              N4-----o-----N3
  !<              |            |
  !<      F4   M4 o      C     o M2   F2
  !<              |            |
  !<              N1-----o-----N2
  !<                     M1
  !<                     F1
  !< N: node or vertex
  !< F: face or edge
  !< M: midpoint of the edge (o)
  !< C: cell or control volume (in finite-volume lingo)

  type quad_cell_t
    real(rk), dimension(4) :: x = 0.0_rk  !< vertex x coords
    real(rk), dimension(4) :: y = 0.0_rk  !< vertex y coords
    real(rk), dimension(2) :: centroid = [0.0_rk, 0.0_rk]  !< centroid x,y coords
    real(rk) :: volume = 0.0_rk !< volume, aka area in 2d

    real(rk), dimension(4, 2, 2) :: edge_norm_vectors = 0.0_rk !< normal vectors at each face (face_id, x0:x1, y0:y1)
    real(rk), dimension(4, 2) :: edge_midpoints = 0.0_rk !< midpoint of each edge (face_id, x, y)
    real(rk), dimension(4) :: edge_lengths = 0.0_rk
  contains
    procedure, public :: initialize
    procedure, private :: calculate_volume
    procedure, private :: calculate_centroid
    procedure, private :: calculate_edge_stats
    procedure, private :: calculate_edge_norm_vectors
  end type quad_cell_t

  integer(ik), parameter :: N1 = 1
  integer(ik), parameter :: M1 = 2
  integer(ik), parameter :: N2 = 3
  integer(ik), parameter :: M2 = 4
  integer(ik), parameter :: N3 = 5
  integer(ik), parameter :: M3 = 6
  integer(ik), parameter :: N4 = 7
  integer(ik), parameter :: M4 = 8

  type :: quad_cell_edge_vector_mapping_t
    integer(ik), dimension(4, 3, 3) :: edge_to_loc
  end type quad_cell_edge_vector_mapping_t

  interface quad_cell_edge_vector_mapping_t
    module procedure mapping_constructor
  end interface

contains

  pure function mapping_constructor() result(map)
    type(quad_cell_edge_vector_mapping_t) :: map

    map%edge_to_loc(1, 1, :) = [N2, N1, N4] ! k=1,1
    map%edge_to_loc(1, 2, :) = [N2, M1, N1] ! k=1,c
    map%edge_to_loc(1, 3, :) = [N3, N2, N1] ! k=1,2

    map%edge_to_loc(2, 1, :) = [N3, N2, N1] ! k=2,1
    map%edge_to_loc(2, 2, :) = [N3, M2, N2] ! k=2,c
    map%edge_to_loc(2, 3, :) = [N4, N3, N2] ! k=2,2

    map%edge_to_loc(3, 1, :) = [N4, N3, N2] ! k=3,1
    map%edge_to_loc(3, 2, :) = [N4, M3, N3] ! k=3,c
    map%edge_to_loc(3, 3, :) = [N1, N4, N3] ! k=3,2

    map%edge_to_loc(4, 1, :) = [N1, N4, N3] ! k=4,1
    map%edge_to_loc(4, 2, :) = [N1, M4, N4] ! k=4,c
    map%edge_to_loc(4, 3, :) = [N2, N1, N4] ! k=4,2
  end function

  subroutine initialize(self, x_coords, y_coords)
    class(quad_cell_t), intent(inout) :: self
    real(rk), intent(in), dimension(4) :: x_coords
    real(rk), intent(in), dimension(4) :: y_coords

    self%x = x_coords
    self%y = y_coords

    call self%calculate_volume()
    call self%calculate_centroid()
    call self%calculate_edge_stats()
    call self%calculate_edge_norm_vectors()

  end subroutine

  subroutine calculate_edge_stats(self)
    !< Find the edge lengths and midpoints of those edges
    class(quad_cell_t), intent(inout) :: self

    associate(x=>self%x, y=>self%y)
      self%edge_lengths(1) = sqrt((x(2) - x(1))**2 + (y(2) - y(1))**2)
      self%edge_lengths(2) = sqrt((x(3) - x(2))**2 + (y(3) - y(2))**2)
      self%edge_lengths(3) = sqrt((x(4) - x(3))**2 + (y(4) - y(3))**2)
      self%edge_lengths(4) = sqrt((x(1) - x(4))**2 + (y(1) - y(4))**2)

      self%edge_midpoints(1, :) = [0.5 * (x(2) + x(1)), 0.5 * (y(2) + y(1))]
      self%edge_midpoints(2, :) = [0.5 * (x(3) + x(2)), 0.5 * (y(3) + y(2))]
      self%edge_midpoints(3, :) = [0.5 * (x(4) + x(3)), 0.5 * (y(4) + y(3))]
      self%edge_midpoints(4, :) = [0.5 * (x(1) + x(4)), 0.5 * (y(1) + y(4))]
    end associate

  end subroutine

  subroutine calculate_centroid(self)
    !< Find the centroid x and y location
    !< This uses the formula for a polygon found here (https://en.wikipedia.org/wiki/Centroid#Of_a_polygon)
    class(quad_cell_t), intent(inout) :: self

    real(rk), dimension(5) :: x, y
    integer(ik) :: i

    x(1:4) = self%x; x(5) = self%x(1)
    y(1:4) = self%y; y(5) = self%y(1)

    self%centroid = 0.0_rk

    associate(v=>self%volume, cx=>self%centroid(1), cy=>self%centroid(2))
      do i = 1, 4
        cx = cx + (x(i) + x(i + 1)) * (x(i) * y(i + 1) - x(i + 1) * y(i))
        cy = cy + (y(i) + y(i + 1)) * (x(i) * y(i + 1) - x(i + 1) * y(i))
      end do
      cx = (1.0_rk / (6.0_rk * v)) * cx
      cy = (1.0_rk / (6.0_rk * v)) * cy
    end associate

  end subroutine

  subroutine calculate_volume(self)
    class(quad_cell_t), intent(inout) :: self

    associate(v=>self%volume, x=>self%x, y=>self%y)
      v = (x(2) - x(1)) * (y(4) - y(1)) - (x(4) - x(1)) * (y(2) - y(1))
    end associate

    if(self%volume <= 0.0_rk) then
      write(*, '(a, 2(g0.4, 1x))') 'N4: ', self%x(4), self%y(4)
      write(*, '(a, 2(g0.4, 1x))') 'N3: ', self%x(3), self%y(3)
      write(*, *) '             N4-----o-----N3'
      write(*, *) '             |            |'
      write(*, *) '             o      C     o'
      write(*, *) '             |            |'
      write(*, *) '             N1-----o-----N2'
      write(*, '(a, 2(g0.4, 1x))') 'N1: ', self%x(1), self%y(1)
      write(*, '(a, 2(g0.4, 1x))') 'N2: ', self%x(2), self%y(2)
      error stop "Negative volume!"
    end if

  end subroutine

  subroutine calculate_edge_norm_vectors(self)
    !< Find the vector normal to the edge originating at the midpoint

    class(quad_cell_t), intent(inout) :: self

    type(vector_2d_t) :: vec, norm_vec
    integer, dimension(4) :: head_idx = [2, 3, 4, 1]
    integer, dimension(4) :: tail_idx = [1, 2, 3, 4]
    integer :: i
    real(rk) :: dx, dy

    do i = 1, 4

      associate(x_tail=>self%x(tail_idx(i)), &
                x_mid=>self%edge_midpoints(i, 1), &
                x_head=>self%x(head_idx(i)), &
                y_tail=>self%y(tail_idx(i)), &
                y_mid=>self%edge_midpoints(i, 2), &
                y_head=>self%y(head_idx(i)), &
                n=>self%edge_norm_vectors)

        dx = x_head - x_mid
        dy = y_head - y_mid

        call vec%initialize(x_coords=[x_mid, x_mid + dy], y_coords=[y_mid, y_mid - dx])
        norm_vec = .unitnorm.vec

        n(i, 1, :) = [norm_vec%x(1), norm_vec%y(1)]
        n(i, 2, :) = [norm_vec%x(2), norm_vec%y(2)]

      end associate

    end do

  end subroutine

end module mod_quad_cell
