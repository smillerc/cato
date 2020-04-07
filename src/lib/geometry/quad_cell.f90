module mod_quad_cell
  !< Summary: Define the quadrilateral cell/control volume object

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_vector, only: vector_t, operator(.unitnorm.), operator(.dot.), operator(.cross.)
  implicit none

  private
  public :: quad_cell_t
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
    real(rk), dimension(4) :: min_dx = 0.0_rk  !< minimum edge length in x
    real(rk), dimension(4) :: min_dy = 0.0_rk  !< minimum edge length in x
    real(rk), dimension(2) :: centroid = [0.0_rk, 0.0_rk]  !< centroid x,y coords
    real(rk) :: volume = 0.0_rk !< volume, aka area in 2d

    real(rk), dimension(2, 4) :: edge_norm_vectors = 0.0_rk !< ((x,y), edge); Normal direction vector of each edge
    real(rk), dimension(2, 4) :: edge_midpoints = 0.0_rk !< ((x,y), edge); Midpoint of each edge
    real(rk), dimension(4) :: edge_lengths = 0.0_rk !< (edge); Length of each edge

  contains
    procedure, public :: initialize
    procedure, public :: get_cell_node_xy_set
    procedure, private :: calculate_volume
    procedure, private :: calculate_centroid
    procedure, private :: calculate_edge_stats
    procedure, private :: calculate_edge_norm_vectors
  end type quad_cell_t

contains

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

      self%edge_midpoints(:, 1) = [0.5 * (x(2) + x(1)), 0.5 * (y(2) + y(1))]
      self%edge_midpoints(:, 2) = [0.5 * (x(3) + x(2)), 0.5 * (y(3) + y(2))]
      self%edge_midpoints(:, 3) = [0.5 * (x(4) + x(3)), 0.5 * (y(4) + y(3))]
      self%edge_midpoints(:, 4) = [0.5 * (x(1) + x(4)), 0.5 * (y(1) + y(4))]

      self%min_dx = min(abs(x(2) - x(1)), &
                        abs(x(3) - x(4)))

      self%min_dy = min(abs(y(4) - y(1)), &
                        abs(y(3) - y(2)))

    end associate

  end subroutine

  subroutine calculate_centroid(self)
    !< Find the centroid x and y location
    !< This uses the formula for a polygon found here (https://en.wikipedia.org/wiki/Centroid#Of_a_polygon)
    class(quad_cell_t), intent(inout) :: self

    real(rk), dimension(5) :: x, y
    integer(ik) :: i

    ! Make the x and y arrays wrap around so that the points are
    ! [1, 2, 3, 4, 1] so as to make the math easy
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

    type(vector_t) :: vec, norm_vec
    integer, dimension(4) :: head_idx = [2, 3, 4, 1]
    integer, dimension(4) :: tail_idx = [1, 2, 3, 4]
    integer :: i
    real(rk) :: dx, dy

    do i = 1, 4

      associate(x_tail=>self%x(tail_idx(i)), &
                x_mid=>self%edge_midpoints(1, i), &
                y_mid=>self%edge_midpoints(2, i), &
                x_head=>self%x(head_idx(i)), &
                y_tail=>self%y(tail_idx(i)), &
                y_head=>self%y(head_idx(i)), &
                n=>self%edge_norm_vectors)

        dx = x_head - x_mid
        dy = y_head - y_mid

        vec = vector_t(x=[x_mid, x_mid + dy], y=[y_mid, y_mid - dx])
        norm_vec = .unitnorm.vec

        n(:, i) = [norm_vec%x, norm_vec%y]

      end associate

    end do

  end subroutine

  pure function get_cell_node_xy_set(self) result(set)
    class(quad_cell_t), intent(in) :: self
    real(rk), dimension(2, 4, 2) :: set  !< ((x,y), (point_1:point_4), (corner=1, midpoint=2))
    set(1, :, 1) = self%x
    set(2, :, 1) = self%y
    set(:, :, 2) = self%edge_midpoints
  end function

end module mod_quad_cell
