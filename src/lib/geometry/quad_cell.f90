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
    real(rk) :: min_dx = 0.0_rk  !< minimum edge length in x
    real(rk) :: min_dy = 0.0_rk  !< minimum edge length in y
    real(rk), dimension(2) :: centroid = [0.0_rk, 0.0_rk]  !< centroid x,y coords
    real(rk) :: volume = 0.0_rk !< volume, aka area in 2d

    real(rk), dimension(2, 4) :: edge_norm_vectors = 0.0_rk !< ((x,y), edge); Normal direction vector of each edge
    real(rk), dimension(2, 4) :: edge_midpoints = 0.0_rk !< ((x,y), edge); Midpoint of each edge
    real(rk), dimension(4) :: edge_lengths = 0.0_rk !< (edge); Length of each edge

  contains
    procedure, public :: initialize
    procedure, public :: get_cell_point_coords
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

    ! round-off checks
    associate(y => self%y, x => self%x)
      if(abs(y(2) - y(1)) < epsilon(1.0_rk)) then
        y(2) = y(1)
      end if
      if(abs(y(4) - y(3)) < epsilon(1.0_rk)) then
        y(4) = y(3)
      end if
      if(abs(x(3) - x(2)) < epsilon(1.0_rk)) then
        x(3) = x(2)
      end if
      if(abs(x(4) - x(1)) < epsilon(1.0_rk)) then
        x(4) = x(1)
      end if
    end associate

    call self%calculate_volume()
    call self%calculate_centroid()
    call self%calculate_edge_stats()
    call self%calculate_edge_norm_vectors()

  end subroutine

  subroutine calculate_edge_stats(self)
    !< Find the edge lengths and midpoints of those edges
    class(quad_cell_t), intent(inout) :: self

    associate(x => self%x, y => self%y)
      self%edge_lengths(1) = sqrt((x(2) - x(1))**2 + (y(2) - y(1))**2)
      self%edge_lengths(2) = sqrt((x(3) - x(2))**2 + (y(3) - y(2))**2)
      self%edge_lengths(3) = sqrt((x(4) - x(3))**2 + (y(4) - y(3))**2)
      self%edge_lengths(4) = sqrt((x(1) - x(4))**2 + (y(1) - y(4))**2)

      if(maxval(self%edge_lengths) - minval(self%edge_lengths) < 2.0_rk * epsilon(1.0_rk)) then
        self%edge_lengths = maxval(self%edge_lengths)
      end if

      self%edge_midpoints(:, 1) = [(x(2) + x(1)) / 2.0_rk,(y(2) + y(1)) / 2.0_rk]
      self%edge_midpoints(:, 2) = [(x(3) + x(2)) / 2.0_rk,(y(3) + y(2)) / 2.0_rk]
      self%edge_midpoints(:, 3) = [(x(4) + x(3)) / 2.0_rk,(y(4) + y(3)) / 2.0_rk]
      self%edge_midpoints(:, 4) = [(x(1) + x(4)) / 2.0_rk,(y(1) + y(4)) / 2.0_rk]

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
    ! integer(ik) :: i

    ! Make the x and y arrays wrap around so that the points are
    ! [1, 2, 3, 4, 1] so as to make the math easy
    x(1:4) = self%x; x(5) = self%x(1)
    y(1:4) = self%y; y(5) = self%y(1)

    ! self%centroid = 0.0_rk

    associate(x => self%x, y => self%y)
      self%centroid = [(x(2) + x(1)) / 2.0_rk,(y(3) + y(2)) / 2.0_rk]
    end associate

    ! associate(v=>self%volume, cx=>self%centroid(1), cy=>self%centroid(2))
    !   do i = 1, 4
    !     !FIXME: truncation/roundoff here!
    !     cx = cx + (x(i) + x(i + 1)) * (x(i) * y(i + 1) - x(i + 1) * y(i))
    !     cy = cy + (y(i) + y(i + 1)) * (x(i) * y(i + 1) - x(i + 1) * y(i))
    !   end do
    !   cx = (1.0_rk / (6.0_rk * v)) * cx
    !   cy = (1.0_rk / (6.0_rk * v)) * cy
    ! end associate

  end subroutine

  subroutine calculate_volume(self)
    class(quad_cell_t), intent(inout) :: self
    real(rk) :: dx1, dy1, dx2, dy2

    associate(v => self%volume, x => self%x, y => self%y)
      dx1 = x(2) - x(1)
      dy1 = y(4) - y(1)
      dx2 = x(4) - x(1)
      dy2 = y(2) - y(1)
      if(abs(dx1) < epsilon(1.0_rk)) dx1 = 0.0_rk
      if(abs(dy1) < epsilon(1.0_rk)) dy1 = 0.0_rk
      if(abs(dx2) < epsilon(1.0_rk)) dx2 = 0.0_rk
      if(abs(dy2) < epsilon(1.0_rk)) dy2 = 0.0_rk

      v = (x(2) - x(1)) * (y(4) - y(1)) - (x(4) - x(1)) * (y(2) - y(1))
    end associate

    if(self%volume <= 0.0_rk) then
      write(*, '(a, es10.3)') "Error! Negative Volume: ", self%volume
      write(*, '(2(a, g0.4, ",", g0.4), a)') &
        'N4: (', self%x(4), self%y(4), ')         N3: (', self%x(3), self%y(3), ')'
      write(*, *) '             N4-----o-----N3'
      write(*, *) '             |            |'
      write(*, *) '             o      C     o'
      write(*, *) '             |            |'
      write(*, *) '             N1-----o-----N2'
      write(*, '(2(a, g0.4, ",", g0.4), a)') &
        'N1: (', self%x(1), self%y(1), ')         N2: (', self%x(2), self%y(2), ')'
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
    real(rk), dimension(2) :: vec_x, vec_y

    do i = 1, 4

      associate(x_tail => self%x(tail_idx(i)), &
                x_mid => self%edge_midpoints(1, i), &
                y_mid => self%edge_midpoints(2, i), &
                x_head => self%x(head_idx(i)), &
                y_tail => self%y(tail_idx(i)), &
                y_head => self%y(head_idx(i)), &
                n => self%edge_norm_vectors)

        dx = x_head - x_mid
        dy = y_head - y_mid

        vec_x = [x_mid, x_mid + dy]
        vec_y = [y_mid, y_mid - dx]
        vec = vector_t(x=vec_x, y=vec_y)
        norm_vec = .unitnorm.vec

        n(:, i) = [norm_vec%x, norm_vec%y]

      end associate

    end do

  end subroutine

  pure subroutine get_cell_point_coords(self, x, y)
    class(quad_cell_t), intent(in) :: self
    real(rk), dimension(8), intent(out) :: x !< x coordinates (c1,m1,c2,m2,c3,m3,c4,m4)
    real(rk), dimension(8), intent(out) :: y !< y coordinates (c1,m1,c2,m2,c3,m3,c4,m4)

    ! Corners
    x(1) = self%x(1); y(1) = self%y(1)
    x(3) = self%x(2); y(3) = self%y(2)
    x(5) = self%x(3); y(5) = self%y(3)
    x(7) = self%x(4); y(7) = self%y(4)

    ! Midpoints
    x(2) = self%edge_midpoints(1, 1); y(2) = self%edge_midpoints(2, 1)
    x(4) = self%edge_midpoints(1, 2); y(4) = self%edge_midpoints(2, 2)
    x(6) = self%edge_midpoints(1, 3); y(6) = self%edge_midpoints(2, 3)
    x(8) = self%edge_midpoints(1, 4); y(8) = self%edge_midpoints(2, 4)

  end subroutine

end module mod_quad_cell
