module mod_mach_cone_abstract
  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64
  use mod_floating_point_utils, only: near_zero, equal
  use mod_vector, only: vector_t, operator(.cross.)
  use mod_geometry, only: super_circle, get_arc_segments_new

  implicit none

  private
  public :: mach_cone_abstract_t

  type, abstract :: mach_cone_abstract_t
    type(vector_t) :: p_prime_vector !< vector from P0 to P' (P0 is shifted to (0,0))

    real(rk) :: pressure_p_prime = 0.0_rk
    !< pressure from the P' cell (due to FVLEG assumptions, this is the reconstructed value at P0 of the cell that contains P')

    real(rk) :: density_p_prime = 0.0_rk
    !< density from the P' cell (due to FVLEG assumptions, this is the reconstructed value at P0 of the cell that contains P')

    real(rk), dimension(2) :: p_prime_xy = 0.0_rk !< (x,y); Location of P'
    real(rk), dimension(2) :: p_xy = 0.0_rk       !< (x,y); Location of P
    integer(ik) :: n_neighbor_cells = 0           !< Number of neighbor cells
    real(rk) :: radius = -1.0_rk                  !< Radius of the cone
    real(rk) :: tau = 1.0e-10_rk                  !< Time evolution increment
    real(rk) :: reference_density = 0.0_rk        !< Reference density (e.g. neighbor averaged)
    real(rk) :: reference_u = 0.0_rk              !< Reference x velocity (e.g. neighbor averaged)
    real(rk) :: reference_v = 0.0_rk              !< Reference y velocity (e.g. neighbor averaged)
    real(rk) :: reference_sound_speed = 0.0_rk    !< Reference sound speed (e.g. neighbor averaged)
    logical :: cone_is_transonic = .false.        !< Flag to enable special treatment for transonic cones
    logical :: cone_is_centered = .false.
    !< Flag signifying if P and P' are collocated -> allows for simplified angles/trig
    !< so a bunch of functions can be skipped for speed (no need to find line/circle intersections)

    character(len=32) :: cone_location = ''       !< What type of cone am I?
  contains
    private
    ! Deferred procedures
    procedure(ref_state_interface), deferred :: get_reference_state
    procedure(basic_interface), deferred :: precompute_trig_angles
    procedure(io_interface), deferred, public :: write
    procedure(arc_angles_interface), deferred :: find_arc_angles

    ! Parent class procedures (used by concrete cone classes)
    procedure, non_overridable :: determine_if_p_prime_is_in_cell
    procedure, non_overridable :: get_arc_segments
    procedure, non_overridable :: get_intersection_angles
    procedure, non_overridable :: find_line_circle_intersections
    procedure, non_overridable :: intersection_angle_from_x_axis

    ! Generic routines
    generic, public :: write(formatted) => write
  end type mach_cone_abstract_t

  abstract interface
    subroutine basic_interface(self)
      import :: mach_cone_abstract_t
      class(mach_cone_abstract_t), intent(inout) :: self
    end subroutine basic_interface

    subroutine io_interface(self, unit, iotype, v_list, iostat, iomsg)
      import :: mach_cone_abstract_t
      class(mach_cone_abstract_t), intent(in) :: self     !< cone class
      integer, intent(in) :: unit           !< input/output unit
      character(*), intent(in) :: iotype    !< LISTDIRECTED or DTxxx
      integer, intent(in) :: v_list(:)      !< parameters from fmt spec.
      integer, intent(out) :: iostat        !< non zero on error, etc.
      character(*), intent(inout) :: iomsg  !< define if iostat non zero.
    end subroutine io_interface

    subroutine ref_stat_interface(self, reconstructed_state)
      import :: mach_cone_abstract_t, rk
      class(mach_cone_abstract_t), intent(inout) :: self
      real(rk), dimension(:, :), intent(in) :: reconstructed_state !< [rho, u, v, a]
    end subroutine ref_stat_interface

    subroutine arc_angles_interface(self, cell_indices, edge_vectors, theta_ib, theta_ie, n_arcs_per_cell)
      import :: mach_cone_abstract_t, ik, rk
      class(mach_cone_abstract_t), intent(inout) :: self
      real(rk), dimension(:, :), intent(in) :: edge_vectors
      integer(ik), dimension(:, :), intent(in) :: cell_indices
      !< ((i,j), cell 1:N_CELLS); set of indices for the neighboring cells -> needed to find P' i,j index
      real(rk), dimension(:, :), intent(out) :: theta_ib
      !< ((arc1, arc2), (cell 1:N_CELLS)); starting angle [rad] for the arc contained in each cell
      real(rk), dimension(:, :), intent(out) :: theta_ie
      !< ((arc1, arc2), (cell 1:N_CELLS)); ending angle [rad] for the arc contained in each cell
      integer(ik), dimension(:), intent(out) :: n_arcs_per_cell !< # of arcs in a given cell
    end subroutine arc_angles_interface
  end interface
contains

  pure logical function determine_if_p_prime_is_in_cell(self, vector_1, vector_2) result(in_cell)
    !< Implementation of whether the P' point is inside the current cell/control volume. This
    !< uses the cross product of 2 vectors in 2d, which gives a scalar

    class(mach_cone_abstract_t), intent(in) :: self
    type(vector_t), intent(in) :: vector_1, vector_2

    if((self%p_prime_vector.cross.edge_vector_1) <= 0.0_rk .and. &
       (edge_vector_2.cross.self%p_prime_vector) <= 0.0_rk) then
      in_cell = .true.
    else
      in_cell = .false.
    end if
  end function determine_if_p_prime_is_in_cell

  pure subroutine get_arc_segments(self, vector_1, vector_2, origin_in_cell, arc_segments, n_arcs)
    !< Given 2 lines and a circle, find their intersections and starting/ending angles for each arc

    class(mach_cone_abstract_t), intent(in) :: self
    type(vector_t), intent(in) :: vector_1, vector_2

    logical, intent(in) :: origin_in_cell  !< is the circle origin in the cell defined by the two lines?

    ! Output
    real(rk), dimension(2, 2), intent(out):: arc_segments !< ((start,end), (arc1, arc2))
    integer(ik), intent(out) :: n_arcs !< Number of arcs with a valid start and end angle

    ! Dummy
    integer(ik) :: n_intersections_per_line !< # intersections for a single line
    integer(ik), dimension(2) :: n_intersections !< # intersections for all lines
    real(rk), dimension(2, 2) :: line !< ((x,y), (tail,head)); Single line point locations
    real(rk), dimension(2) :: intersection_angles_per_line !< intersection angles
    real(rk), dimension(2, 2) :: intersection_angles !< intersection angles for all lines
    integer(ik) :: i
    logical, dimension(2) :: valid_intersections
    logical, dimension(2, 2) :: total_valid_intersections !< ((intersection 1, intersection 2), (line 1, line 2)); valid intersections for each line

    valid_intersections = .false.
    total_valid_intersections = .false.
    n_intersections_per_line = 0
    n_intersections = 0
    intersection_angles_per_line = 0.0_rk
    intersection_angles = 0.0_rk

    line(:, 1) = origin ! both lines share the same origin, so only do this once

    ! Line 1
    line(:, 2) = vec_1_head
    ! For a given line & circle intersection, find the angle with respect to the x-axis for each intersection point
    call get_intersection_angles(line_xy=line, self%p_prime_xy=self%p_prime_xy, self%radius=self%radius, &
                                 arc_angles=intersection_angles_per_line, valid_intersections=valid_intersections)
    intersection_angles(1, :) = intersection_angles_per_line
    n_intersections(1) = count(valid_intersections)
    total_valid_intersections(:, 1) = valid_intersections

    ! Line 2
    line(:, 2) = vec_2_head
    ! For a given line & circle intersection, find the angle with respect to the x-axis for each intersection point
    call get_intersection_angles(line_xy=line, self%p_prime_xy=self%p_prime_xy, self%radius=self%radius, &
                                 arc_angles=intersection_angles_per_line, valid_intersections=valid_intersections)
    intersection_angles(2, :) = intersection_angles_per_line
    n_intersections(2) = count(valid_intersections)
    total_valid_intersections(:, 2) = valid_intersections

    ! Find arc starting and ending angles
    call get_theta_start_end(thetas=intersection_angles, origin_in_cell=origin_in_cell, &
                             valid_intersections=total_valid_intersections, &
                             n_intersections=n_intersections, &
                             theta_start_end=arc_segments, n_arcs=n_arcs)
  end subroutine get_arc_segments

  pure subroutine get_intersection_angles(self, line_xy, arc_angles, valid_intersections)
    !< Given a arbitrary line from (x1,y1) to (x2,y2) and a circle at (x,y) with a given radius, find
    !< the angle that a the vector from the circle's center to the intersection point(s) has with respect
    !< to the x-axis

    class(mach_cone_abstract_t), intent(in) :: self
    real(rk), dimension(2, 2), intent(in) :: line_xy !< ((x,y), (point_1, point_2))

    ! Output
    logical, dimension(2), intent(out) :: valid_intersections !< (point_1, point_2); .true. or .false.
    real(rk), dimension(2), intent(out):: arc_angles

    ! Dummy
    integer(ik) :: i
    real(rk), dimension(2, 2) :: intersection_xy !< ((x,y), (point_1, point_2))

    intersection_xy = 0.0_rk
    valid_intersections = .false.

    ! Find the intersection (x,y) locations (if any)
    call find_line_circle_intersections(line_xy, intersection_xy, valid_intersections)

    arc_angles = 0.0_rk
    do i = 1, 2
      if(valid_intersections(i)) then
        arc_angles(i) = intersection_angle_from_x_axis(self%p_prime_xy, intersection_xy(:, i))
      end if
    end do
  end subroutine get_intersection_angles

  pure subroutine find_line_circle_intersections(self, line_xy, intersection_xy, valid_intersection)
    !< Find the intersections between an arbitrary line and circle. There can be 0, 1, or 2 intersections.

    class(mach_cone_abstract_t), intent(in) :: self
    real(rk), dimension(2, 2), intent(in) :: line_xy !< ((x,y) (point_1, point_2)); Line location
    real(rk), dimension(2, 2), intent(out) :: intersection_xy !< ((x,y), (point_1, point_2)); Intersection point
    logical, dimension(2), intent(out) :: valid_intersection !< (point_1, point_2); Is there an intersection or not?

    real(rk) :: discriminiant  !< term under the square root in the quadratic formula
    real(rk) :: sqrt_discriminiant  !< term under the square root in the quadratic formula
    integer(ik) :: i
    real(rk), dimension(2) :: t !< scale factor (should be between 0 and 1) of where the intersection point is along the line
    real(rk) :: a, b, c  !< quadratic formula variables

    t = 0.0_rk
    discriminiant = 0.0_rk
    a = 0.0_rk
    b = 0.0_rk
    c = 0.0_rk

    associate(x0=>line_xy(1, 1), x1=>line_xy(1, 2), y0=>line_xy(2, 1), y1=>line_xy(2, 2), &
              r=>self%radius, h=>self%p_prime_xy(1), k=>self%p_prime_xy(2))
      a = (x1 - x0)**2 + (y1 - y0)**2
      b = 2 * (x1 - x0) * (x0 - h) + 2 * (y1 - y0) * (y0 - k)
      c = (x0 - h)**2 + (y0 - k)**2 - r**2
    end associate

    discriminiant = b**2 - 4 * a * c

    if(discriminiant > 0.0_rk) then
      sqrt_discriminiant = sqrt(discriminiant)

      ! This used the alternative quadratic formula better suited for floating point operations
      if(near_zero(-b + sqrt_discriminiant)) then
        t(1) = 0.0_rk ! t_1 -> 0 when intersection is at the vector start point
      else
        t(1) = (2 * c) / (-b + sqrt_discriminiant)
      end if

      if(near_zero(-b - sqrt_discriminiant)) then
        t(2) = 1.0_rk  ! t_2 -> 1 when intersection is at the vector end point
      else
        t(2) = (2 * c) / (-b - sqrt_discriminiant)
      end if
    end if

    ! The scale factor t must be between 0 and 1, otherwise there is no intersection
    valid_intersection = .false.
    if(discriminiant > 0.0_rk) then
      do i = 1, 2 ! TODO: Can we use a where block?
        if(t(i) < 0.0_rk .or. t(i) > 1.0_rk) then
          valid_intersection(i) = .false.
        else
          valid_intersection(i) = .true.
        end if
      end do
    end if

    intersection_xy = 0.0_rk
    if(discriminiant > 0.0_rk) then
      associate(x0=>line_xy(1, 1), x1=>line_xy(1, 2), &
                y0=>line_xy(2, 1), y1=>line_xy(2, 2))
        ! x_t
        intersection_xy(1, :) = (x1 - x0) * t + x0

        ! y_t
        intersection_xy(2, :) = (y1 - y0) * t + y0
      end associate
    end if
  end subroutine find_line_circle_intersections

  pure real(rk) function intersection_angle_from_x_axis(self, intersection_xy) result(angle)
    !< Given the the intersection point and origin of the circle it intersected, determine the
    !< angle with respect to the x axis from 0 to 2pi

    class(mach_cone_abstract_t), intent(in) :: self
    real(rk), dimension(2), intent(in) :: intersection_xy !< ((x,y)

    angle = atan2(y=intersection_xy(2) - self%p_prime_xy(2), &
                  x=intersection_xy(1) - self%p_prime_xy(1))
  end function intersection_angle_from_x_axis

end module mod_mach_cone_abstract
