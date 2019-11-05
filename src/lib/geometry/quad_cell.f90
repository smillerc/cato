module mod_quad_cell
  !< Summary: Define the quadrilateral cell/control volume object

  use iso_fortran_env, only : int32, real64

  implicit none
  

  !  Numbering convention for a 2D quadrilateral cell
  !  
  !                     F3
  !                     M3
  !              N4-----o-----N3
  !              |            |  
  !      F4   M4 o      C     o M2   F2
  !              |            |
  !              N1-----o-----N2
  !                     M1
  !                     F1 
  ! N: node or vertex
  ! F: face or edge
  ! M: midpoint of the edge (o)
  ! C: cell or control volume (in finite-volume lingo)

  type quad_cell_t
    real(real64), dimension(4) :: x = 0.0_real64  !< vertex x coords
    real(real64), dimension(4) :: y = 0.0_real64  !< vertex y coords
    real(real64) :: volume = 0.0_real64 !< volume, aka area in 2d

    real(real64), dimension(4,2,2) :: edge_norm_vectors = 0.0_real64 !< normal vectors at each face (face_id, x0:x1, y0:y1)
    real(real64), dimension(4,2) :: edge_midpoints = 0.0_real64 !< midpoint of each edge (face_id, x, y)
    real(real64), dimension(4) :: edge_lengths = 0.0_real64
  contains
    procedure, public :: initialize
    procedure, private :: calculate_volume
    procedure, private :: calculate_edge_stats
    procedure, private :: calculate_edge_norm_vectors
  end type quad_cell_t

contains

  subroutine initialize(self, x_coords, y_coords)
    class(quad_cell_t), intent(inout) :: self
    real(real64), intent(in), dimension(4) :: x_coords
    real(real64), intent(in), dimension(4) :: y_coords

    self%x = x_coords
    self%y = y_coords

    call self%calculate_volume()
    call self%calculate_edge_stats()
    call self%calculate_edge_norm_vectors()

  end subroutine

  subroutine calculate_edge_stats(self)
    class(quad_cell_t), intent(inout) :: self

    associate(x => self%x, y => self%y)
      self%edge_lengths(1) = sqrt((x(2) - x(1))**2 + (y(2) - y(1))**2)
      self%edge_lengths(2) = sqrt((x(3) - x(2))**2 + (y(3) - y(2))**2)
      self%edge_lengths(3) = sqrt((x(4) - x(3))**2 + (y(4) - y(3))**2)
      self%edge_lengths(4) = sqrt((x(1) - x(4))**2 + (y(1) - y(4))**2)

      self%edge_midpoints(1,:) = [0.5*(x(2) + x(1)), 0.5*(y(2) + y(1))]
      self%edge_midpoints(2,:) = [0.5*(x(3) + x(2)), 0.5*(y(3) + y(2))]
      self%edge_midpoints(3,:) = [0.5*(x(4) + x(3)), 0.5*(y(4) + y(3))]
      self%edge_midpoints(4,:) = [0.5*(x(1) + x(4)), 0.5*(y(1) + y(4))]
    end associate

  end subroutine

  subroutine calculate_volume(self)
    class(quad_cell_t), intent(inout) :: self

    associate(v => self%volume, x => self%x, y => self%y)
      v = (x(2) - x(1)) * (y(4) - y(1)) - (x(4) - x(1)) * (y(2) - y(1))
    end associate

    if (self%volume <= 0.0_real64) then
      error stop "Error: Negative volume (found at) "! // __LINE__ // " in " // __FILE__
    end if

  end subroutine

  subroutine calculate_edge_norm_vectors(self)
    !< Find the vector normal to the edge originating at the midpoint 

    class(quad_cell_t), intent(inout) :: self
    integer, dimension(4) :: head_idx = [2,3,4,1]
    integer, dimension(4) :: tail_idx = [1,2,3,4]
    integer :: i
    real(real64) :: dx, dy


    do i = 1,4

      associate(x_tail => self%x(tail_idx(i)), &
                x_mid =>self%edge_midpoints(i,1), &
                x_head => self%x(head_idx(i)), &
                y_tail => self%y(tail_idx(i)), &
                y_mid =>self%edge_midpoints(i,2), &
                y_head => self%y(head_idx(i)), &
                n => self%edge_norm_vectors)

        dx = x_head - x_mid
        dy = y_head - y_mid

        n(i,1,:) = [x_mid, y_mid]
        n(i,2,:) = [x_mid + dy, y_mid - dx]

        print*, 'i', i
        ! print*, n(i,1,:)
        print*, n(i,2,:)
        print*


      end associate

    end do
    
  end subroutine

end module mod_quad_cell