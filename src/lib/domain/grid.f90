module mod_grid
  use iso_fortran_env, only: ik => int32, rk => real64
  use mod_quad_cell, only: quad_cell_t
  use mod_input, only: input_t

  implicit none

  private
  public :: grid_t

  type :: grid_t
    !< Summary: The grid_t type holds all of the geometry info relevant to the grid.
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

    real(rk), dimension(:, :), allocatable :: x
    real(rk), dimension(:, :), allocatable :: y

    ! Rather than keep an array of element types, make arrays that hold the
    ! element information. Use the element type to populate this
    ! once at the beginning of the simulation or on demand

    real(rk), dimension(:, :), allocatable :: cell_volumes !< (i,j)
    real(rk), dimension(:, :, :), allocatable :: cell_centroids !< (i, j, x:y)
    real(rk), dimension(:, :, :), allocatable :: cell_edge_lengths  !< (i, j, face1:face4)
    real(rk), dimension(:, :, :, :), allocatable :: cell_edge_midpoints !< (i, j, face1:face4, x:y)
    real(rk), dimension(:, :, :, :, :), allocatable :: cell_edge_norm_vectors  !< (i, j, face1:face4, x:1:x2, y1:y2)

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
    ! procedure, public :: get_x
    ! procedure, public :: get_y
    ! procedure, public :: get_cell_volumes
    ! procedure, public :: get_cell_centroids
    ! procedure, public :: get_cell_edge_lengths
    ! procedure, public :: get_cell_edge_midpoints
    ! procedure, public :: get_cell_edge_norm_vectors
    final :: finalize
  end type grid_t

contains

  subroutine initialize(self, input)

    class(grid_t), intent(inout) :: self
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

    allocate(self%x(self%ilo:self%ihi, self%jlo:self%jhi), &
             stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate grid_t%x"

    do i = self%ilo, self%ihi
      self%x(i, :) = (i - 1) * self%dx
    end do

    allocate(self%y(self%ilo:self%ihi, self%jlo:self%jhi), &
             stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate grid_t%y"

    do j = self%jlo, self%jhi
      self%y(:, j) = (j - 1) * self%dy
    end do

    allocate(self%cell_volumes(self%ilo:self%ihi - 1, self%jlo:self%jhi - 1), &
             stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate grid_t%cell_volumes"
    self%cell_volumes = 0.0_rk

    allocate(self%cell_centroids(self%ilo:self%ihi - 1, self%jlo:self%jhi - 1, 2), &
             stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate grid_t%cell_centroids"
    self%cell_centroids = 0.0_rk

    allocate(self%cell_edge_lengths(self%ilo:self%ihi - 1, self%jlo:self%jhi - 1, 4), &
             stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate grid_t%cell_edge_lengths"
    self%cell_edge_lengths = 0.0_rk

    allocate(self%cell_edge_midpoints(self%ilo:self%ihi - 1, self%jlo:self%jhi - 1, 4, 2), &
             stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate grid_t%cell_edge_midpoints"
    self%cell_edge_midpoints = 0.0_rk

    allocate(self%cell_edge_norm_vectors(self%ilo:self%ihi - 1, self%jlo:self%jhi - 1, 4, 2, 2), &
             stat=alloc_status)
    if(alloc_status /= 0) error stop "Unable to allocate grid_t%cell_edge_norm_vectors"
    self%cell_edge_norm_vectors = 0.0_rk

    call self%populate_element_specifications()
  end subroutine

  subroutine populate_element_specifications(self)
    !< Summary: Fill the element arrays up with the geometric information
    !< This seemed to be better for memory access patterns elsewhere in the code. Fortran prefers
    !< and structure of arrays rather than an array of structures

    class(grid_t), intent(inout) :: self
    class(quad_cell_t), allocatable :: quad

    integer(ik) :: i, j

    do j = self%jlo, self%jhi - 1
      do i = self%ilo, self%ihi - 1

        allocate(quad_cell_t :: quad)

        associate(x=>self%x, y=>self%y)
          call quad%initialize(x_coords=[x(i, j), x(i + 1, j), x(i + 1, j + 1), x(i, j + 1)], &
                               y_coords=[y(i, j), y(i + 1, j), y(i + 1, j + 1), y(i, j + 1)])
        end associate

        self%cell_volumes(i, j) = quad%volume
        self%cell_centroids(i, j, :) = quad%centroid
        self%cell_edge_lengths(i, j, :) = quad%edge_lengths
        self%cell_edge_midpoints(i, j, :, :) = quad%edge_midpoints
        self%cell_edge_norm_vectors(i, j, :, :, :) = quad%edge_norm_vectors

        deallocate(quad)

      end do
    end do

  end subroutine

  subroutine finalize(self)
    type(grid_t), intent(inout) :: self

    print *, 'Finalizing grid_t'
    if(allocated(self%x)) deallocate(self%x)
    if(allocated(self%y)) deallocate(self%y)
    if(allocated(self%cell_volumes)) deallocate(self%cell_volumes)
    if(allocated(self%cell_centroids)) deallocate(self%cell_centroids)
    if(allocated(self%cell_edge_lengths)) deallocate(self%cell_edge_lengths)
    if(allocated(self%cell_edge_midpoints)) deallocate(self%cell_edge_midpoints)
    if(allocated(self%cell_edge_norm_vectors)) deallocate(self%cell_edge_norm_vectors)
    print *, 'Done'
  end subroutine

  pure function get_ihi(self) result(ihi)
    !< Public interface to get ihi
    class(grid_t), intent(in) :: self
    integer(ik) :: ihi
    ihi = self%ihi
  end function

  pure function get_ilo(self) result(ilo)
    !< Public interface to get ilo
    class(grid_t), intent(in) :: self
    integer(ik) :: ilo
    ilo = self%ilo
  end function

  pure function get_ni(self) result(ni)
    !< Public interface to get ni
    class(grid_t), intent(in) :: self
    integer(ik) :: ni
    ni = self%ni
  end function

  pure function get_nj(self) result(nj)
    !< Public interface to get nj
    class(grid_t), intent(in) :: self
    integer(ik) :: nj
    nj = self%nj
  end function

  pure function get_jlo(self) result(jlo)
    !< Public interface to get jlo
    class(grid_t), intent(in) :: self
    integer(ik) :: jlo
    jlo = self%jlo
  end function

  pure function get_jhi(self) result(jhi)
    !< Public interface to get jhi
    class(grid_t), intent(in) :: self
    integer(ik) :: jhi
    jhi = self%jhi
  end function

  pure function get_xmin(self) result(xmin)
    !< Public interface to get xmin
    class(grid_t), intent(in) :: self
    real(rk) :: xmin
    xmin = self%xmin
  end function

  pure function get_xmax(self) result(xmax)
    !< Public interface to get xmax
    class(grid_t), intent(in) :: self
    real(rk) :: xmax
    xmax = self%xmax
  end function

  pure function get_ymin(self) result(ymin)
    !< Public interface to get ymin
    class(grid_t), intent(in) :: self
    real(rk) :: ymin
    ymin = self%ymin
  end function

  pure function get_ymax(self) result(ymax)
    !< Public interface to get ymax
    class(grid_t), intent(in) :: self
    real(rk) :: ymax
    ymax = self%ymax
  end function

  pure function get_x_length(self) result(x_length)
    !< Public interface to get x_length
    class(grid_t), intent(in) :: self
    real(rk) :: x_length
    x_length = self%x_length
  end function

  pure function get_y_length(self) result(y_length)
    !< Public interface to get y_length
    class(grid_t), intent(in) :: self
    real(rk) :: y_length
    y_length = self%y_length
  end function

  ! pure subroutine get_x(self,x)
  !   !< Public interface to get x
  !   class(grid_t), intent(in) :: self
  !   real(rk) :: x
  !   x = self%x
  ! end subroutine

  ! pure subroutine get_y(self,y)
  !   !< Public interface to get y
  !   class(grid_t), intent(in) :: self
  !   real(rk) :: y
  !   y = self%y
  ! end subroutine

  ! pure subroutine get_cell_volumes(self,cell_volumes)
  !   !< Public interface to get cell_volumes
  !   class(grid_t), intent(in) :: self
  !   real(rk), dimension(:,:) :: cell_volumes
  !   cell_volumes = self%cell_volumes
  ! end subroutine

  ! pure subroutine get_cell_centroids(self,cell_centroids)
  !   !< Public interface to get cell_centroids
  !   class(grid_t), intent(in) :: self
  !   real(rk), dimension(:,:,:) :: cell_centroids
  !   cell_centroids = self%cell_centroids
  ! end subroutine

  ! pure subroutine get_cell_edge_lengths(self,cell_edge_lengths)
  !   !< Public interface to get cell_edge_lengths
  !   class(grid_t), intent(in) :: self
  !   real(rk), dimension(:,:,:) :: cell_edge_lengths
  !   cell_edge_lengths = self%cell_edge_lengths
  ! end subroutine

  ! pure subroutine get_cell_edge_midpoints(self,cell_edge_midpoints)
  !   !< Public interface to get cell_edge_midpoints
  !   class(grid_t), intent(in) :: self
  !   real(rk), dimension(:,:,:,:) :: cell_edge_midpoints
  !   cell_edge_midpoints = self%cell_edge_midpoints
  ! end subroutine

  ! pure subroutine get_cell_edge_norm_vectors(self,cell_edge_norm_vectors)
  !   !< Public interface to get cell_edge_norm_vectors
  !   class(grid_t), intent(in) :: self
  !   real(rk), dimension(:,:,:,:) :: cell_edge_norm_vectors
  !   cell_edge_norm_vectors = self%cell_edge_norm_vectors
  ! end subroutine

end module mod_grid
