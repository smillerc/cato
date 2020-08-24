! MIT License
! Copyright (c) 2020 Sam Miller
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

module mod_field
  !< Summary: Provide the base 2D field class. This originates from https://github.com/modern-fortran/tsunami
  !< Date: 08/18/2020
  !< Author: Milan Curcic, Sam Miller (minor mods)
  !< Notes:
  !< References:
  !      [1] Milan Curcic, "Modern Fortran: Building efficient parallel applications", 2020

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  ! use mod_parallel, only: tile_indices, tile_neighbors_2d

  implicit none

  private
  public :: field_t

  ! Neighbor image indices
  integer(ik), parameter :: left = 1  !< left neighbor
  integer(ik), parameter :: right = 2 !< right neighbor
  integer(ik), parameter :: down = 3  !< neighbor below
  integer(ik), parameter :: up = 4    !< neighbor above

  type :: field_t
    !< A base class that encapsulates 2D data with knowledge of its parallel neighbors and bounds

    real(rk), allocatable :: data(:, :)           !< (i, j); field data
    character(:), allocatable :: name             !< name, e.g. 'rho'
    character(:), allocatable :: long_name        !< long name, e.g. 'density'
    character(:), allocatable :: units            !< physical units, e.g. 'g/cc'
    character(:), allocatable :: descrip          !< description, e.g. 'Cell-centered density'
    real(rk) :: to_nondim = 1.0_rk                !< scaling factor to convert to non-dimensional form
    real(rk) :: to_dim = 1.0_rk                   !< scaling factor to convert to dimensional form
    integer(ik) :: n_halo_cells = 0               !< # of halo cells used in this block; typically 2 or 3 depending on spatial order
    integer(ik), dimension(2) :: lbounds = 0      !< (i, j); lower bounds (not including halo cells)
    integer(ik), dimension(2) :: ubounds = 0      !< (i, j); upper bounds (not including halo cells)
    integer(ik), dimension(2) :: lbounds_halo = 0 !< (i, j); lower bounds (including halo cells)
    integer(ik), dimension(2) :: ubounds_halo = 0 !< (i, j); upper bounds (including halo cells)
    integer(ik), dimension(2) :: dims = 0         !< (i, j); dimensions
    integer(ik), dimension(4) :: neighbors = 0    !< (left, right, down, up); parallel neighbor tiles/images
    integer(ik) :: edge_size = 0                  !< max number of cells on the edge

    ! Intel compiler alignment hints
    !dir$ attributes align:__ALIGNBYTES__ :: data
  contains
    private

    ! Private methods
    procedure :: assign_field, assign_real_scalar
    procedure :: field_mul_array, field_mul_real, field_mul_field
    procedure :: field_div_real
    procedure :: field_add_field, field_add_real
    procedure :: field_sub_field, field_sub_real

    ! Public methods
    procedure, public :: zero_out_halo
    procedure, public :: gather
    procedure, public :: sync_edges

    generic :: assignment(=) => assign_field, assign_real_scalar
    generic :: operator(+) => field_add_field, field_add_real
    generic :: operator(-) => field_sub_field, field_sub_real
    generic :: operator(*) => field_mul_array, field_mul_real, field_mul_field
    generic :: operator(/) => field_div_real

    ! Finalization
    final :: finalize

  end type field_t

  ! Constructor interface
  interface field_t
    module procedure :: field_constructor
  end interface field_t

contains

  type(field_t) function field_constructor(name, long_name, descrip, units, dims) result(self)
    !< Construct the field_t object
    character(len=*), intent(in) :: name          !< name, e.g. 'rho'
    character(len=*), intent(in) :: long_name     !< long name, e.g. 'density'
    character(len=*), intent(in) :: units         !< physical units, e.g. 'g/cc'
    character(len=*), intent(in) :: descrip       !< description, e.g. 'Cell-centered density'
    integer(ik), dimension(2), intent(in) :: dims !< (i,j); domain size in x and y

    ! Locals
    integer(ik) :: indices(4)

    self%name = name
    self%long_name = long_name
    self%units = units
    self%descrip = descrip
    self%dims = dims

    ! indices = tile_indices(dims)

    ! self%lbounds = indices([1, 3])
    ! self%ubounds = indices([2, 4])
    ! allocate(self%data(self%lbounds(1) - 1:self%ubounds(1) + 1, &
    !                    self%lbounds(2) - 1:self%ubounds(2) + 1))
    ! self%data = 0
    ! self%neighbors = tile_neighbors_2d(periodic=.true.)
    ! self%edge_size = max(self%ubounds(1) - self%lbounds(1) + 1, &
    !                      self%ubounds(2) - self%lbounds(2) + 1)
    ! call co_max(self%edge_size)
  end function field_constructor

  subroutine assign_field(self, f)
    !< Implementation of the (=) operator from another field_t
    class(field_t), intent(in out) :: self
    class(field_t), intent(in) :: f
    call from_field(self, f)
    call self%sync_edges()
  end subroutine assign_field

  pure subroutine assign_real_scalar(self, a)
    !< Implementation of the (=) operator from a real scalar
    class(field_t), intent(in out) :: self
    real(rk), intent(in) :: a
    self%data = a
  end subroutine assign_real_scalar

  pure subroutine from_field(target, source)
    !< Initializes field_t instance target using components
    !< from field_t instance source. Used to initialize a
    !< field_t from another field_t without invoking the
    !< assignment operator.
    type(field_t), intent(in out) :: target
    type(field_t), intent(in) :: source

    target%data = source%data
    target%name = source%name
    target%long_name = source%long_name
    target%units = source%units
    target%descrip = source%descrip
    target%to_nondim = source%to_nondim
    target%to_dim = source%to_dim
    target%n_halo_cells = source%n_halo_cells
    target%lbounds = source%lbounds
    target%ubounds = source%ubounds
    target%lbounds_halo = source%lbounds_halo
    target%ubounds_halo = source%ubounds_halo
    target%dims = source%dims
    target%neighbors = source%neighbors
    target%edge_size = source%edge_size

  end subroutine from_field

  function gather(self, image)
    !< Performs a gather of field data to image.
    class(field_t), intent(in) :: self
    integer(ik), intent(in) :: image
    real(rk) :: gather(self%dims(1), self%dims(2))

    ! real(rk), allocatable :: gather_coarray(:, :)[:]
    ! allocate(gather_coarray(self%dims(1), self%dims(2))[*])
    ! associate(is => self%lbounds(1), ie => self%ubounds(1), &
    !           js => self%lbounds(2), je => self%ubounds(2))
    !   gather_coarray(is:ie, js:je)[image] = self%data(is:ie, js:je)
    !   sync all
    !   if(this_image() == image) gather = gather_coarray
    ! end associate
    ! deallocate(gather_coarray)
  end function gather

  pure type(field_t) function field_add_field(self, f) result(res)
    !< Implementation of the field_t + field_t operation
    class(field_t), intent(in) :: self, f
    call from_field(res, self)
    res%data = self%data + f%data
  end function field_add_field

  pure type(field_t) function field_add_real(self, x) result(res)
    !< Implementation of the field_t + real64 operation
    class(field_t), intent(in) :: self
    real(rk), dimension(:, :), intent(in) :: x
    call from_field(res, self)
    res%data = self%data + x
  end function field_add_real

  pure type(field_t) function field_div_real(self, x) result(res)
    !< Implementation of the field_t / real64 operation
    class(field_t), intent(in) :: self
    real(rk), intent(in) :: x
    call from_field(res, self)
    res%data = self%data / x
  end function field_div_real

  pure type(field_t) function field_mul_array(self, x) result(res)
    !< Implementation of the field_t * array operation
    class(field_t), intent(in) :: self
    real(rk), dimension(:, :), intent(in) :: x
    call from_field(res, self)
    res%data = self%data * x
  end function field_mul_array

  pure type(field_t) function field_mul_real(self, x) result(res)
    !< Implementation of the field_t * real64 operation
    class(field_t), intent(in) :: self
    real(rk), intent(in) :: x
    call from_field(res, self)
    res%data = self%data * x
  end function field_mul_real

  pure type(field_t) function field_mul_field(self, f) result(res)
    !< Implementation of the field_t * field_t operation
    class(field_t), intent(in) :: self, f
    call from_field(res, self)
    res%data = self%data * f%data
  end function field_mul_field

  pure type(field_t) function field_sub_real(self, x) result(res)
    !< Implementation of the field_t - real64 operation
    class(field_t), intent(in) :: self
    real(rk), dimension(:, :), intent(in) :: x
    call from_field(res, self)
    res%data = self%data - x
  end function field_sub_real

  pure type(field_t) function field_sub_field(self, f) result(res)
    !< Implementation of the field_t - field_t operation
    class(field_t), intent(in) :: self, f
    call from_field(res, self)
    res%data = self%data - f%data
  end function field_sub_field

  pure subroutine zero_out_halo(self)
    !< Zero out all of the halo (boundary) cells
    class(field_t), intent(inout) :: self

    associate(ilo_s => self%lbounds_halo(1), ilo_e => self%lbounds(1) - 1, &
              ihi_s => self%ubounds(1) + 1, ihi_e => self%ubounds_halo(1), &
              jlo_s => self%lbounds_halo(2), jlo_e => self%lbounds(2) - 1, &
              jhi_s => self%ubounds(2) + 1, jhi_e => self%ubounds_halo(2))

      self%data(ilo_s:ilo_e, :) = 0.0_rk ! lower i cells
      self%data(ihi_s:ihi_e, :) = 0.0_rk ! upper i cells
      self%data(:, jlo_s:jlo_e) = 0.0_rk ! lower j cells
      self%data(:, jhi_s:jhi_e) = 0.0_rk ! upper j cells
    end associate
  end subroutine zero_out_halo

  subroutine sync_edges(self)
    class(field_t), intent(inout) :: self
    ! real(rk), allocatable, save :: edge(:, :)[:]

    ! associate(is => self%lbounds(1), ie => self%ubounds(1), &
    !           js => self%lbounds(2), je => self%ubounds(2), &
    !           nh => self%n_halo_cells, &
    !           neighbors => self%neighbors)

    !   if(.not. allocated(edge)) allocate(edge(self%edge_size, 4)[*])

    !   sync images(set(neighbors))

    !   ! copy data into coarray buffer
    !   edge(1:je - js + 1, left)[neighbors(left)] = self%data(is, js:je) ! send left
    !   edge(1:je - js + 1, right)[neighbors(right)] = self%data(ie, js:je) ! send right
    !   edge(1:ie - is + 1, down)[neighbors(down)] = self%data(is:ie, js) ! send down
    !   edge(1:ie - is + 1, up)[neighbors(up)] = self%data(is:ie, je) ! send up

    !   sync images(set(neighbors))

    !   ! copy from halo buffer into array
    !   self%data(is - 1, js:je) = edge(1:je - js + 1, right) ! from left
    !   self%data(ie + 1, js:je) = edge(1:je - js + 1, left) ! from right
    !   self%data(is:ie, js - 1) = edge(1:ie - is + 1, up) ! from down
    !   self%data(is:ie, je + 1) = edge(1:ie - is + 1, down) ! from up

    ! end associate

  end subroutine sync_edges

  pure recursive function set(a) result(res)
    integer, intent(in) :: a(:)
    integer, allocatable :: res(:)
    if(size(a) > 1) then
      res = [a(1), set(pack(a(2:),.not. a(2:) == a(1)))]
    else
      res = a
    end if
  end function set

  pure subroutine finalize(self)
    !< Cleanup the field_t object
    type(field_t), intent(inout) :: self
    if(allocated(self%data)) deallocate(self%data)
    if(allocated(self%units)) deallocate(self%units)
    if(allocated(self%name)) deallocate(self%name)
    if(allocated(self%long_name)) deallocate(self%long_name)
  end subroutine finalize
end module mod_field
