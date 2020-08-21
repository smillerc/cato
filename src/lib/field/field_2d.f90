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
  !> Summary: Provide the base 2D field class. This originates from https://github.com/modern-fortran/tsunami
  !> Date: 08/18/2020
  !> Author: Milan Curcic, Sam Miller (minor mods)
  !> Notes:
  !> References:
  !      [1] Milan Curcic, "Modern Fortran: Building efficient parallel applications", 2020

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use, intrinsic :: ieee_arithmetic
  use mod_error, only: error_msg
  use mod_globals, only: enable_debug_print, debug_print, n_ghost_layers
  ! use mod_parallel, only: tile_indices, tile_neighbors_2d

  implicit none

  private
  public :: field_2d_t

  ! Neighbor image indices
  integer(ik), parameter :: left = 1  !< left neighbor
  integer(ik), parameter :: right = 2 !< right neighbor
  integer(ik), parameter :: down = 3  !< neighbor below
  integer(ik), parameter :: up = 4    !< neighbor above

  type :: field_2d_t
    !< A base class that encapsulates 2D data with knowledge of its parallel neighbors and bounds

    private

    ! Core data
    real(rk), allocatable, public :: data(:, :)           !< (i, j); field data

    ! Intel compiler alignment hints
    !dir$ attributes align:__ALIGNBYTES__ :: data

    ! Attributes/metadata
    character(:), allocatable, public :: name      !< name, e.g. 'rho'
    character(:), allocatable, public :: long_name !< long name, e.g. 'density'
    character(:), allocatable, public :: units     !< physical units, e.g. 'g/cc'
    character(:), allocatable, public :: descrip   !< description, e.g. 'Cell-centered density'

    ! Non-dimensionalization factors
    real(rk) :: to_nondim = 1.0_rk  !< scaling factor to convert to non-dimensional form
    real(rk) :: to_dim = 1.0_rk     !< scaling factor to convert to dimensional form
    logical :: is_dim = .true.      !< the current field is in dimensional form
    logical :: is_nondim = .false.  !< the current field is in non-dimensional form

    ! Bounds information
    integer(ik), public :: n_halo_cells = 0               !< # of halo cells used in this block; typically 2 or 3 depending on spatial order
    integer(ik), dimension(2), public :: dims = 0         !< (i, j); dimensions
    integer(ik), dimension(2), public :: lbounds = 0      !< (i, j); lower bounds (not including halo cells)
    integer(ik), dimension(2), public :: ubounds = 0      !< (i, j); upper bounds (not including halo cells)
    integer(ik), dimension(2), public :: lbounds_halo = 0 !< (i, j); lower bounds (including halo cells)
    integer(ik), dimension(2), public :: ubounds_halo = 0 !< (i, j); upper bounds (including halo cells)

    integer(ik), public :: ilo = 0
    integer(ik), public :: ilo_halo = 0
    integer(ik), public :: ihi = 0
    integer(ik), public :: ihi_halo = 0
    integer(ik), public :: jlo = 0
    integer(ik), public :: jlo_halo = 0
    integer(ik), public :: jhi = 0
    integer(ik), public :: jhi_halo = 0

    ! Boundary condition checks
    logical, public :: on_plus_x_bc = .true.  !< does this field live on the +x boundary?
    logical, public :: on_minus_x_bc = .true. !< does this field live on the -x boundary?
    logical, public :: on_plus_y_bc = .true.  !< does this field live on the +y boundary?
    logical, public :: on_minus_y_bc = .true. !< does this field live on the -y boundary?

    ! Parallel neighbor information
    integer(ik), dimension(4) :: neighbors = 0    !< (left, right, down, up); parallel neighbor tiles/images
    integer(ik) :: edge_size = 0                  !< max number of cells on the edge
  contains
    private

    ! Private methods
    procedure, pass(lhs) :: field_add_field, field_sub_field, field_mul_field, field_div_field
    procedure, pass(lhs) :: field_add_real_1d, field_add_real_2d
    procedure, pass(rhs) :: real_1d_add_field, real_2d_add_field
    procedure, pass(lhs) :: field_sub_real_1d, field_sub_real_2d
    procedure, pass(rhs) :: real_1d_sub_field, real_2d_sub_field
    procedure, pass(lhs) :: field_div_real_1d, field_div_real_2d
    procedure, pass(rhs) :: real_1d_div_field, real_2d_div_field
    procedure, pass(lhs) :: field_mul_real_2d, field_mul_real_1d
    procedure, pass(rhs) :: real_1d_mul_field, real_2d_mul_field
    procedure, pass(lhs) :: assign_real_scalar, assign_field
    procedure, pass(self) :: write_field

    ! Public methods
    procedure, public :: make_non_dimensional
    procedure, public :: make_dimensional
    procedure, public :: zero_out_halo
    procedure, public :: has_nans
    procedure, public :: has_negatives
    procedure, public :: gather
    procedure, public :: sync_edges

    generic :: write(formatted) => write_field
    generic, public :: assignment(=) => assign_field, assign_real_scalar
    generic, public :: operator(+) => field_add_field, field_add_real_1d, field_add_real_2d, real_1d_add_field, real_2d_add_field
    generic, public :: operator(-) => field_sub_field, field_sub_real_1d, field_sub_real_2d, real_1d_sub_field, real_2d_sub_field
    generic, public :: operator(*) => field_mul_field, field_mul_real_2d, field_mul_real_1d, real_1d_mul_field, real_2d_mul_field
    generic, public :: operator(/) => field_div_field, field_div_real_1d, field_div_real_2d, real_1d_div_field, real_2d_div_field

    ! Finalization
    final :: finalize

  end type field_2d_t

  ! Constructor interface
  interface field_2d_t
    module procedure :: field_constructor
  end interface field_2d_t

contains

  type(field_2d_t) function field_constructor(name, long_name, descrip, units, dims) result(self)
    !< Construct the field_2d_t object
    character(len=*), intent(in) :: name          !< name, e.g. 'rho'
    character(len=*), intent(in) :: long_name     !< long name, e.g. 'density'
    character(len=*), intent(in) :: units         !< physical units, e.g. 'g/cc'
    character(len=*), intent(in) :: descrip       !< description, e.g. 'Cell-centered density'
    integer(ik), dimension(2), intent(in) :: dims !< (i,j); domain size in x and y

    ! Locals
    integer(ik) :: indices(4)
    integer(ik) :: alloc_status

    if(enable_debug_print) call debug_print('Running field_2d_t%field_constructor() for "'//name//'"', __FILE__, __LINE__)

    self%name = name
    self%long_name = long_name
    self%units = units
    self%descrip = descrip
    self%dims = dims
    self%n_halo_cells = n_ghost_layers

    ! Coarray domain-decomposition
    ! indices = tile_indices(dims)
    ! self%lbounds = indices([1, 3])
    ! self%ubounds = indices([2, 4])

    ! Non-coarray version
    self%lbounds = [1, 1]
    self%ubounds = [dims(1), dims(2)]

    self%lbounds_halo = self%lbounds - self%n_halo_cells
    self%ubounds_halo = self%ubounds + self%n_halo_cells

    self%ilo = 1
    self%jlo = 1
    self%ihi = dims(1)
    self%jhi = dims(2)

    self%ilo_halo = self%ilo - self%n_halo_cells
    self%jlo_halo = self%jlo - self%n_halo_cells
    self%ihi_halo = self%ihi + self%n_halo_cells
    self%jhi_halo = self%jhi + self%n_halo_cells

    associate(ilo => self%lbounds_halo(1), ihi => self%ubounds_halo(1), &
              jlo => self%lbounds_halo(2), jhi => self%ubounds_halo(2))

      ! Allocate w/ halo cells
      allocate(self%data(ilo:ihi, jlo:jhi), stat=alloc_status)
      if(alloc_status /= 0) then
        call error_msg(module_name='mod_fluid', class_name='fluid_t', procedure_name='fluid_mul_real', &
                       message="Unable to allocate field_2d_t%data for '"//self%name//"'", &
                       file_name=__FILE__, line_number=__LINE__)
      end if
    end associate

    self%data = 0.0_rk
    ! self%neighbors = tile_neighbors_2d(periodic=.true.)
    ! self%edge_size = max(self%ubounds(1) - self%lbounds(1) + 1, &
    !                      self%ubounds(2) - self%lbounds(2) + 1)
    ! call co_max(self%edge_size)
  end function field_constructor

  subroutine assign_field(lhs, f)
    !< Implementation of the (=) operator from another field_2d_t
    class(field_2d_t), intent(in out) :: lhs !< left-hand-side of the operation
    class(field_2d_t), intent(in) :: f
    call from_field(lhs, f)
    call lhs%sync_edges()
  end subroutine assign_field

  pure subroutine assign_real_scalar(lhs, a)
    !< Implementation of the (=) operator from a real scalar
    class(field_2d_t), intent(in out) :: lhs !< left-hand-side of the operation
    real(rk), intent(in) :: a
    lhs%data = a
  end subroutine assign_real_scalar

  pure subroutine from_field(target, source)
    !< Initializes field_2d_t instance target using components
    !< from field_2d_t instance source. Used to initialize a
    !< field_2d_t from another field_2d_t without invoking the
    !< assignment operator.
    type(field_2d_t), intent(in out) :: target
    type(field_2d_t), intent(in) :: source

    if(.not. allocated(source%data)) error stop "Error in field_2d_t%from_field(): source%data isn't allocated"
    target%data = source%data

    if(.not. allocated(source%name)) error stop "Error in field_2d_t%from_field(): source%name isn't allocated"
    target%name = source%name

    if(.not. allocated(source%long_name)) error stop "Error in field_2d_t%from_field(): source%long_name isn't allocated"
    target%long_name = source%long_name

    if(.not. allocated(source%units)) error stop "Error in field_2d_t%from_field(): source%units isn't allocated"
    target%units = source%units

    if(.not. allocated(source%descrip)) error stop "Error in field_2d_t%from_field(): source%descrip isn't allocated"
    target%descrip = source%descrip

    target%to_nondim = source%to_nondim
    target%to_dim = source%to_dim
    target%is_dim = source%is_dim
    target%is_nondim = source%is_nondim

    target%on_plus_x_bc = source%on_plus_x_bc
    target%on_minus_x_bc = source%on_minus_x_bc
    target%on_plus_y_bc = source%on_plus_y_bc
    target%on_minus_y_bc = source%on_minus_y_bc

    target%n_halo_cells = source%n_halo_cells
    target%lbounds = source%lbounds
    target%ubounds = source%ubounds

    target%ilo = source%ilo
    target%ihi = source%ihi
    target%jlo = source%jlo
    target%jhi = source%jhi
    target%ilo_halo = source%ilo_halo
    target%ihi_halo = source%ihi_halo
    target%jlo_halo = source%jlo_halo
    target%jhi_halo = source%jhi_halo

    target%lbounds_halo = source%lbounds_halo
    target%ubounds_halo = source%ubounds_halo
    target%dims = source%dims
    target%neighbors = source%neighbors
    target%edge_size = source%edge_size

  end subroutine from_field

  function gather(self, image)
    !< Performs a gather of field data to image.
    class(field_2d_t), intent(in) :: self
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

  pure type(field_2d_t) function field_add_field(lhs, f) result(res)
    !< Implementation of the field_2d_t + field_2d_t operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    class(field_2d_t), intent(in) :: f

    if(.not. field_shapes_match(lhs, f)) then
      error stop "Error in field_2d_t%field_add_field(): field_2d_t shapes don't match"
    end if

    call from_field(res, lhs)
    res%data(res%ilo:res%ihi, res%jlo:res%jhi) = lhs%data(lhs%ilo:lhs%ihi, lhs%jlo:lhs%jhi) + &
                                                 f%data(f%ilo:f%ihi, f%jlo:f%jhi)
  end function field_add_field

  pure type(field_2d_t) function field_sub_field(lhs, f) result(res)
    !< Implementation of the field_2d_t + field_2d_t operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    class(field_2d_t), intent(in) :: f

    if(.not. field_shapes_match(lhs, f)) then
      error stop "Error in field_2d_t%field_sub_field(): field_2d_t shapes don't match"
    end if

    call from_field(res, lhs)
    res%data(res%ilo:res%ihi, res%jlo:res%jhi) = lhs%data(lhs%ilo:lhs%ihi, lhs%jlo:lhs%jhi) - &
                                                 f%data(f%ilo:f%ihi, f%jlo:f%jhi)

  end function field_sub_field

  pure type(field_2d_t) function field_mul_field(lhs, f) result(res)
    !< Implementation of the field_2d_t * field_2d_t operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    class(field_2d_t), intent(in) :: f

    if(.not. field_shapes_match(lhs, f)) then
      error stop "Error in field_2d_t%field_mul_field(): field_2d_t shapes don't match"
    end if

    call from_field(res, lhs)
    res%data(res%ilo:res%ihi, res%jlo:res%jhi) = lhs%data(lhs%ilo:lhs%ihi, lhs%jlo:lhs%jhi) * &
                                                 f%data(f%ilo:f%ihi, f%jlo:f%jhi)
  end function field_mul_field

  pure type(field_2d_t) function field_div_field(lhs, f) result(res)
    !< Implementation of the field_2d_t * field_2d_t operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    class(field_2d_t), intent(in) :: f

    if(.not. field_shapes_match(lhs, f)) then
      error stop "Error in field_2d_t%field_div_field(): field_2d_t shapes don't match"
    end if

    call from_field(res, lhs)
    res%data(res%ilo:res%ihi, res%jlo:res%jhi) = lhs%data(lhs%ilo:lhs%ihi, lhs%jlo:lhs%jhi) / &
                                                 f%data(f%ilo:f%ihi, f%jlo:f%jhi)
  end function field_div_field

  pure type(field_2d_t) function field_add_real_1d(lhs, x) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), intent(in) :: x
    integer(ik) :: ilo, ihi, jlo, jhi

    call from_field(res, lhs)
    res%data = lhs%data + x
  end function field_add_real_1d

  pure type(field_2d_t) function field_add_real_2d(lhs, x) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), dimension(:, :), intent(in) :: x
    integer(ik) :: ilo, ihi, jlo, jhi

    ilo = lbound(x, dim=1) + lhs%n_halo_cells
    ihi = ubound(x, dim=1) - lhs%n_halo_cells
    jlo = lbound(x, dim=2) + lhs%n_halo_cells
    jhi = ubound(x, dim=2) - lhs%n_halo_cells

    if(.not. shapes_match(lhs, x)) then
      error stop "Error in field_2d_t%field_add_real_2d(): field_2d_t and shapes don't match"
    end if
    call from_field(res, lhs)

    res%data(res%ilo:res%ihi, res%jlo:res%jhi) = lhs%data(lhs%ilo:lhs%ihi, lhs%jlo:lhs%jhi) + x(ilo:ihi, jlo:jhi)
  end function field_add_real_2d

  pure type(field_2d_t) function real_1d_add_field(x, rhs) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), intent(in) :: x
    integer(ik) :: ilo, ihi, jlo, jhi

    call from_field(res, rhs)
    res%data = rhs%data + x
  end function real_1d_add_field

  pure type(field_2d_t) function real_2d_add_field(x, rhs) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), dimension(:, :), intent(in) :: x
    integer(ik) :: ilo, ihi, jlo, jhi

    ilo = lbound(x, dim=1) + rhs%n_halo_cells
    ihi = ubound(x, dim=1) - rhs%n_halo_cells
    jlo = lbound(x, dim=2) + rhs%n_halo_cells
    jhi = ubound(x, dim=2) - rhs%n_halo_cells

    if(.not. shapes_match(rhs, x)) then
      error stop "Error in field_2d_t%real_2d_add_field(): field_2d_t and x shapes don't match"
    end if

    call from_field(res, rhs)
    res%data(res%ilo:res%ihi, res%jlo:res%jhi) = rhs%data(rhs%ilo:rhs%ihi, rhs%jlo:rhs%jhi) + x(ilo:ihi, jlo:jhi)
  end function real_2d_add_field

  pure type(field_2d_t) function field_sub_real_1d(lhs, x) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), intent(in) :: x
    call from_field(res, lhs)
    res%data = lhs%data - x
  end function field_sub_real_1d

  pure type(field_2d_t) function field_sub_real_2d(lhs, x) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), dimension(:, :), intent(in) :: x
    integer(ik) :: ilo, ihi, jlo, jhi

    ilo = lbound(x, dim=1) + lhs%n_halo_cells
    ihi = ubound(x, dim=1) - lhs%n_halo_cells
    jlo = lbound(x, dim=2) + lhs%n_halo_cells
    jhi = ubound(x, dim=2) - lhs%n_halo_cells

    if(.not. shapes_match(lhs, x)) then
      error stop "Error in field_2d_t%field_sub_real_2d(): field_2d_t and x shapes don't match"
    end if
    call from_field(res, lhs)

    res%data(res%ilo:res%ihi, res%jlo:res%jhi) = lhs%data(lhs%ilo:lhs%ihi, lhs%jlo:lhs%jhi) - x(ilo:ihi, jlo:jhi)
  end function field_sub_real_2d

  pure type(field_2d_t) function real_1d_sub_field(x, rhs) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), intent(in) :: x
    call from_field(res, rhs)
    res%data = rhs%data - x
  end function real_1d_sub_field

  pure type(field_2d_t) function real_2d_sub_field(x, rhs) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), dimension(:, :), intent(in) :: x
    integer(ik) :: ilo, ihi, jlo, jhi

    ilo = lbound(x, dim=1) + rhs%n_halo_cells
    ihi = ubound(x, dim=1) - rhs%n_halo_cells
    jlo = lbound(x, dim=2) + rhs%n_halo_cells
    jhi = ubound(x, dim=2) - rhs%n_halo_cells

    if(.not. shapes_match(rhs, x)) then
      error stop "Error in field_2d_t%real_2d_sub_field(): field_2d_t and x shapes don't match"
    end if

    call from_field(res, rhs)
    res%data(res%ilo:res%ihi, res%jlo:res%jhi) = rhs%data(rhs%ilo:rhs%ihi, rhs%jlo:rhs%jhi) - x(ilo:ihi, jlo:jhi)
  end function real_2d_sub_field

  pure type(field_2d_t) function field_div_real_1d(lhs, x) result(res)
    !< Implementation of the field_2d_t / real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), intent(in) :: x

    if(abs(x) < tiny(1.0_rk)) then
      error stop "Error in field_2d_t%field_div_real_1d: attempting a divide by 0 (abs(x) < tiny(1.0_rk))"
    end if

    call from_field(res, lhs)
    res%data = lhs%data / x
  end function field_div_real_1d

  pure type(field_2d_t) function field_div_real_2d(lhs, x) result(res)
    !< Implementation of the field_2d_t / real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), dimension(:, :), intent(in) :: x
    integer(ik) :: ilo, ihi, jlo, jhi

    ilo = lbound(x, dim=1) + lhs%n_halo_cells
    ihi = ubound(x, dim=1) - lhs%n_halo_cells
    jlo = lbound(x, dim=2) + lhs%n_halo_cells
    jhi = ubound(x, dim=2) - lhs%n_halo_cells

    if(.not. shapes_match(lhs, x)) then
      error stop "Error in field_2d_t%field_div_real_2d(): field_2d_t and x shapes don't match"
    end if
    call from_field(res, lhs)

    res%data(res%ilo:res%ihi, res%jlo:res%jhi) = lhs%data(lhs%ilo:lhs%ihi, lhs%jlo:lhs%jhi) / x(ilo:ihi, jlo:jhi)
  end function field_div_real_2d

  pure type(field_2d_t) function real_1d_div_field(x, rhs) result(res)
    !< Implementation of the field_2d_t / real64 operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), intent(in) :: x
    call from_field(res, rhs)
    res%data = x / rhs%data
  end function real_1d_div_field

  pure type(field_2d_t) function real_2d_div_field(x, rhs) result(res)
    !< Implementation of the field_2d_t / real64 operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), dimension(:, :), intent(in) :: x
    integer(ik) :: ilo, ihi, jlo, jhi

    ilo = lbound(x, dim=1) + rhs%n_halo_cells
    ihi = ubound(x, dim=1) - rhs%n_halo_cells
    jlo = lbound(x, dim=2) + rhs%n_halo_cells
    jhi = ubound(x, dim=2) - rhs%n_halo_cells

    if(.not. shapes_match(rhs, x)) then
      error stop "Error in field_2d_t%real_2d_div_field(): field_2d_t and x shapes don't match"
    end if
    call from_field(res, rhs)

    res%data(res%ilo:res%ihi, res%jlo:res%jhi) = x(ilo:ihi, jlo:jhi) / rhs%data(rhs%ilo:rhs%ihi, rhs%jlo:rhs%jhi)
  end function real_2d_div_field

  pure type(field_2d_t) function field_mul_real_2d(lhs, x) result(res)
    !< Implementation of the field_2d_t * array operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), dimension(:, :), intent(in) :: x
    integer(ik) :: ilo, ihi, jlo, jhi

    ilo = lbound(x, dim=1) + lhs%n_halo_cells
    ihi = ubound(x, dim=1) - lhs%n_halo_cells
    jlo = lbound(x, dim=2) + lhs%n_halo_cells
    jhi = ubound(x, dim=2) - lhs%n_halo_cells

    if(.not. shapes_match(lhs, x)) then
      error stop "Error in field_2d_t%field_mul_real_2d(): field_2d_t and x shapes don't match"
    end if
    call from_field(res, lhs)

    res%data(res%ilo:res%ihi, res%jlo:res%jhi) = lhs%data(lhs%ilo:lhs%ihi, lhs%jlo:lhs%jhi) * x(ilo:ihi, jlo:jhi)
  end function field_mul_real_2d

  pure type(field_2d_t) function field_mul_real_1d(lhs, x) result(res)
    !< Implementation of the field_2d_t * real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), intent(in) :: x
    call from_field(res, lhs)
    res%data = lhs%data * x
  end function field_mul_real_1d

  pure type(field_2d_t) function real_1d_mul_field(x, rhs) result(res)
    !< Implementation of the real64 * field_2d_t operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), intent(in) :: x
    call from_field(res, rhs)
    res%data = rhs%data * x
  end function real_1d_mul_field

  pure type(field_2d_t) function real_2d_mul_field(x, rhs) result(res)
    !< Implementation of the real64 * field_2d_t operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), dimension(:, :), intent(in) :: x
    integer(ik) :: ilo, ihi, jlo, jhi

    ilo = lbound(x, dim=1) + rhs%n_halo_cells
    ihi = ubound(x, dim=1) - rhs%n_halo_cells
    jlo = lbound(x, dim=2) + rhs%n_halo_cells
    jhi = ubound(x, dim=2) - rhs%n_halo_cells

    if(.not. shapes_match(rhs, x)) then
      error stop "Error in field_2d_t%field_mul_real_2d(): field_2d_t and x shapes don't match"
    end if

    call from_field(res, rhs)
    res%data(res%ilo:res%ihi, res%jlo:res%jhi) = rhs%data(rhs%ilo:rhs%ihi, rhs%jlo:rhs%jhi) * x(ilo:ihi, jlo:jhi)
  end function real_2d_mul_field

  pure subroutine zero_out_halo(self)
    !< Zero out all of the halo (boundary) cells
    class(field_2d_t), intent(inout) :: self

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

  subroutine make_dimensional(self)
    !< Convert to the dimensional form
    class(field_2d_t), intent(inout) :: self

    if(self%is_nondim) then

      if(enable_debug_print) then
        call debug_print('Running field_2d_t%make_dimensional() for "'//self%name//'"', __FILE__, __LINE__)
      end if

      self%data = self%data * self%to_dim
      self%is_nondim = .false.
      self%is_dim = .true.
    end if
  end subroutine make_dimensional

  subroutine make_non_dimensional(self, non_dim_factor)
    !< Convert to the non-dimensional form
    class(field_2d_t), intent(inout) :: self
    real(rk), intent(in) :: non_dim_factor

    if(abs(non_dim_factor) < tiny(1.0_rk)) error stop "Invalid non-dimensionalization factor"

    if(self%is_dim) then

      if(enable_debug_print) then
        call debug_print('Running field_2d_t%make_non_dimensional() for "'//self%name//'"', __FILE__, __LINE__)
      end if

      self%to_nondim = 1.0_rk / non_dim_factor
      self%to_dim = non_dim_factor
      self%data = self%data * self%to_nondim
      self%is_nondim = .true.
      self%is_dim = .false.
    end if
  end subroutine make_non_dimensional

  logical function has_nans(self)
    !< Check for NaNs
    class(field_2d_t), intent(in) :: self

    ! Locals
    integer(ik) :: i, j
    character(len=:), allocatable :: err_message

    has_nans = .false.

    do j = lbound(self%data, dim=2), ubound(self%data, dim=2)
      do i = lbound(self%data, dim=1), ubound(self%data, dim=1)
        if(ieee_is_nan(self%data(i, j))) then
          write(err_message, '(a, i0, ", ", i0, a)') "NaN "//self%name//" found at (", i, j, ")"
          call error_msg(module_name='mod_field', class_name='field_2d_t', procedure_name='check_for_nans', &
                         message=err_message, file_name=__FILE__, line_number=__LINE__, error_stop=.false.)
          has_nans = .true.
        end if
      end do
    end do
  end function has_nans

  logical function has_negatives(self)
    !< Check for negative numbers
    class(field_2d_t), intent(in) :: self

    ! Locals
    integer(ik) :: i, j
    character(len=:), allocatable :: err_message

    has_negatives = .false.

    do j = lbound(self%data, dim=2), ubound(self%data, dim=2)
      do i = lbound(self%data, dim=1), ubound(self%data, dim=1)
        if(self%data(i, j) < 0.0_rk) then
          write(err_message, '(a, i0, ", ", i0, a)') "Negative "//self%name//" found at (", i, j, ")"
          call error_msg(module_name='mod_field', class_name='field_2d_t', procedure_name='check_for_negatives', &
                         message=err_message, file_name=__FILE__, line_number=__LINE__, error_stop=.false.)
          has_negatives = .true.
        end if
      end do
    end do
  end function has_negatives

  pure logical function field_shapes_match(lhs, rhs)
    class(field_2d_t), intent(in) :: lhs
    class(field_2d_t), intent(in) :: rhs
    field_shapes_match = .false.
    if(lhs%dims(1) == rhs%dims(1) .and. lhs%dims(2) == rhs%dims(2)) field_shapes_match = .true.
  end function

  pure logical function shapes_match(lhs, rhs)
    class(field_2d_t), intent(in) :: lhs
    real(rk), dimension(:, :), intent(in) :: rhs
    shapes_match = .false.
    if(lhs%dims(1) == size(rhs, dim=1) .and. lhs%dims(2) == size(rhs, dim=2)) shapes_match = .true.
  end function

  subroutine sync_edges(self)
    class(field_2d_t), intent(inout) :: self
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

  subroutine write_field(self, unit, iotype, v_list, iostat, iomsg)
    !< Implementation of write(*,*) field_2d_t
    class(field_2d_t), intent(in)  :: self      ! Object to write.
    integer, intent(in)            :: unit      ! Internal unit to write to.
    character(*), intent(in)       :: iotype    ! LISTDIRECTED or DTxxx
    integer, intent(in)            :: v_list(:) ! parameters from fmt spec.
    integer, intent(out)           :: iostat    ! non zero on error, etc.
    character(*), intent(inout)    :: iomsg     ! define if iostat non zero.

    write(unit, '(a)') new_line('')
    write(unit, '(a)') "<field_2d_t>: '"//self%name//"'"//new_line('')
    write(unit, '(a)') "long_name: "//self%long_name//new_line('')
    write(unit, '(a)') "units: "//self%units//new_line('')
    write(unit, '(a)') "descrip: "//self%descrip//new_line('')
    write(unit, '(a, l1, a)') "data allocated?: ", allocated(self%data), new_line('')

    if(allocated(self%data)) then
      write(unit, '(a, f0.4, a)') "data minval: ", minval(self%data(self%lbounds(1):self%ubounds(1), self%lbounds(2):self%ubounds(2))), new_line('')
      write(unit, '(a, f0.4, a)') "data maxval: ", maxval(self%data(self%lbounds(1):self%ubounds(1), self%lbounds(2):self%ubounds(2))), new_line('')
    end if

    write(unit, '(a, 2(i0, 1x), a)') "dims: ", self%dims, new_line('')
    write(unit, '(a, 2(i0, 1x), a)') "lbounds: ", self%lbounds, new_line('')
    write(unit, '(a, 2(i0, 1x), a)') "ubounds: ", self%ubounds, new_line('')
    write(unit, '(a, 2(i0, 1x), a)') "lbounds_halo: ", self%lbounds_halo, new_line('')
    write(unit, '(a, 2(i0, 1x), a)') "ubounds_halo: ", self%ubounds_halo, new_line('')
    write(unit, '(a, l1, a)') "On the global +x boundary: ", self%on_plus_x_bc, new_line('')
    write(unit, '(a, l1, a)') "On the global -x boundary: ", self%on_minus_x_bc, new_line('')
    write(unit, '(a, l1, a)') "On the global +y boundary: ", self%on_plus_y_bc, new_line('')
    write(unit, '(a, l1, a)') "On the global -y boundary: ", self%on_minus_y_bc, new_line('')
  end subroutine write_field

  subroutine finalize(self)
    !< Cleanup the field_2d_t object
    type(field_2d_t), intent(inout) :: self

    if(enable_debug_print) call debug_print('Running field_2d_t%finalize() for "'//self%name//'"', __FILE__, __LINE__)

    if(allocated(self%data)) deallocate(self%data)
    if(allocated(self%units)) deallocate(self%units)
    if(allocated(self%name)) deallocate(self%name)
    if(allocated(self%long_name)) deallocate(self%long_name)
  end subroutine finalize
end module mod_field
