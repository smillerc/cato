module mod_field
  !< Summary: Provide the base 2D field class. This originates from https://github.com/modern-fortran/tsunami
  !< Date: 08/18/2020
  !< Author: Milan Curcic, Sam Miller (minor mods)
  !< Notes:
  !< References:
  !      [1] Milan Curcic, "Modern Fortran: Building efficient parallel applications", 2020

  use, intrinsic :: iso_fortran_env, only: ik => int32, rk => real64, std_err => error_unit
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic
  use h5fortran, only: hdf5_file, hsize_t
  use mod_error, only: error_msg
  use mod_parallel, only: tile_indices, tile_neighbors_2d
  use mod_parallel, only: jlo_neighbor => DOWN, &
                          ilo_neighbor => LEFT, &
                          ihi_neighbor => RIGHT, &
                          jhi_neighbor => UP, &
                          ilo_jlo_neighbor => LOWER_LEFT, &
                          ihi_jlo_neighbor => LOWER_RIGHT, &
                          ilo_jhi_neighbor => UPPER_LEFT, &
                          ihi_jhi_neighbor => UPPER_RIGHT
  use mod_globals, only: enable_debug_print, debug_print

  implicit none

  private
  public :: field_2d_t, field_2d, field_shapes_match

  type :: field_2d_t
    !< A base class that encapsulates 2D data with knowledge of its parallel neighbors and bounds

    ! private

    ! Core data
    real(rk), allocatable, public :: data(:, :)   !< (i, j); field data
    ! type(fclDeviceFloat) :: device_data                 !< field data on the OpenCL device (GPU)

    ! Intel compiler alignment hints
    !dir$ attributes align:__ALIGNBYTES__ :: data

    ! Attributes/metadata
    character(:), allocatable, public :: name      !< name, e.g. 'rho'
    character(:), allocatable, public :: long_name !< long name, e.g. 'density'
    character(:), allocatable, public :: units     !< physical units, e.g. 'g/cc'
    character(:), allocatable, public :: description   !< description, e.g. 'Cell-centered density'

    ! Non-dimensionalization factors
    real(rk) :: to_nondim = 1.0_rk  !< scaling factor to convert to non-dimensional form
    real(rk) :: to_dim = 1.0_rk     !< scaling factor to convert to dimensional form
    logical :: is_dim = .true.      !< the current field is in dimensional form
    logical :: is_nondim = .false.  !< the current field is in non-dimensional form

    ! Bounds information
    integer(ik), public :: n_halo_cells = 0               !< # of halo cells used in this block; typically 2 or 3 depending on spatial order
    integer(ik), dimension(2), public :: global_dims = 0         !< (i, j); # of array indices in each dimension
    integer(ik), dimension(2), public :: domain_shape = 0 !< (i, j); # of array indices in each dimension
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
    logical, public :: on_ilo_bc = .false. !< does this field live on the +i boundary?
    logical, public :: on_ihi_bc = .false. !< does this field live on the -i boundary?
    logical, public :: on_jlo_bc = .false. !< does this field live on the +j boundary?
    logical, public :: on_jhi_bc = .false. !< does this field live on the -j boundary?

    ! Sync settings: VERY IMPORTANT, check your algorithms to see if you need this to sync on every
    ! call to =, otherwise be selective about it. This will reduce communication, which is a good thing
    logical :: sync_on_equal = .false. !< Call sync_edges whenever the (=) operator is called (very crucial in some cases)

    ! Parallel neighbor information
    integer(ik), dimension(8) :: neighbors = 0
    integer(ik), allocatable, dimension(:) :: unique_neighbors
    !< (ilo_jlo, jlo, ihi_jlo, ilo, ihi, ilo_jhi, jhi, ihi_jhi); parallel neighbor image indices

    integer(ik) :: host_image_id = 0 !< what image owns this field instance?

    ! OpenCL stuff
    logical :: use_gpu = .false.
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

    ! Public methods
    procedure, public :: make_non_dimensional
    procedure, public :: make_dimensional
    procedure, public :: zero_out_halo
    procedure, public :: check_for_nans
    procedure, public :: check_for_negatives
    procedure, public :: write_field
    procedure, public :: gather
    procedure, public :: sync_edges
    procedure, public :: read_from_h5
    procedure, public :: write_to_h5

    procedure, public :: maxval => field_maxval
    procedure, public :: maxloc => field_maxloc
    procedure, public :: minval => field_minval
    procedure, public :: minloc => field_minloc
    procedure, public :: sum => field_sum

    generic :: write(formatted) => write_field
    generic, public :: assignment(=) => assign_field, assign_real_scalar
    generic, public :: operator(+) => field_add_field, field_add_real_1d, field_add_real_2d, real_1d_add_field, real_2d_add_field
    generic, public :: operator(-) => field_sub_field, field_sub_real_1d, field_sub_real_2d, real_1d_sub_field, real_2d_sub_field
    generic, public :: operator(*) => field_mul_field, field_mul_real_2d, field_mul_real_1d, real_1d_mul_field, real_2d_mul_field
    generic, public :: operator(/) => field_div_field, field_div_real_1d, field_div_real_2d, real_1d_div_field, real_2d_div_field

    ! Finalization
    final :: finalize

  endtype field_2d_t

  interface field_2d
    module procedure :: new_field
  endinterface field_2d

  ! --------------------------------------------------------------------
  ! Arithemetic operators (high level; these farm it out to cpu or gpu)
  ! --------------------------------------------------------------------
  interface
    ! CPU operators
    module subroutine field_add_field_cpu(lhs, f, res)
      !< Implementation of the field_2d_t + field_2d_t operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      class(field_2d_t), intent(in) :: f
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_add_field_cpu

    module subroutine field_sub_field_cpu(lhs, f, res)
      !< Implementation of the field_2d_t + field_2d_t operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      class(field_2d_t), intent(in) :: f
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_sub_field_cpu

    module subroutine field_mul_field_cpu(lhs, f, res)
      !< Implementation of the field_2d_t * field_2d_t operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      class(field_2d_t), intent(in) :: f
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_mul_field_cpu

    module subroutine field_div_field_cpu(lhs, f, res)
      !< Implementation of the field_2d_t * field_2d_t operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      class(field_2d_t), intent(in) :: f
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_div_field_cpu

    module subroutine field_add_real_1d_cpu(lhs, x, res)
      !< Implementation of the field_2d_t + real64 operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_add_real_1d_cpu

    module subroutine field_add_real_2d_cpu(lhs, x, res)
      !< Implementation of the field_2d_t + real64 operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_add_real_2d_cpu

    module subroutine field_sub_real_1d_cpu(lhs, x, res)
      !< Implementation of the field_2d_t + real64 operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_sub_real_1d_cpu

    module subroutine field_sub_real_2d_cpu(lhs, x, res)
      !< Implementation of the field_2d_t + real64 operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_sub_real_2d_cpu

    module subroutine real_1d_sub_field_cpu(x, rhs, res)
      !< Implementation of the field_2d_t + real64 operation
      class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
      real(rk), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine real_1d_sub_field_cpu

    module subroutine real_2d_sub_field_cpu(x, rhs, res)
      !< Implementation of the field_2d_t + real64 operation
      class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
      real(rk), dimension(rhs%ilo:, rhs%jlo:), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine real_2d_sub_field_cpu

    module subroutine field_div_real_1d_cpu(lhs, x, res)
      !< Implementation of the field_2d_t / real64 operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_div_real_1d_cpu

    module subroutine field_div_real_2d_cpu(lhs, x, res)
      !< Implementation of the field_2d_t / real64 operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_div_real_2d_cpu

    module subroutine real_1d_div_field_cpu(x, rhs, res)
      !< Implementation of the field_2d_t / real64 operation
      class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
      real(rk), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine real_1d_div_field_cpu

    module subroutine real_2d_div_field_cpu(x, rhs, res)
      !< Implementation of the field_2d_t / real64 operation
      class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
      real(rk), dimension(rhs%ilo:, rhs%jlo:), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine real_2d_div_field_cpu

    module subroutine field_mul_real_2d_cpu(lhs, x, res)
      !< Implementation of the field_2d_t * array operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_mul_real_2d_cpu

    module subroutine field_mul_real_1d_cpu(lhs, x, res)
      !< Implementation of the field_2d_t * real64 operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_mul_real_1d_cpu

    module function field_maxval_cpu(f) result(res)
      class(field_2d_t), intent(in) :: f
      real(rk) :: res
    endfunction field_maxval_cpu

    module function field_maxloc_cpu(f) result(res)
      class(field_2d_t), intent(in) :: f
      integer(ik), dimension(2) :: res
    endfunction field_maxloc_cpu

    module function field_minval_cpu(f) result(res)
      class(field_2d_t), intent(in) :: f
      real(rk) :: res
    endfunction field_minval_cpu

    module function field_minloc_cpu(f) result(res)
      class(field_2d_t), intent(in) :: f
      integer(ik), dimension(2) :: res
    endfunction field_minloc_cpu

    module function field_sum_cpu(f) result(res)
      class(field_2d_t), intent(in) :: f
      real(rk) :: res
    endfunction field_sum_cpu

    ! GPU operators
    module subroutine field_add_field_gpu(lhs, f, res)
      !< Implementation of the field_2d_t + field_2d_t operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      class(field_2d_t), intent(in) :: f
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_add_field_gpu

    module subroutine field_sub_field_gpu(lhs, f, res)
      !< Implementation of the field_2d_t + field_2d_t operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      class(field_2d_t), intent(in) :: f
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_sub_field_gpu

    module subroutine field_mul_field_gpu(lhs, f, res)
      !< Implementation of the field_2d_t * field_2d_t operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      class(field_2d_t), intent(in) :: f
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_mul_field_gpu

    module subroutine field_div_field_gpu(lhs, f, res)
      !< Implementation of the field_2d_t * field_2d_t operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      class(field_2d_t), intent(in) :: f
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_div_field_gpu

    module subroutine field_add_real_1d_gpu(lhs, x, res)
      !< Implementation of the field_2d_t + real64 operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_add_real_1d_gpu

    module subroutine field_add_real_2d_gpu(lhs, x, res)
      !< Implementation of the field_2d_t + real64 operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_add_real_2d_gpu

    module subroutine field_sub_real_1d_gpu(lhs, x, res)
      !< Implementation of the field_2d_t + real64 operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_sub_real_1d_gpu

    module subroutine field_sub_real_2d_gpu(lhs, x, res)
      !< Implementation of the field_2d_t + real64 operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_sub_real_2d_gpu

    module subroutine real_1d_sub_field_gpu(x, rhs, res)
      !< Implementation of the field_2d_t + real64 operation
      class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
      real(rk), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine real_1d_sub_field_gpu

    module subroutine real_2d_sub_field_gpu(x, rhs, res)
      !< Implementation of the field_2d_t + real64 operation
      class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
      real(rk), dimension(rhs%ilo:, rhs%jlo:), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine real_2d_sub_field_gpu

    module subroutine field_div_real_1d_gpu(lhs, x, res)
      !< Implementation of the field_2d_t / real64 operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_div_real_1d_gpu

    module subroutine field_div_real_2d_gpu(lhs, x, res)
      !< Implementation of the field_2d_t / real64 operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_div_real_2d_gpu

    module subroutine real_1d_div_field_gpu(x, rhs, res)
      !< Implementation of the field_2d_t / real64 operation
      class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
      real(rk), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine real_1d_div_field_gpu

    module subroutine real_2d_div_field_gpu(x, rhs, res)
      !< Implementation of the field_2d_t / real64 operation
      class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
      real(rk), dimension(rhs%ilo:, rhs%jlo:), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine real_2d_div_field_gpu

    module subroutine field_mul_real_2d_gpu(lhs, x, res)
      !< Implementation of the field_2d_t * array operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_mul_real_2d_gpu

    module subroutine field_mul_real_1d_gpu(lhs, x, res)
      !< Implementation of the field_2d_t * real64 operation
      class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
      real(rk), intent(in) :: x
      type(field_2d_t), intent(inout) :: res
    endsubroutine field_mul_real_1d_gpu

    module function field_maxval_gpu(f) result(res)
      class(field_2d_t), intent(in) :: f
      real(rk) :: res
    endfunction field_maxval_gpu

    module function field_maxloc_gpu(f) result(res)
      class(field_2d_t), intent(in) :: f
      integer(ik), dimension(2) :: res
    endfunction field_maxloc_gpu

    module function field_minval_gpu(f) result(res)
      class(field_2d_t), intent(in) :: f
      real(rk) :: res
    endfunction field_minval_gpu

    module function field_minloc_gpu(f) result(res)
      class(field_2d_t), intent(in) :: f
      integer(ik), dimension(2) :: res
    endfunction field_minloc_gpu

    module function field_sum_gpu(f) result(res)
      class(field_2d_t), intent(in) :: f
      real(rk) :: res
    endfunction field_sum_gpu
  endinterface

contains

  type(field_2d_t) function new_field(name, long_name, descrip, units, global_dims, n_halo_cells) result(self)
    !< Construct the field_2d_t object
    character(len=*), intent(in) :: name          !< name, e.g. 'rho'
    character(len=*), intent(in) :: long_name     !< long name, e.g. 'density'
    character(len=*), intent(in) :: units         !< physical units, e.g. 'g/cc'
    character(len=*), intent(in) :: descrip       !< description, e.g. 'Cell-centered density'
    integer(ik), dimension(2), intent(in) :: global_dims !< (i,j); domain size in x and y (not including halo)
    integer(ik), intent(in) :: n_halo_cells

    ! Locals
    integer(ik) :: indices(4)
    integer(ik) :: alloc_status

    if(enable_debug_print) call debug_print('Running field_2d_t%new_field() for "'//name//'"', __FILE__, __LINE__)

    self%name = name
    self%long_name = long_name
    self%units = units
    self%description = descrip
    self%global_dims = global_dims
    self%n_halo_cells = n_halo_cells

    ! Coarray domain-decomposition. The indices are tiled based on the number of available images
    indices = tile_indices(global_dims)
    self%lbounds = indices([1, 3])
    self%ubounds = indices([2, 4])

    ! The domain shape is exent in each direction
    self%domain_shape(1) = self%ubounds(1) - self%lbounds(1) + 1 ! # of cells in the i direction
    self%domain_shape(2) = self%ubounds(2) - self%lbounds(2) + 1 ! # of cells in the j direction
    self%host_image_id = this_image()

    ! Determine if this field (or subdomain) has an edge on one of the boundaries
    ! When data is synced between images / subdomains, those on the boundary must
    ! have the procedure done via the boundary condition classes. Those on the interior
    ! of the domain just exchange info with their neighbors
    if(self%lbounds(1) == 1) self%on_ilo_bc = .true.
    if(self%lbounds(2) == 1) self%on_jlo_bc = .true.
    if(self%ubounds(1) == global_dims(1)) self%on_ihi_bc = .true.
    if(self%ubounds(2) == global_dims(2)) self%on_jhi_bc = .true.

    self%lbounds_halo = self%lbounds - self%n_halo_cells
    self%ubounds_halo = self%ubounds + self%n_halo_cells

    self%ilo = self%lbounds(1)
    self%jlo = self%lbounds(2)
    self%ihi = self%ubounds(1)
    self%jhi = self%ubounds(2)

    self%ilo_halo = self%ilo - self%n_halo_cells
    self%jlo_halo = self%jlo - self%n_halo_cells
    self%ihi_halo = self%ihi + self%n_halo_cells
    self%jhi_halo = self%jhi + self%n_halo_cells

    associate(ilo => self%lbounds_halo(1), ihi => self%ubounds_halo(1), &
              jlo => self%lbounds_halo(2), jhi => self%ubounds_halo(2))

      ! Allocate w/ halo cells
      allocate(self%data(ilo:ihi, jlo:jhi), stat=alloc_status)
      if(alloc_status /= 0) error stop "Unable to allocate field%data"
    endassociate

    self%data = 0.0_rk
    self%neighbors = tile_neighbors_2d(is_periodic=.true.)
    self%unique_neighbors = set(self%neighbors)
  endfunction new_field

  function gather(self, image, dimensionalize, unit_conversion_factor)
    !< Performs a gather of field data to image.
    class(field_2d_t), intent(in) :: self
    integer(ik), intent(in) :: image                         !< Which image to gather to, e.g. 1 ?
    logical, intent(in), optional :: dimensionalize          !< Apply dimensionalization before gathering?
    real(rk), intent(in), optional :: unit_conversion_factor ! < Apply unit conversion?

    integer(ik) :: ilo, ihi, jlo, jhi
    real(rk), dimension(self%global_dims(1), self%global_dims(2)) :: gather
    real(rk), dimension(self%lbounds(1):self%ubounds(1), self%lbounds(2):self%ubounds(2)) :: dimensional_data
    real(rk), allocatable :: gather_coarray(:, :)[:]

    if(.not. allocated(gather_coarray)) allocate(gather_coarray(self%global_dims(1), self%global_dims(2))[*])

    ilo = self%lbounds(1)
    ihi = self%ubounds(1)
    jlo = self%lbounds(2)
    jhi = self%ubounds(2)

    if(present(dimensionalize)) then
      if(dimensionalize) then
        dimensional_data = self%data(ilo:ihi, jlo:jhi) * self%to_dim
      else
        dimensional_data = self%data(ilo:ihi, jlo:jhi)
      endif
    else
      dimensional_data = self%data(ilo:ihi, jlo:jhi)
    endif

    if(present(unit_conversion_factor)) dimensional_data = dimensional_data * unit_conversion_factor

    gather_coarray(ilo:ihi, jlo:jhi)[image] = dimensional_data(ilo:ihi, jlo:jhi)

    sync all

    if(this_image() == image) gather = gather_coarray
    deallocate(gather_coarray)
  endfunction gather

  pure logical function field_shapes_match(lhs, rhs)
    class(field_2d_t), intent(in) :: lhs
    class(field_2d_t), intent(in) :: rhs
    field_shapes_match = .false.

    if(lhs%domain_shape(1) == rhs%domain_shape(1) .and. &
       lhs%domain_shape(2) == rhs%domain_shape(2)) then
      field_shapes_match = .true.
    endif
  endfunction field_shapes_match

  logical function shapes_match(lhs, rhs)
    class(field_2d_t), intent(in) :: lhs
    real(rk), dimension(:, :), intent(in) :: rhs
    shapes_match = .false.

    if(lhs%domain_shape(1) == size(rhs, dim=1) .and. &
       lhs%domain_shape(2) == size(rhs, dim=2)) then
      shapes_match = .true.
    endif
  endfunction shapes_match

  subroutine sync_edges(self)
    !< Syncronize the edges, e.g. exchange halo data, with neighboring fields / sub-domains. This will NOT sync the edge of
    !< a field if that edge is on the boundary. The boundary halo data is handled separately by boundary condition types.

    class(field_2d_t), intent(inout) :: self

    integer(ik) :: sync_stat           !< syncronization status
    character(len=200) :: sync_err_msg !< syncronization error message (if any)
    real(rk), dimension(:, :), allocatable, save :: ilo_edge[:]   !< (i, j)[image]; Coarray buffer to copy left neighbor halo data
    real(rk), dimension(:, :), allocatable, save :: ihi_edge[:]  !< (i, j)[image]; Coarray buffer to copy right neighbor halo data
    real(rk), dimension(:, :), allocatable, save :: jhi_edge[:]    !< (i, j)[image]; Coarray buffer to copy top neighbor halo data
    real(rk), dimension(:, :), allocatable, save :: jlo_edge[:] !< (i, j)[image]; Coarray buffer to copy bottom neighbor halo data

    real(rk), dimension(:, :), allocatable, save :: ilo_jhi_corner[:]  !< (i, j)[image]; Coarray buffer to copy bottom neighbor halo data
    real(rk), dimension(:, :), allocatable, save :: ihi_jhi_corner[:] !< (i, j)[image]; Coarray buffer to copy bottom neighbor halo data
    real(rk), dimension(:, :), allocatable, save :: ilo_jlo_corner[:]  !< (i, j)[image]; Coarray buffer to copy bottom neighbor halo data
    real(rk), dimension(:, :), allocatable, save :: ihi_jlo_corner[:] !< (i, j)[image]; Coarray buffer to copy bottom neighbor halo data

    integer(ik) :: ilo, ihi, jlo, jhi
    integer(ik) :: ilo_halo, ihi_halo, jlo_halo, jhi_halo
    integer(ik) :: nh, ni, nj
    if(enable_debug_print) call debug_print('Running field_2d_t()%sync_edges ', __FILE__, __LINE__)
    sync_stat = 0
    sync_err_msg = ''

    ! Only allocate once, b/c this will cause an implicit sync all due to the coarray index
    ilo = self%lbounds(1)
    ihi = self%ubounds(1)
    jlo = self%lbounds(2)
    jhi = self%ubounds(2)
    ilo_halo = self%lbounds_halo(1)
    ihi_halo = self%ubounds_halo(1)
    jlo_halo = self%lbounds_halo(2)
    jhi_halo = self%ubounds_halo(2)
    nh = self%n_halo_cells
    ni = self%domain_shape(1)
    nj = self%domain_shape(2)

    if(.not. allocated(ilo_edge)) allocate(ilo_edge(nh, nj)[*])
    if(.not. allocated(ihi_edge)) allocate(ihi_edge(nh, nj)[*])
    if(.not. allocated(jhi_edge)) allocate(jhi_edge(ni, nh)[*])
    if(.not. allocated(jlo_edge)) allocate(jlo_edge(ni, nh)[*])

    if(.not. allocated(ilo_jhi_corner)) allocate(ilo_jhi_corner(nh, nh)[*])
    if(.not. allocated(ilo_jlo_corner)) allocate(ilo_jlo_corner(nh, nh)[*])
    if(.not. allocated(ihi_jhi_corner)) allocate(ihi_jhi_corner(nh, nh)[*])
    if(.not. allocated(ihi_jlo_corner)) allocate(ihi_jlo_corner(nh, nh)[*])

    ilo_jhi_corner = 0.0_rk
    ihi_jhi_corner = 0.0_rk
    ilo_jlo_corner = 0.0_rk
    ihi_jlo_corner = 0.0_rk
    ilo_edge = 0.0_rk
    ihi_edge = 0.0_rk
    jhi_edge = 0.0_rk
    jlo_edge = 0.0_rk

    ! These were originally inside of an associate block, but testing shows that associate blocks
    ! and coarray syntax don't always give the right behabior at runtime

    ! sync images(set(self%neighbors), stat=sync_stat, errmsg=sync_err_msg)
    sync images(self%unique_neighbors, stat=sync_stat, errmsg=sync_err_msg)
    if(sync_stat /= 0) then
      write(std_err, '(a, i0, a)') "Error: unable to sync images in "//__FILE__//":", &
        __LINE__, " sync_err_msg: '"//trim(sync_err_msg)//"'"
      error stop "Unable to sync images, see standard error unit for more information"
    endif

    ! Copy data into the coarray buffer. Data is not copied if the field is on the boundary, b/c this
    ! transfer of data must be handled by boundary condition classes
    ! We are transfering the cells from inside the real domain onto the halo cells of the neighbor domain. Since
    ! the # of halo cells >= 1, we need to account for variable sizes.

    ! send the ihi edge data to the ilo edge of the neighbor on the ihi side
    ilo_edge(:, :)[self%neighbors(ihi_neighbor)] = self%data(ihi - nh + 1:ihi, jlo:jhi)
    ihi_edge(:, :)[self%neighbors(ilo_neighbor)] = self%data(ilo:ilo + nh - 1, jlo:jhi)
    jhi_edge(:, :)[self%neighbors(jlo_neighbor)] = self%data(ilo:ihi, jlo:jlo + nh - 1)
    jlo_edge(:, :)[self%neighbors(jhi_neighbor)] = self%data(ilo:ihi, jhi - nh + 1:jhi)

    ilo_jhi_corner(:, :)[self%neighbors(ihi_jlo_neighbor)] = self%data(ihi - nh + 1:ihi, jlo:jlo + nh - 1)
    ilo_jlo_corner(:, :)[self%neighbors(ihi_jhi_neighbor)] = self%data(ihi - nh + 1:ihi, jhi - nh + 1:jhi)
    ihi_jhi_corner(:, :)[self%neighbors(ilo_jlo_neighbor)] = self%data(ilo:ilo + nh - 1, jlo:jlo + nh - 1)
    ihi_jlo_corner(:, :)[self%neighbors(ilo_jhi_neighbor)] = self%data(ilo:ilo + nh - 1, jhi - nh + 1:jhi)

    ! sync images(set(self%neighbors), stat=sync_stat, errmsg=sync_err_msg)
    sync images(self%unique_neighbors, stat=sync_stat, errmsg=sync_err_msg)
    if(sync_stat /= 0) then
      write(std_err, '(a, i0, a)') "Error: unable to sync images in "//__FILE__//":", &
        __LINE__, " sync_err_msg: '"//trim(sync_err_msg)//"'"
      error stop "Unable to sync images, see standard error unit for more information"
    endif

    ! Now copy the edge data to the halo cells of the current image
    if(.not. self%on_ilo_bc) self%data(ilo_halo:ilo - 1, jlo:jhi) = ilo_edge
    if(.not. self%on_ihi_bc) self%data(ihi + 1:ihi_halo, jlo:jhi) = ihi_edge
    if(.not. self%on_jlo_bc) self%data(ilo:ihi, jlo_halo:jlo - 1) = jlo_edge
    if(.not. self%on_jhi_bc) self%data(ilo:ihi, jhi + 1:jhi_halo) = jhi_edge

    if(.not. self%on_ihi_bc .and. .not. self%on_jhi_bc) self%data(ihi + 1:ihi_halo, jhi + 1:jhi_halo) = ihi_jhi_corner
    if(.not. self%on_ilo_bc .and. .not. self%on_jhi_bc) self%data(ilo_halo:ilo - 1, jhi + 1:jhi_halo) = ilo_jhi_corner
    if(.not. self%on_ihi_bc .and. .not. self%on_jlo_bc) self%data(ihi + 1:ihi_halo, jlo_halo:jlo - 1) = ihi_jlo_corner
    if(.not. self%on_ilo_bc .and. .not. self%on_jlo_bc) self%data(ilo_halo:ilo - 1, jlo_halo:jlo - 1) = ilo_jlo_corner

  endsubroutine sync_edges

  pure recursive function set(a) result(res)
    integer, intent(in) :: a(:)
    integer, allocatable :: res(:)
    if(size(a) > 1) then
      res = [a(1), set(pack(a(2:),.not. a(2:) == a(1)))]
    else
      res = a
    endif
  endfunction set

  subroutine write_field(self, unit, iotype, v_list, iostat, iomsg)
    !< Implementation of write(*,*) field_2d_t
    class(field_2d_t), intent(in)  :: self      ! Object to write.
    integer, intent(in)            :: unit      ! Internal unit to write to.
    character(*), intent(in)       :: iotype    ! LISTDIRECTED or DTxxx
    integer, intent(in)            :: v_list(:) ! parameters from fmt spec.
    integer, intent(out)           :: iostat    ! non zero on error, etc.
    character(*), intent(inout)    :: iomsg     ! define if iostat non zero.

    character(len=1), parameter :: nl = new_line('')

    iostat = 0
    write(unit, '(a)') nl
    write(unit, '(a)') "<field_2d_t>: '"//self%name//"'"//nl
    write(unit, '(a)') "long_name: "//self%long_name//nl
    write(unit, '(a)') "units: "//self%units//nl
    write(unit, '(a)') "description: "//self%description//nl
    write(unit, '(a, i0, a)') "host image: ", self%host_image_id, nl

    write(unit, '(a, l1, a)') "data allocated?: ", allocated(self%data), nl

    if(allocated(self%data)) then
      write(unit, '(a, f0.4, a)') "data minval: ", minval(self%data(self%lbounds(1):self%ubounds(1), &
                                                                    self%lbounds(2):self%ubounds(2))), nl
      write(unit, '(a, f0.4, a)') "data maxval: ", maxval(self%data(self%lbounds(1):self%ubounds(1), &
                                                                    self%lbounds(2):self%ubounds(2))), nl
    endif

    write(unit, '(a, 2(i0, 1x), a)') "global_dims: ", self%global_dims, nl
    write(unit, '(a, 2(i0, 1x), a)') "domain_shape: ", self%domain_shape, nl
    write(unit, '(a, 2(i0, 1x), a)') "data_shape: ", shape(self%data), nl
    write(unit, '(a, i0, a)') "left neighbor : ", self%neighbors(ilo_neighbor), nl
    write(unit, '(a, i0, a)') "right neighbor: ", self%neighbors(ihi_neighbor), nl
    write(unit, '(a, i0, a)') "down neighbor : ", self%neighbors(jlo_neighbor), nl
    write(unit, '(a, i0, a)') "up neighbor   : ", self%neighbors(jhi_neighbor), nl

    write(unit, '(a, 2(i0, 1x), a)') "lbounds: ", self%lbounds, nl
    write(unit, '(a, 2(i0, 1x), a)') "ubounds: ", self%ubounds, nl
    write(unit, '(a, 2(i0, 1x), a)') "lbounds_halo: ", self%lbounds_halo, nl
    write(unit, '(a, 2(i0, 1x), a)') "ubounds_halo: ", self%ubounds_halo, nl
    write(unit, '(a, l1, a)') "On the global +x boundary: ", self%on_ihi_bc, nl
    write(unit, '(a, l1, a)') "On the global -x boundary: ", self%on_ilo_bc, nl
    write(unit, '(a, l1, a)') "On the global +y boundary: ", self%on_jhi_bc, nl
    write(unit, '(a, l1, a)') "On the global -y boundary: ", self%on_jlo_bc, nl
  endsubroutine write_field

  ! --------------------------------------------------------------------
  ! I/O for HDF5
  ! --------------------------------------------------------------------
  subroutine read_from_h5(self, filename, dataset)
    !< Read in the data from an hdf5 file. This will read from the
    !< file and only grab the per-image data, e.g. each grid block will
    !< only read from the indices it has been assigned with respect
    !< to the global domain.
    class(field_2d_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dataset

    type(hdf5_file) :: h5
    integer(ik) :: alloc_status, ilo, ihi, jlo, jhi
    logical :: file_exists
    character(32) :: str_buff = ''
    character(300) :: msg = ''
    real(rk), allocatable, dimension(:, :) :: data

    integer(hsize_t), allocatable :: dims(:)

    file_exists = .false.
    inquire(file=filename, exist=file_exists)

    if(.not. file_exists) then
      call error_msg(module_name='mod_field', class_name='field_2d_t', &
                     procedure_name='read_from_h5', &
                     message='file not found: "'//filename//'"', &
                     file_name=__FILE__, line_number=__LINE__)
    endif

    call h5%initialize(filename=trim(filename), status='old', action='r')
    ! I couldn't get the slice read to work, so we read the whole array
    ! in and extract the slice later
    call h5%shape(dataset, dims)
    allocate(data(dims(1), dims(2)))
    call h5%read(dname=dataset, value=data)
    self%data = data(self%lbounds_halo(1) + self%n_halo_cells:self%ubounds_halo(1) + self%n_halo_cells, &
                     self%lbounds_halo(2) + self%n_halo_cells:self%ubounds_halo(2) + self%n_halo_cells)
    call h5%finalize()
  endsubroutine read_from_h5

  subroutine write_to_h5(self, filename, dataset)
    !< Write the global data to an hdf5 file. This gathers all to image 1 and
    !< writes to a single file. Technically, this could be done in parallel
    !< in the future, but this is the simple case
    class(field_2d_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: dataset
    type(hdf5_file) :: h5

    ! Gather all to the master image and write to file
    call h5%initialize(filename=trim(filename), status='new', action='w')
    call h5%write(dname=dataset, value=self%gather(image=1))
    call h5%finalize()
  endsubroutine write_to_h5

  ! --------------------------------------------------------------------
  ! Utitlities
  ! --------------------------------------------------------------------
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
    endassociate
  endsubroutine zero_out_halo

  subroutine make_dimensional(self)
    !< Convert to the dimensional form
    class(field_2d_t), intent(inout) :: self

    if(self%is_nondim) then

      if(enable_debug_print) then
        call debug_print('Running field_2d_t%make_dimensional() for "'//self%name//'"', __FILE__, __LINE__)
      endif

      self%data = self%data * self%to_dim
      self%is_nondim = .false.
      self%is_dim = .true.
    endif
  endsubroutine make_dimensional

  subroutine make_non_dimensional(self, non_dim_factor)
    !< Convert to the non-dimensional form
    class(field_2d_t), intent(inout) :: self
    real(rk), intent(in) :: non_dim_factor

    if(abs(non_dim_factor) < tiny(1.0_rk)) error stop "Invalid non-dimensionalization factor"

    if(self%is_dim) then

      if(enable_debug_print) then
        call debug_print('Running field_2d_t%make_non_dimensional() for "'//self%name//'"', __FILE__, __LINE__)
      endif

      self%to_nondim = non_dim_factor
      self%to_dim = 1.0_rk / non_dim_factor
      self%data = self%data * self%to_nondim
      self%is_nondim = .true.
      self%is_dim = .false.
    endif
  endsubroutine make_non_dimensional

  subroutine check_for_nans(self)
    !< Check for NaNs
    class(field_2d_t), intent(in) :: self

    ! Locals
    integer(ik) :: i, j
    character(len=128) :: err_message = ''

    if(enable_debug_print) call debug_print('Running "'//trim(self%name)//'" field_2d_t%has_nans() ', __FILE__, __LINE__)

    do j = self%lbounds(2), self%ubounds(2)
      do i = self%lbounds(1), self%ubounds(1)
        if(ieee_is_nan(self%data(i, j))) then
          write(err_message, '(a, i0, ", ", i0, a)') "NaN "//self%name//" found at (", i, j, ")"
          call error_msg(module_name='mod_field', class_name='field_2d_t', procedure_name='check_for_nans', &
                         message=trim(err_message), file_name=__FILE__, line_number=__LINE__, error_stop=.true.)
        endif
      enddo
    enddo
  endsubroutine check_for_nans

  subroutine check_for_negatives(self)
    !< Check for negative numbers
    class(field_2d_t), intent(in) :: self

    ! Locals
    integer(ik) :: i, j
    character(len=128) :: err_message = ''

    if(enable_debug_print) then
      call debug_print('Running "'//trim(self%name)//'" field_2d_t%check_for_negatives() ', __FILE__, __LINE__)
    endif

    do j = self%lbounds(2), self%ubounds(2)
      do i = self%lbounds(1), self%ubounds(1)
        if(self%data(i, j) < 0.0_rk) then
          write(err_message, '(a, i0, ", ", i0, a, es16.6)') "Negative "//self%name//" found at (", i, j, ")", self%data(i, j)
          call error_msg(module_name='mod_field', class_name='field_2d_t', procedure_name='check_for_negatives', &
                         message=trim(err_message), file_name=__FILE__, line_number=__LINE__, error_stop=.true.)
        endif
      enddo
    enddo
  endsubroutine check_for_negatives

  ! --------------------------------------------------------------------
  ! Clear-up all the memory
  ! --------------------------------------------------------------------
  subroutine finalize(self)
    !< Cleanup the field_2d_t object
    type(field_2d_t), intent(inout) :: self

    if(enable_debug_print) call debug_print('Running "'//trim(self%name)//'" field_2d_t%finalize()', __FILE__, __LINE__)

    if(allocated(self%data)) deallocate(self%data)
    if(allocated(self%unique_neighbors)) deallocate(self%unique_neighbors)
    if(allocated(self%units)) deallocate(self%units)
    if(allocated(self%name)) deallocate(self%name)
    if(allocated(self%long_name)) deallocate(self%long_name)
  endsubroutine finalize

  ! --------------------------------------------------------------------
  ! Initializer operators
  ! --------------------------------------------------------------------
  pure subroutine from_field(new, old)
    !< Initializes field_2d_t instance new using components
    !< from field_2d_t instance old. Used to initialize a
    !< field_2d_t from another field_2d_t without invoking the
    !< assignment operator.
    type(field_2d_t), intent(inout) :: new
    type(field_2d_t), intent(in) :: old

    new%data = old%data
    new%name = old%name
    new%long_name = old%long_name
    new%units = old%units
    new%description = old%description
    new%to_nondim = old%to_nondim
    new%to_dim = old%to_dim
    new%is_dim = old%is_dim
    new%is_nondim = old%is_nondim
    new%n_halo_cells = old%n_halo_cells
    new%global_dims = old%global_dims
    new%domain_shape = old%domain_shape
    new%lbounds = old%lbounds
    new%ubounds = old%ubounds
    new%lbounds_halo = old%lbounds_halo
    new%ubounds_halo = old%ubounds_halo
    new%ilo = old%ilo
    new%ihi = old%ihi
    new%jlo = old%jlo
    new%jhi = old%jhi
    new%ihi_halo = old%ihi_halo
    new%ilo_halo = old%ilo_halo
    new%jlo_halo = old%jlo_halo
    new%jhi_halo = old%jhi_halo
    new%on_ilo_bc = old%on_ilo_bc
    new%on_ihi_bc = old%on_ihi_bc
    new%on_jlo_bc = old%on_jlo_bc
    new%on_jhi_bc = old%on_jhi_bc
    new%neighbors = old%neighbors
    new%unique_neighbors = old%unique_neighbors
    new%host_image_id = old%host_image_id
  endsubroutine from_field

  ! --------------------------------------------------------------------
  ! Assignment operators
  ! --------------------------------------------------------------------
  subroutine assign_field(lhs, f)
    !< Implementation of the (=) operator from another field_2d_t
    class(field_2d_t), intent(inout) :: lhs !< left-hand-side of the operation
    class(field_2d_t), intent(in) :: f

    if(enable_debug_print) then
      call debug_print('Running "'//trim(lhs%name)//'" field_2d_t%assign_field() ', __FILE__, __LINE__)
    endif

    call from_field(lhs, f)

    if(lhs%sync_on_equal) call lhs%sync_edges()
  endsubroutine assign_field

  subroutine assign_real_scalar(lhs, a)
    !< Implementation of the (=) operator from a real scalar
    class(field_2d_t), intent(inout) :: lhs !< left-hand-side of the operation
    real(rk), intent(in) :: a

    if(enable_debug_print) call debug_print('Running field_2d_t%assign_real_scalar() ', __FILE__, __LINE__)
    lhs%data = a
    if(lhs%sync_on_equal) call lhs%sync_edges()
  endsubroutine assign_real_scalar

  ! --------------------------------------------------------------------
  ! Arithmetic operators
  ! --------------------------------------------------------------------
  function field_add_field(lhs, f) result(res)
    !< Implementation of the field_2d_t + field_2d_t operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    class(field_2d_t), intent(in) :: f
    type(field_2d_t) :: res

    if(.not. field_shapes_match(lhs, f)) then
      error stop "Error in field_2d_t%field_add_field(): field_2d_t shapes don't match"
    endif

    if(lhs%use_gpu) then
      call field_add_field_gpu(lhs, f, res)
    else
      call from_field(res, lhs)
      call field_add_field_cpu(lhs, f, res)
    endif

  endfunction field_add_field

  function field_sub_field(lhs, f) result(res)
    !< Implementation of the field_2d_t + field_2d_t operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    class(field_2d_t), intent(in) :: f
    type(field_2d_t) :: res

    if(.not. field_shapes_match(lhs, f)) then
      error stop "Error in field_2d_t%field_add_field(): field_2d_t shapes don't match"
    endif

    if(lhs%use_gpu) then
      call field_sub_field_gpu(lhs, f, res)
    else
      call from_field(res, lhs)
      call field_sub_field_cpu(lhs, f, res)
    endif
  endfunction field_sub_field

  function field_mul_field(lhs, f) result(res)
    !< Implementation of the field_2d_t * field_2d_t operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    class(field_2d_t), intent(in) :: f
    type(field_2d_t) :: res

    if(.not. field_shapes_match(lhs, f)) then
      error stop "Error in field_2d_t%field_add_field(): field_2d_t shapes don't match"
    endif

    if(lhs%use_gpu) then
      call field_mul_field_gpu(lhs, f, res)
    else
      call from_field(res, lhs)
      call field_mul_field_cpu(lhs, f, res)
    endif

  endfunction field_mul_field

  function field_div_field(lhs, f) result(res)
    !< Implementation of the field_2d_t * field_2d_t operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    class(field_2d_t), intent(in) :: f
    type(field_2d_t) :: res

    if(.not. field_shapes_match(lhs, f)) then
      error stop "Error in field_2d_t%field_add_field(): field_2d_t shapes don't match"
    endif

    if(lhs%use_gpu) then
      call field_div_field_gpu(lhs, f, res)
    else
      call from_field(res, lhs)
      call field_div_field_cpu(lhs, f, res)
    endif

  endfunction field_div_field

  function field_add_real_1d(lhs, x) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), intent(in) :: x
    type(field_2d_t) :: res

    if(lhs%use_gpu) then
      call field_add_real_1d_gpu(lhs, x, res)
    else
      call from_field(res, lhs)
      call field_add_real_1d_cpu(lhs, x, res)
    endif

  endfunction field_add_real_1d

  function field_add_real_2d(lhs, x) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
    type(field_2d_t) :: res

    if(.not. shapes_match(lhs, x)) then
      error stop "Error in field_2d_t%field_add_real_2d(): field_2d_t and shapes don't match"
    endif

    if(lhs%use_gpu) then
      call field_add_real_2d_gpu(lhs, x, res)
    else
      call from_field(res, lhs)
      call field_add_real_2d_cpu(lhs, x, res)
    endif
  endfunction field_add_real_2d

  function real_1d_add_field(x, rhs) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), intent(in) :: x
    type(field_2d_t) :: res

    ! Reuse the field + x_1d version
    if(rhs%use_gpu) then
      call field_add_real_1d_gpu(rhs, x, res)
    else
      call from_field(res, rhs)
      call field_add_real_1d_cpu(rhs, x, res)
    endif

  endfunction real_1d_add_field

  function real_2d_add_field(x, rhs) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), dimension(rhs%ilo:, rhs%jlo:), intent(in) :: x
    type(field_2d_t) :: res

    if(.not. shapes_match(rhs, x)) then
      error stop "Error in field_2d_t%real_2d_add_field(): field_2d_t and x shapes don't match"
    endif

    ! Reuse the field + 2d version
    if(rhs%use_gpu) then
      call field_add_real_2d_gpu(rhs, x, res)
    else
      call from_field(res, rhs)
      call field_add_real_2d_cpu(rhs, x, res)
    endif
  endfunction real_2d_add_field

  function field_sub_real_1d(lhs, x) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), intent(in) :: x
    type(field_2d_t) :: res

    if(lhs%use_gpu) then
      call field_sub_real_1d_gpu(lhs, x, res)
    else
      call from_field(res, lhs)
      call field_sub_real_1d_cpu(lhs, x, res)
    endif

  endfunction field_sub_real_1d

  function field_sub_real_2d(lhs, x) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
    type(field_2d_t) :: res

    if(.not. shapes_match(lhs, x)) then
      error stop "Error in field_2d_t%field_sub_real_2d(): field_2d_t and x shapes don't match"
    endif

    if(lhs%use_gpu) then
      call field_sub_real_2d_gpu(lhs, x, res)
    else
      call from_field(res, lhs)
      call field_sub_real_2d_cpu(lhs, x, res)
    endif
  endfunction field_sub_real_2d

  function real_1d_sub_field(x, rhs) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), intent(in) :: x
    type(field_2d_t) :: res

    if(rhs%use_gpu) then
      call real_1d_sub_field_gpu(x, rhs, res)
    else
      call from_field(res, rhs)
      call real_1d_sub_field_cpu(x, rhs, res)
    endif

  endfunction real_1d_sub_field

  function real_2d_sub_field(x, rhs) result(res)
    !< Implementation of the field_2d_t + real64 operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), dimension(rhs%ilo:, rhs%jlo:), intent(in) :: x
    type(field_2d_t) :: res

    if(.not. shapes_match(rhs, x)) then
      error stop "Error in field_2d_t%real_2d_sub_field(): field_2d_t and x shapes don't match"
    endif

    if(rhs%use_gpu) then
      call real_2d_sub_field_gpu(x, rhs, res)
    else
      call from_field(res, rhs)
      call real_2d_sub_field_cpu(x, rhs, res)
    endif

  endfunction real_2d_sub_field

  function field_div_real_1d(lhs, x) result(res)
    !< Implementation of the field_2d_t / real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), intent(in) :: x
    type(field_2d_t) :: res

    real(rk) :: one_over_x

    one_over_x = 1.0_rk / x

    if(abs(x) < tiny(1.0_rk)) then
      error stop "Error in field_2d_t%field_div_real_1d: attempting a divide by 0 (abs(x) < tiny(1.0_rk))"
    endif

    if(lhs%use_gpu) then
      call field_mul_real_1d_gpu(lhs, one_over_x, res)
    else
      call from_field(res, lhs)
      call field_mul_real_1d_cpu(lhs, one_over_x, res)
    endif

  endfunction field_div_real_1d

  function field_div_real_2d(lhs, x) result(res)
    !< Implementation of the field_2d_t / real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
    type(field_2d_t) :: res

    if(.not. shapes_match(lhs, x)) then
      error stop "Error in field_2d_t%field_div_real_2d(): field_2d_t and x shapes don't match"
    endif

    if(lhs%use_gpu) then
      call field_div_real_2d_gpu(lhs, x, res)
    else
      call from_field(res, lhs)
      call field_div_real_2d_cpu(lhs, x, res)
    endif

  endfunction field_div_real_2d

  function real_1d_div_field(x, rhs) result(res)
    !< Implementation of the real64 / field_2d_t operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), intent(in) :: x
    type(field_2d_t) :: res

    if(rhs%use_gpu) then
      call real_1d_div_field_gpu(x, rhs, res)
    else
      call from_field(res, rhs)
      call real_1d_div_field_cpu(x, rhs, res)
    endif

  endfunction real_1d_div_field

  function real_2d_div_field(x, rhs) result(res)
    !< Implementation of the field_2d_t / real64 operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), dimension(rhs%ilo:, rhs%jlo:), intent(in) :: x
    type(field_2d_t) :: res

    if(.not. shapes_match(rhs, x)) then
      error stop "Error in field_2d_t%real_2d_div_field(): field_2d_t and x shapes don't match"
    endif

    if(rhs%use_gpu) then
      call real_2d_div_field_gpu(x, rhs, res)
    else
      call from_field(res, rhs)
      call real_2d_div_field_cpu(x, rhs, res)
    endif

  endfunction real_2d_div_field

  function field_mul_real_2d(lhs, x) result(res)
    !< Implementation of the field_2d_t * array operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), dimension(lhs%ilo:, lhs%jlo:), intent(in) :: x
    type(field_2d_t) :: res

    if(enable_debug_print) call debug_print('Running field_mul_real_2d_cpu ', __FILE__, __LINE__)
    if(.not. shapes_match(lhs, x)) then
      error stop "Error in field_2d_t%field_mul_real_2d(): field_2d_t and x shapes don't match"
    endif

    if(lhs%use_gpu) then
      call field_mul_real_2d_gpu(lhs, x, res)
    else
      call from_field(res, lhs)
      call field_mul_real_2d_cpu(lhs, x, res)
    endif
  endfunction field_mul_real_2d

  function field_mul_real_1d(lhs, x) result(res)
    !< Implementation of the field_2d_t * real64 operation
    class(field_2d_t), intent(in) :: lhs !< left-hand-side of the operation
    real(rk), intent(in) :: x
    type(field_2d_t) :: res

    call from_field(res, lhs)
    if(lhs%use_gpu) then
      call field_mul_real_1d_gpu(lhs, x, res)
    else
      call from_field(res, lhs)
      call field_mul_real_1d_cpu(lhs, x, res)
    endif

  endfunction field_mul_real_1d

  function real_1d_mul_field(x, rhs) result(res)
    !< Implementation of the real64 * field_2d_t operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), intent(in) :: x
    type(field_2d_t) :: res

    ! Reuse the field * 1d operator
    if(rhs%use_gpu) then
      call field_mul_real_1d_gpu(rhs, x, res)
    else
      call from_field(res, rhs)
      call field_mul_real_1d_cpu(rhs, x, res)
    endif

  endfunction real_1d_mul_field

  function real_2d_mul_field(x, rhs) result(res)
    !< Implementation of the real64 * field_2d_t operation
    class(field_2d_t), intent(in) :: rhs !< right-hand-side of the operation
    real(rk), dimension(rhs%ilo:, rhs%jlo:), intent(in) :: x
    type(field_2d_t) :: res

    if(.not. shapes_match(rhs, x)) then
      error stop "Error in field_2d_t%real_2d_mul_field(): field_2d_t and x shapes don't match"
    endif

    ! Reuse the field * 2d operator
    if(rhs%use_gpu) then
      call field_mul_real_2d_gpu(rhs, x, res)
    else
      call from_field(res, rhs)
      call field_mul_real_2d_cpu(rhs, x, res)
    endif

  endfunction real_2d_mul_field

  ! --------------------------------------------------------------------
  ! Reduction operators (min/max/sum val and loc)
  ! --------------------------------------------------------------------

  function field_maxval(self) result(res)
    class(field_2d_t), intent(in) :: self
    real(rk) :: res

    if(self%use_gpu) then
      res = field_maxval_gpu(self)
    else
      res = field_maxval_cpu(self)
    endif
  endfunction field_maxval

  function field_maxloc(self) result(res)
    class(field_2d_t), intent(in) :: self
    integer(ik), dimension(2) :: res

    if(self%use_gpu) then
      res = field_maxloc_gpu(self)
    else
      res = field_maxloc_cpu(self)
    endif
  endfunction field_maxloc

  function field_minval(self) result(res)
    class(field_2d_t), intent(in) :: self
    real(rk) :: res

    if(self%use_gpu) then
      res = field_minval_gpu(self)
    else
      res = field_minval_cpu(self)
    endif
  endfunction field_minval

  function field_minloc(self) result(res)
    class(field_2d_t), intent(in) :: self
    integer(ik), dimension(2) :: res

    if(self%use_gpu) then
      res = field_minloc_gpu(self)
    else
      res = field_minloc_cpu(self)
    endif
  endfunction field_minloc

  function field_sum(self) result(res)
    class(field_2d_t), intent(in) :: self
    real(rk) :: res

    if(self%use_gpu) then
      res = field_sum_gpu(self)
    else
      ! res = field_sum_cpu(self)
      res = sum(self%data(self%ilo:self%ihi, self%jlo:self%jhi))
    endif
  endfunction field_sum

endmodule mod_field
