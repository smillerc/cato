module mod_xdmf

  use, intrinsic :: iso_fortran_env, only: int32, real32, real64, stderr => error_unit
  use mod_contour_writer, only: base_io_t
  use hdf5_interface, only: hdf5_file

  implicit none

  type, extends(base_io_t) :: xdmf_writer_t
    private
    type(hdf5_file) :: hdf5_file
    character(len=:), allocatable :: format !< xdmf or just plain hdf5
    character(len=:), allocatable :: hdf5_filename
    character(len=:), allocatable :: xdmf_filename

  contains
    private
    ! Public methods
    procedure, public :: initialize
    procedure, public :: finalize

    ! Private methods
    procedure :: add_field_int_2d
    procedure :: add_field_real32_2d
    procedure :: add_field_real64_2d

    ! Generic interface
    generic, public :: add => add_field_int_2d, add_field_real32_2d, add_field_real64_2d

    ! Internal
    final :: finalize_xdmf
  end type xdmf_writer_t

contains

  subroutine add_field_int_2d(self, name, data)
    class(xdmf_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: name
    integer(int32), dimension(:, :), intent(in) :: data
  end subroutine add_field_int_2d

  subroutine add_field_real32_2d(self, name, data)
    class(xdmf_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(real32), dimension(:, :), intent(in) :: data
  end subroutine add_field_real32_2d

  subroutine add_field_real64_2d(self, name, data)
    class(xdmf_writer_t), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(real64), dimension(:, :), intent(in) :: data
  end subroutine add_field_real64_2d

end module mod_xdmf
