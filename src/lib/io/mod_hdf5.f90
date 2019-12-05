! BSD 3-Clause License
!
! Copyright (c) 2018, Michael Hirsch, Ph.D.
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! * Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! * Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! * Neither the name of the copyright holder nor the names of its
! contributors may be used to endorse or promote products derived from
! this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! HDF5 object oriented interface
! Author: Michael Hirsch, https://www.scivision.dev/
! https://github.com/scivision/oo_hdf5_fortran

module hdf5_interface

  use, intrinsic :: iso_fortran_env, only: real32, real64, stderr => error_unit
  use H5LT

  implicit none

  public :: hdf5_file, toLower

  private

  type :: hdf5_file

    character(:), allocatable :: filename
    integer(HID_T) :: lid, &   !! location identifier
                      gid, &    !! group identifier
                      glid, &   !! group location identifier
                      sid, did, pid
    integer :: comp_lvl = 0 !! compression level (1-9)  0: disable compression
    integer(HSIZE_T) :: chunk_size(6) = [64, 64, 1, 1, 1, 1]  !! chunk size per dimension
    logical :: verbose = .false.

  contains
    !! initialize HDF5 file
    procedure, public :: initialize => hdf_initialize, &
      finalize => hdf_finalize, writeattr, &
      open => hdf_open_group, &
      close => hdf_close_group

    !! add group or dataset integer/real
    generic, public :: add => hdf_add_group, hdf_add_int, hdf_add_int1d, hdf_add_int2d, hdf_add_int3d, &
      hdf_add_real32, hdf_add_real32_1d, hdf_add_real32_2d, hdf_add_real32_3d, &
      hdf_add_real32_4d, hdf_add_real32_5d, hdf_add_real32_6d, &
      hdf_add_real64, hdf_add_real64_1d, hdf_add_real64_2d, hdf_add_real64_3d, &
      hdf_add_real64_4d, hdf_add_real64_5d, hdf_add_real64_6d, &
      hdf_add_string

    !! get dataset integer/real
    generic, public :: get => hdf_get_int, hdf_get_int1d, hdf_get_int2d, hdf_get_int3d, &
      hdf_get_real32, hdf_get_real32_1d, hdf_get_real32_2d, hdf_get_real32_3d, &
      hdf_get_real64, hdf_get_real64_1d, hdf_get_real64_2d, hdf_get_real64_3d, hdf_get_real64_4d, &
      hdf_get_string

    !! private methods
    procedure, private :: hdf_add_group, &
      hdf_add_int, hdf_add_int1d, hdf_add_int2d, hdf_add_int3d, &
      hdf_get_int, hdf_get_int1d, hdf_get_int2d, hdf_get_int3d, &
      hdf_add_real32, hdf_add_real32_1d, hdf_add_real32_2d, hdf_add_real32_3d, &
      hdf_add_real32_4d, hdf_add_real32_5d, hdf_add_real32_6d, &
      hdf_add_real64, hdf_add_real64_1d, hdf_add_real64_2d, hdf_add_real64_3d, &
      hdf_add_real64_4d, hdf_add_real64_5d, hdf_add_real64_6d, &
      hdf_get_real32, hdf_get_real32_1d, hdf_get_real32_2d, hdf_get_real32_3d, &
      hdf_get_real64, hdf_get_real64_1d, hdf_get_real64_2d, hdf_get_real64_3d, hdf_get_real64_4d, &
      hdf_add_string, hdf_get_string

  end type hdf5_file

contains
  !=============================================================================
  subroutine hdf_initialize(self, filename, status, action, comp_lvl)
    !! Opens hdf5 file

    class(hdf5_file), intent(inout) :: self
    character(*), intent(in) :: filename
    character(*), intent(in), optional :: status
    character(*), intent(in), optional :: action
    integer, intent(in), optional :: comp_lvl

    character(:), allocatable :: lstatus, laction
    integer :: ierr

    self%filename = filename

    if(present(comp_lvl)) self%comp_lvl = comp_lvl

    !! Initialize FORTRAN interface.
    call h5open_f(ierr)
    if(ierr /= 0) error stop 'Error: HDF5 library initialize Failed!'

    lstatus = 'old'  ! not merge() due to unequal character length
    if(present(status)) lstatus = toLower(status)

    laction = 'rw'  ! not merge() due to unequal character length
    if(present(action)) laction = toLower(action)

    select case(lstatus)
    case('old')
      select case(laction)
      case('read', 'r')  !! Open an existing file.
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, self%lid, ierr)
      case('write', 'readwrite', 'w', 'rw')
        call h5fopen_f(filename, H5F_ACC_RDWR_F, self%lid, ierr)
      case default
        print *, 'Error: Unsupported action ->'//laction
        error stop
      endselect
    case('new', 'replace')
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, self%lid, ierr)
    case default
      print *, 'Error: Unsupported status ->'//lstatus
      error stop
    endselect

    if(ierr /= 0) then
      print *, 'Error: HDF5 open/create failed: '//filename
      error stop
    endif
  end subroutine hdf_initialize
  !=============================================================================
  subroutine hdf_finalize(self)
    class(hdf5_file), intent(in) :: self

    integer :: ierr

    !! close hdf5 file
    call h5fclose_f(self%lid, ierr)

    !!  Close FORTRAN interface.
    call h5close_f(ierr)

    if(ierr /= 0) then
      print *, 'Error: HDF5 finalization: '//self%filename
      error stop
    end if

  end subroutine hdf_finalize
  !=============================================================================
  subroutine hdf_add_group(self, gname)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: gname    !! relative path to group

    integer(HID_T) :: gid

    integer :: ierr, sp, ep, sl
    logical :: gexist

    sl = len(gname)
    sp = 1
    ep = 0

    do
      ep = index(gname(sp + 1:sl), "/")

      ! no subgroup found
      if(ep == 0) exit

      ! check subgroup exists
      sp = sp + ep
      call h5lexists_f(self%lid, gname(1:sp - 1), gexist, ierr)
      if(ierr /= 0) then
        print *, 'problem finding group '//gname
        error stop
      endif

      if(.not. gexist) then
        call h5gcreate_f(self%lid, gname(1:sp - 1), gid, ierr)
        if(ierr /= 0) then
          print *, 'problem creating group '//gname
          error stop
        end if

        call h5gclose_f(gid, ierr)
        if(ierr /= 0) then
          print *, 'problem closing group '//gname
          error stop
        end if
      endif
    end do

  end subroutine hdf_add_group
  !=============================================================================
  subroutine hdf_open_group(self, gname)
    class(hdf5_file), intent(inout) :: self
    character(*), intent(in) :: gname

    integer :: ierr

    call h5gopen_f(self%lid, gname, self%gid, ierr)
    if(ierr /= 0) then
      print *, 'problem opening group '//gname
      error stop
    end if

    self%glid = self%lid
    self%lid = self%gid

  end subroutine hdf_open_group
  !=============================================================================
  subroutine hdf_close_group(self)
    class(hdf5_file), intent(inout) :: self

    integer :: ierr

    call h5gclose_f(self%gid, ierr)
    if(ierr /= 0) then
      print *, 'problem closing group '//self%filename
      error stop
    end if

    self%lid = self%glid

  end subroutine hdf_close_group
  !=============================================================================
  subroutine hdf_set_deflate(self, dims)
    class(hdf5_file), intent(inout) :: self
    integer(HSIZE_T), intent(in) :: dims(:)

    integer :: ierr, ndims, i
    integer(HSIZE_T), allocatable :: chunk_size(:)

    ndims = size(dims)
    allocate(chunk_size(ndims))

    do concurrent(i=1:ndims)
      chunk_size(i) = min(self%chunk_size(i), dims(i))
    enddo

    if(self%verbose) print *, 'dims: ', dims, 'chunk size: ', chunk_size

    call h5pcreate_f(H5P_DATASET_CREATE_F, self%pid, ierr)
    if(ierr /= 0) then
      print *, 'error creating property '//self%filename
      error stop
    end if

    call h5pset_chunk_f(self%pid, ndims, chunk_size, ierr)
    if(ierr /= 0) then
      print *, 'error setting chunk '//self%filename
      error stop
    end if

    if(self%comp_lvl < 1 .or. self%comp_lvl > 9) return

    call h5pset_shuffle_f(self%pid, ierr)
    if(ierr /= 0) then
      print *, 'error enabling Shuffle '//self%filename
      error stop
    end if

    call h5pset_deflate_f(self%pid, self%comp_lvl, ierr)
    if(ierr /= 0) then
      print *, 'error enabling Deflate compression '//self%filename
      error stop
    end if

  end subroutine hdf_set_deflate

  subroutine hdf_setup_write(self, dname, dtype, dims, chunk_size)
    class(hdf5_file), intent(inout) :: self
    character(*), intent(in) :: dname
    integer(HID_T), intent(in) :: dtype
    integer(HSIZE_T), intent(in) :: dims(:)
    integer, intent(in), optional :: chunk_size(:)

    integer :: ierr

    call self%add(dname)

    if(present(chunk_size)) self%chunk_size(:size(dims)) = chunk_size

    call hdf_set_deflate(self, dims)

    call h5screate_simple_f(size(dims), dims, self%sid, ierr)
    if(ierr /= 0) then
      print *, 'error on dataspace '//dname//' '//self%filename
      error stop
    end if

    call h5dcreate_f(self%lid, dname, dtype, self%sid, self%did, ierr, self%pid)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' '//self%filename
      error stop
    end if

  end subroutine hdf_setup_write

  subroutine hdf_wrapup(self)
    class(hdf5_file), intent(in) :: self
    integer :: ierr

    call h5sclose_f(self%sid, ierr)
    call h5pclose_f(self%pid, ierr)
    call h5dclose_f(self%did, ierr)
    if(ierr /= 0) then
      print *, 'error on closing dataset '//self%filename
      error stop
    end if

  end subroutine hdf_wrapup

  subroutine writeattr(self, dname, attr, attrval)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname, attr, attrval

    integer :: ierr
    logical :: exists

    call self%add(dname)

    call h5ltpath_valid_f(self%lid, dname, .true., exists, ierr)
    if(ierr /= 0) then
      print *, 'problem checking existence: '//dname//' file '//self%filename
      error stop
    end if

    if(.not. exists) then
      write(stderr, *) 'WARNING: variable '//dname//' must be created before writing '//attr
      return
    endif

    call h5ltset_attribute_string_f(self%lid, dname, attr, attrval, ierr)
    if(ierr /= 0) then
      print *, 'problem writing attribute '//attr//' to '//dname//' file '//self%filename
      error stop
    end if

  end subroutine writeattr
  !===================================
  subroutine hdf_add_int(self, dname, value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    integer, intent(in) :: value

    integer(HID_T) :: sid, did
    integer :: ierr

    call self%add(dname)

    ! HDF5 >= 1.10
    !call h5ltmake_dataset_f(self%lid, dname, &
    !  rank(value), int(shape(value),HSIZE_T), h5kind_to_type(kind(value),H5_INTEGER_KIND), value, ierr)
    !if (ierr /= 0) then
    !  print*, 'error on dataset '//dname//' write '//self%filename
    !  stop
    !end if

    !  HDF5 1.8 compatbility below:
    !! create dataspace
    call h5screate_f(H5S_SCALAR_F, sid, ierr)
    if(ierr /= 0) then
      print *, 'error create dataspace '//dname//' write '//self%filename
      stop
    end if

    !! create dataset
    call h5dcreate_f(self%lid, dname, h5kind_to_type(kind(value), H5_INTEGER_KIND), sid, did, ierr)
    if(ierr /= 0) then
      print *, 'error create dataspace '//dname//' write '//self%filename
      stop
    end if

    !! write dataset
    call h5dwrite_f(did, h5kind_to_type(kind(value), H5_INTEGER_KIND), value, int(shape(value), HSIZE_T), ierr)
    if(ierr /= 0) then
      print *, 'error write dataspace '//dname//' write '//self%filename
      error stop
    end if

    !! close space and dataset
    call h5dclose_f(did, ierr)
    if(ierr /= 0) then
      print *, 'error close dataset '//dname//' write '//self%filename
      error stop
    end if
    call h5sclose_f(sid, ierr)
    if(ierr /= 0) then
      print *, 'error close dataspace '//dname//' write '//self%filename
      error stop
    end if

  end subroutine hdf_add_int

  subroutine hdf_add_int1d(self, dname, value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    integer, intent(in) :: value(:)

    integer :: ierr

    call self%add(dname)

    call h5ltmake_dataset_f(self%lid, dname, &
                            rank(value), int(shape(value), HSIZE_T), h5kind_to_type(kind(value), H5_INTEGER_KIND), value, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

  end subroutine hdf_add_int1d

  subroutine hdf_add_int2d(self, dname, value, chunk_size)
    class(hdf5_file), intent(inout) :: self
    character(*), intent(in) :: dname
    integer, intent(in) :: value(:, :)
    integer, intent(in), optional :: chunk_size(:)

    integer :: ierr
    integer(HID_T) :: dtype
    integer(HSIZE_T) :: dims(rank(value))

    dims = shape(value)
    dtype = h5kind_to_type(kind(value), H5_INTEGER_KIND)

    call hdf_setup_write(self, dname, dtype, dims, chunk_size)

    call h5dwrite_f(self%did, dtype, value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

    call hdf_wrapup(self)

  end subroutine hdf_add_int2d

  subroutine hdf_add_int3d(self, dname, value, chunk_size)
    class(hdf5_file), intent(inout) :: self
    character(*), intent(in) :: dname
    integer, intent(in) :: value(:, :, :)
    integer, intent(in), optional :: chunk_size(:)

    integer :: ierr
    integer(HID_T) :: dtype
    integer(HSIZE_T) :: dims(rank(value))

    dims = shape(value)
    dtype = h5kind_to_type(kind(value), H5_INTEGER_KIND)

    call hdf_setup_write(self, dname, dtype, dims, chunk_size)

    call h5dwrite_f(self%did, dtype, value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

    call hdf_wrapup(self)

  end subroutine hdf_add_int3d

  !========================================================================================

  subroutine hdf_add_real32(self, dname, value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(real32), intent(in) :: value

    integer(HID_T) :: sid, did
    integer :: ierr

    call self%add(dname)

    ! HDF5 >= 1.10
    !call h5ltmake_dataset_f(self%lid, dname, &
    !   rank(value), int(shape(value),HSIZE_T), h5kind_to_type(kind(value),H5_REAL_KIND), value, ierr)
    !if (ierr /= 0) error print*, 'error on dataset '//dname//' write '//self%filename

    ! HDF5 1.8 compatbility below:
    !! create dataspace
    call h5screate_f(H5S_SCALAR_F, sid, ierr)
    if(ierr /= 0) then
      print *, 'error create dataspace '//dname//' write '//self%filename
      error stop
    end if

    !! create dataset
    call h5dcreate_f(self%lid, dname, h5kind_to_type(kind(value), H5_REAL_KIND), sid, did, ierr)
    if(ierr /= 0) then
      print *, 'error create dataset '//dname//' write '//self%filename
      error stop
    end if

    !! write dataset
    call h5dwrite_f(did, h5kind_to_type(kind(value), H5_REAL_KIND), value, int(shape(value), HSIZE_T), ierr)
    if(ierr /= 0) then
      print *, 'error write dataset '//dname//' write '//self%filename
      error stop
    end if

    !! close space and dataset
    call h5dclose_f(did, ierr)
    call h5sclose_f(sid, ierr)
    if(ierr /= 0) then
      print *, 'error close dataspace '//dname//' write '//self%filename
      error stop
    end if

  end subroutine hdf_add_real32

  subroutine hdf_add_real32_1d(self, dname, value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(real32), intent(in) :: value(:)

    integer :: ierr

    call self%add(dname)

    call h5ltmake_dataset_f(self%lid, dname, &
                            rank(value), int(shape(value), HSIZE_T), h5kind_to_type(kind(value), H5_REAL_KIND), value, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

  end subroutine hdf_add_real32_1d

  subroutine hdf_add_real32_2d(self, dname, value, chunk_size)
    class(hdf5_file), intent(inout) :: self
    character(*), intent(in) :: dname
    real(real32), intent(in) :: value(:, :)
    integer, intent(in), optional :: chunk_size(:)

    integer :: ierr
    integer(HID_T) :: dtype
    integer(HSIZE_T) :: dims(rank(value))

    dims = shape(value)
    dtype = h5kind_to_type(kind(value), H5_REAL_KIND)

    call hdf_setup_write(self, dname, dtype, dims, chunk_size)

    call h5dwrite_f(self%did, dtype, value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

    call hdf_wrapup(self)

  end subroutine hdf_add_real32_2d

  subroutine hdf_add_real32_3d(self, dname, value, chunk_size)
    class(hdf5_file), intent(inout) :: self
    character(*), intent(in) :: dname
    real(real32), intent(in) :: value(:, :, :)
    integer, intent(in), optional :: chunk_size(:)

    integer :: ierr
    integer(HID_T) :: dtype
    integer(HSIZE_T) :: dims(rank(value))

    dims = shape(value)
    dtype = h5kind_to_type(kind(value), H5_REAL_KIND)

    call hdf_setup_write(self, dname, dtype, dims, chunk_size)

    call h5dwrite_f(self%did, dtype, value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

    call hdf_wrapup(self)

  end subroutine hdf_add_real32_3d

  subroutine hdf_add_real32_4d(self, dname, value, chunk_size)
    class(hdf5_file), intent(inout) :: self
    character(*), intent(in) :: dname
    real(real32), intent(in) :: value(:, :, :, :)
    integer, intent(in), optional :: chunk_size(:)

    integer :: ierr
    integer(HID_T) :: dtype
    integer(HSIZE_T) :: dims(rank(value))

    dims = shape(value)
    dtype = h5kind_to_type(kind(value), H5_REAL_KIND)

    call hdf_setup_write(self, dname, dtype, dims, chunk_size)

    call h5dwrite_f(self%did, dtype, value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

    call hdf_wrapup(self)

  end subroutine hdf_add_real32_4d

  subroutine hdf_add_real32_5d(self, dname, value, chunk_size)
    class(hdf5_file), intent(inout) :: self
    character(*), intent(in) :: dname
    real(real32), intent(in) :: value(:, :, :, :, :)
    integer, intent(in), optional :: chunk_size(:)

    integer :: ierr
    integer(HID_T) :: dtype
    integer(HSIZE_T) :: dims(rank(value))

    dims = shape(value)
    dtype = h5kind_to_type(kind(value), H5_REAL_KIND)

    call hdf_setup_write(self, dname, dtype, dims, chunk_size)

    call h5dwrite_f(self%did, dtype, value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

    call hdf_wrapup(self)

  end subroutine hdf_add_real32_5d

  subroutine hdf_add_real32_6d(self, dname, value, chunk_size)
    class(hdf5_file), intent(inout) :: self
    character(*), intent(in) :: dname
    real(real32), intent(in) :: value(:, :, :, :, :, :)
    integer, intent(in), optional :: chunk_size(:)

    integer :: ierr
    integer(HID_T) :: dtype
    integer(HSIZE_T) :: dims(rank(value))

    dims = shape(value)
    dtype = h5kind_to_type(kind(value), H5_REAL_KIND)

    call hdf_setup_write(self, dname, dtype, dims, chunk_size)

    call h5dwrite_f(self%did, dtype, value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

    call hdf_wrapup(self)

  end subroutine hdf_add_real32_6d

  subroutine hdf_add_real64(self, dname, value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(real64), intent(in) :: value

    integer(HID_T) :: sid, did
    integer :: ierr

    call self%add(dname)

    ! HDF5 >= 1.10
    !call h5ltmake_dataset_f(self%lid, dname, &
    !   rank(value), int(shape(value),HSIZE_T), h5kind_to_type(kind(value),H5_REAL_KIND), value, ierr)
    !if (ierr /= 0) then
    !      print*, 'error on dataset '//dname//' write '//self%filename
    !      stop
    !    end if

    ! HDF5 1.8 compatbility below:
    !! create dataspace
    call h5screate_f(H5S_SCALAR_F, sid, ierr)
    if(ierr /= 0) then
      print *, 'error create dataspace '//dname//' write '//self%filename
      error stop
    end if

    !! create dataset
    call h5dcreate_f(self%lid, dname, h5kind_to_type(kind(value), H5_REAL_KIND), sid, did, ierr)
    if(ierr /= 0) then
      print *, 'error create dataset '//dname//' write '//self%filename
      error stop
    end if

    !! write dataset
    call h5dwrite_f(did, h5kind_to_type(kind(value), H5_REAL_KIND), value, int(shape(value), HSIZE_T), ierr)
    if(ierr /= 0) then
      print *, 'error write dataset '//dname//' write '//self%filename
      error stop
    end if

    !! close space and dataset
    call h5dclose_f(did, ierr)
    call h5sclose_f(sid, ierr)
    if(ierr /= 0) then
      print *, 'error close dataspace '//dname//' write '//self%filename
      error stop
    end if

  end subroutine hdf_add_real64

  subroutine hdf_add_real64_1d(self, dname, value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(real64), intent(in) :: value(:)

    integer :: ierr

    call self%add(dname)

    call h5ltmake_dataset_f(self%lid, dname, &
                            rank(value), int(shape(value), HSIZE_T), h5kind_to_type(kind(value), H5_REAL_KIND), value, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

  end subroutine hdf_add_real64_1d

  subroutine hdf_add_real64_2d(self, dname, value, chunk_size)
    class(hdf5_file), intent(inout) :: self
    character(*), intent(in) :: dname
    real(real64), intent(in) :: value(:, :)
    integer, intent(in), optional :: chunk_size(:)

    integer :: ierr
    integer(HID_T) :: dtype
    integer(HSIZE_T) :: dims(rank(value))

    dims = shape(value)
    dtype = h5kind_to_type(kind(value), H5_REAL_KIND)

    call hdf_setup_write(self, dname, dtype, dims, chunk_size)

    call h5dwrite_f(self%did, dtype, value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

    call hdf_wrapup(self)

  end subroutine hdf_add_real64_2d

  subroutine hdf_add_real64_3d(self, dname, value, chunk_size)
    class(hdf5_file), intent(inout) :: self
    character(*), intent(in) :: dname
    real(real64), intent(in) :: value(:, :, :)
    integer, intent(in), optional :: chunk_size(:)

    integer :: ierr
    integer(HID_T) :: dtype
    integer(HSIZE_T) :: dims(rank(value))

    dims = shape(value)
    dtype = h5kind_to_type(kind(value), H5_REAL_KIND)

    call hdf_setup_write(self, dname, dtype, dims, chunk_size)

    call h5dwrite_f(self%did, dtype, value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

    call hdf_wrapup(self)

  end subroutine hdf_add_real64_3d

  subroutine hdf_add_real64_4d(self, dname, value, chunk_size)
    class(hdf5_file), intent(inout) :: self
    character(*), intent(in) :: dname
    real(real64), intent(in) :: value(:, :, :, :)
    integer, intent(in), optional :: chunk_size(:)

    integer :: ierr
    integer(HID_T) :: dtype
    integer(HSIZE_T) :: dims(rank(value))

    dims = shape(value)
    dtype = h5kind_to_type(kind(value), H5_REAL_KIND)

    call hdf_setup_write(self, dname, dtype, dims, chunk_size)

    call h5dwrite_f(self%did, dtype, value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

    call hdf_wrapup(self)

  end subroutine hdf_add_real64_4d

  subroutine hdf_add_real64_5d(self, dname, value, chunk_size)
    class(hdf5_file), intent(inout) :: self
    character(*), intent(in) :: dname
    real(real64), intent(in) :: value(:, :, :, :, :)
    integer, intent(in), optional :: chunk_size(:)

    integer :: ierr
    integer(HID_T) :: dtype
    integer(HSIZE_T) :: dims(rank(value))

    dims = shape(value)
    dtype = h5kind_to_type(kind(value), H5_REAL_KIND)

    call hdf_setup_write(self, dname, dtype, dims, chunk_size)

    call h5dwrite_f(self%did, dtype, value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

    call hdf_wrapup(self)

  end subroutine hdf_add_real64_5d

  subroutine hdf_add_real64_6d(self, dname, value, chunk_size)
    class(hdf5_file), intent(inout) :: self
    character(*), intent(in) :: dname
    real(real64), intent(in) :: value(:, :, :, :, :, :)
    integer, intent(in), optional :: chunk_size(:)

    integer :: ierr
    integer(HID_T) :: dtype
    integer(HSIZE_T) :: dims(rank(value))

    dims = shape(value)
    dtype = h5kind_to_type(kind(value), H5_REAL_KIND)

    call hdf_setup_write(self, dname, dtype, dims, chunk_size)

    call h5dwrite_f(self%did, dtype, value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

    call hdf_wrapup(self)

  end subroutine hdf_add_real64_6d

  !================================================================================

  subroutine hdf_add_string(self, dname, value)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname, value
    integer :: ierr

    call h5ltmake_dataset_string_f(self%lid, dname, value, ierr)
    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' write '//self%filename
      stop
    end if

  end subroutine hdf_add_string

  !====== READ =====================================

  subroutine hdf_get_string(self, dname, value)
    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    character(*), intent(out) :: value

    integer :: ierr, i

    call h5ltread_dataset_string_f(self%lid, dname, value, ierr)

    ! gets rid of the ^@ character (null)
    do i = 1, len(value)
      if(value(i:i) == char(0)) value(i:i) = ' '
    end do

    if(ierr /= 0) then
      print *, 'error on dataset '//dname//' read '//self%filename
      error stop
    end if

  end subroutine hdf_get_string

  subroutine hdf_get_int(self, dname, value)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    integer, intent(out) :: value

    integer(HID_T) :: did
    integer :: ierr

    ! open dataset
    call h5dopen_f(self%lid, dname, did, ierr)
    if(ierr /= 0) then
      print *, 'error open dataset '//dname//' read '//self%filename
      error stop
    end if

    ! read dataset
    call h5dread_f(did, h5kind_to_type(kind(value), H5_INTEGER_KIND), value, int(shape(value), HSIZE_T), ierr)
    if(ierr /= 0) then
      print *, 'error read dataset '//dname//' read '//self%filename
      error stop
    end if

    ! close dataset
    call h5dclose_f(did, ierr)
    if(ierr /= 0) then
      print *, 'error close dataset '//dname//' read '//self%filename
      error stop
    end if

  end subroutine hdf_get_int

  subroutine hdf_get_int1d(self, dname, value)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    integer, intent(out), allocatable :: value(:)

    integer(SIZE_T) :: dims(1), dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
    if(ierr /= 0) then
      print *, 'error open dataset '//dname//' read '//self%filename
      error stop
    end if

    allocate(value(dims(1)))

    call h5ltread_dataset_f(self%lid, dname, h5kind_to_type(kind(value), H5_INTEGER_KIND), value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error read dataset '//dname//' read '//self%filename
      error stop
    end if

  end subroutine hdf_get_int1d

  subroutine hdf_get_int2d(self, dname, value)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    integer, intent(out), allocatable :: value(:, :)

    integer(SIZE_T) :: dims(2), dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
    if(ierr /= 0) then
      print *, 'error open dataset '//dname//' read '//self%filename
      error stop
    end if

    allocate(value(dims(1), dims(2)))

    call h5ltread_dataset_f(self%lid, dname, h5kind_to_type(kind(value), H5_INTEGER_KIND), value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error read dataset '//dname//' read '//self%filename
      error stop
    end if

  end subroutine hdf_get_int2d

  subroutine hdf_get_int3d(self, dname, value)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    integer, intent(out), allocatable :: value(:, :, :)

    integer(SIZE_T) :: dims(3), dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
    if(ierr /= 0) then
      print *, 'error open dataset '//dname//' read '//self%filename
      error stop
    end if

    allocate(value(dims(1), dims(2), dims(3)))

    call h5ltread_dataset_f(self%lid, dname, h5kind_to_type(kind(value), H5_INTEGER_KIND), value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error read dataset '//dname//' read '//self%filename
      error stop
    end if

  end subroutine hdf_get_int3d
  !=============================================================================
  subroutine hdf_get_real32(self, dname, value)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(real32), intent(out) :: value

    integer(HID_T) :: set_id
    integer :: ierr

    ! open dataset
    call h5dopen_f(self%lid, dname, set_id, ierr)

    ! read dataset
    call h5dread_f(set_id, h5kind_to_type(kind(value), H5_REAL_KIND), value, int(shape(value), HSIZE_T), ierr)

    ! close dataset
    call h5dclose_f(set_id, ierr)

  end subroutine hdf_get_real32

  subroutine hdf_get_real32_1d(self, dname, value)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(real32), intent(out), allocatable :: value(:)

    integer(SIZE_T) :: dims(1), dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
    if(ierr /= 0) then
      print *, 'error open dataset '//dname//' read '//self%filename
      error stop
    end if

    allocate(value(dims(1)))

    call h5ltread_dataset_f(self%lid, dname, h5kind_to_type(kind(value), H5_REAL_KIND), value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error read dataset '//dname//' read '//self%filename
      error stop
    end if

  end subroutine hdf_get_real32_1d

  subroutine hdf_get_real32_2d(self, dname, value)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(real32), intent(out), allocatable :: value(:, :)

    integer(SIZE_T) :: dims(2), dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
    if(ierr /= 0) then
      print *, 'error open dataset '//dname//' read '//self%filename
      error stop
    end if

    allocate(value(dims(1), dims(2)))

    call h5ltread_dataset_f(self%lid, dname, h5kind_to_type(kind(value), H5_REAL_KIND), value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error read dataset '//dname//' read '//self%filename
      error stop
    end if

  end subroutine hdf_get_real32_2d

  subroutine hdf_get_real32_3d(self, dname, value)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(real32), intent(out), allocatable :: value(:, :, :)

    integer(SIZE_T) :: dims(3), dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
    if(ierr /= 0) then
      print *, 'error open dataset '//dname//' read '//self%filename
      error stop
    end if

    allocate(value(dims(1), dims(2), dims(3)))

    call h5ltread_dataset_f(self%lid, dname, h5kind_to_type(kind(value), H5_REAL_KIND), value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error read dataset '//dname//' read '//self%filename
      error stop
    end if

  end subroutine hdf_get_real32_3d

  subroutine hdf_get_real64(self, dname, value)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(real64), intent(out) :: value

    integer(HID_T) :: set_id
    integer :: ierr

    ! open dataset
    call h5dopen_f(self%lid, dname, set_id, ierr)

    ! read dataset
    call h5dread_f(set_id, h5kind_to_type(kind(value), H5_REAL_KIND), value, int(shape(value), HSIZE_T), ierr)

    ! close dataset
    call h5dclose_f(set_id, ierr)

  end subroutine hdf_get_real64

  subroutine hdf_get_real64_1d(self, dname, value)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(real64), intent(out), allocatable :: value(:)

    integer(SIZE_T) :: dims(1), dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
    if(ierr /= 0) then
      print *, 'error open dataset '//dname//' read '//self%filename
      error stop
    end if

    allocate(value(dims(1)))

    call h5ltread_dataset_f(self%lid, dname, h5kind_to_type(kind(value), H5_REAL_KIND), value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error read dataset '//dname//' read '//self%filename
      error stop
    end if

  end subroutine hdf_get_real64_1d

  subroutine hdf_get_real64_2d(self, dname, value)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(real64), intent(out), allocatable :: value(:, :)

    integer(SIZE_T) :: dims(2), dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
    if(ierr /= 0) then
      print *, 'error open dataset '//dname//' read '//self%filename
      error stop
    end if

    allocate(value(dims(1), dims(2)))

    call h5ltread_dataset_f(self%lid, dname, h5kind_to_type(kind(value), H5_REAL_KIND), value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error read dataset '//dname//' read '//self%filename
      error stop
    end if

  end subroutine hdf_get_real64_2d

  subroutine hdf_get_real64_3d(self, dname, value)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(real64), intent(out), allocatable :: value(:, :, :)

    integer(SIZE_T) :: dims(3), dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
    if(ierr /= 0) then
      print *, 'error open dataset '//dname//' read '//self%filename
      error stop
    end if

    allocate(value(dims(1), dims(2), dims(3)))

    call h5ltread_dataset_f(self%lid, dname, h5kind_to_type(kind(value), H5_REAL_KIND), value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error read dataset '//dname//' read '//self%filename
      error stop
    end if

  end subroutine hdf_get_real64_3d

  subroutine hdf_get_real64_4d(self, dname, value)

    class(hdf5_file), intent(in) :: self
    character(*), intent(in) :: dname
    real(real64), intent(out), allocatable :: value(:, :, :, :)

    integer(SIZE_T) :: dims(4), dsize
    integer :: ierr, dtype

    call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ierr)
    if(ierr /= 0) then
      print *, 'error open dataset '//dname//' read '//self%filename
      error stop
    end if

    allocate(value(dims(1), dims(2), dims(3), dims(4)))

    call h5ltread_dataset_f(self%lid, dname, h5kind_to_type(kind(value), H5_REAL_KIND), value, dims, ierr)
    if(ierr /= 0) then
      print *, 'error read dataset '//dname//' read '//self%filename
      error stop
    end if

  end subroutine hdf_get_real64_4d

  !----- Helper functions

  elemental function toLower(str)
    ! can be trivially extended to non-ASCII
    character(*), intent(in) :: str
    character(len(str)) :: toLower
    character(*), parameter :: lower = "abcdefghijklmnopqrstuvwxyz", &
                               upper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    integer :: i, j

    toLower = str

    do concurrent(i=1:len(str))
      j = index(upper, str(i:i))
      if(j > 0) toLower(i:i) = lower(j:j)
    end do

  end function toLower

end module hdf5_interface
