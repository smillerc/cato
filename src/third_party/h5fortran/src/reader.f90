! BSD 3-Clause License

! Copyright (c) 2018-2020, Michael Hirsch, Ph.D.
! All rights reserved.

! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:

! * Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer.

! * Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.

! * Neither the name of the copyright holder nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.

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

submodule(h5fortran:read) reader
!! This submodule is for reading 0-D..7-D data

use hdf5, only: h5dread_f
use h5lt, only: h5ltread_dataset_string_f

implicit none(type, external)

contains

module procedure hdf_read_scalar

integer(HSIZE_T) :: dims(rank(value))
integer(hid_t) :: did, sid
integer :: ier

if(.not. self%is_open) error stop 'h5fortran:reader: file handle is not open'

sid = 0

if(.not. self%exist(dname)) then
  write(stderr, *) 'h5fortran:ERROR: '//dname//' does not exist in '//self%filename
  error stop
endif

call h5dopen_f(self%lid, dname, did, ier)
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR: '//dname//' could not be opened in '//self%filename
  error stop
endif

select type(value)
type is(character(*))
  call hdf_wrapup(did, sid, ier)  !< FIXME: till character is treated same as other types
  block
    character(len(value)) :: buf
    call h5ltread_dataset_string_f(self%lid, dname, buf, ier)
    value = buf
  end block
  return
type is(real(real64))
  call h5dread_f(did, H5T_NATIVE_DOUBLE, value, dims, ier)
type is(real(real32))
  call h5dread_f(did, H5T_NATIVE_REAL, value, dims, ier)
type is(integer(int32))
  call h5dread_f(did, H5T_NATIVE_INTEGER, value, dims, ier)
class default
  error stop 'h5fortran:reader: incorrect data type'
end select

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_read_scalar

module procedure hdf_read_1d
integer(HSIZE_T) :: dims(rank(value))
integer(HID_T) :: did, sid, mem_sid
integer :: ier

did = 0 !< sentinel
sid = H5S_ALL_F
mem_sid = H5S_ALL_F
dims = shape(value)

if(present(istart) .and. present(iend)) then
  if(present(stride)) then
    !! necessary to use this present check for Intel and GCC
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend, stride)
  else
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend)
  endif
else
  call hdf_shape_check(self, dname, dims)
  call h5dopen_f(self%lid, dname, did, ier)
  if(ier /= 0) then
    write(stderr, *) 'h5fortran:ERROR:reader: could not setup read ', dname, ' from ', self%filename
    error stop
  endif
endif

select type(value)
type is(real(real64))
  call h5dread_f(did, H5T_NATIVE_DOUBLE, value, dims, ier, mem_sid, sid)
type is(real(real32))
  call h5dread_f(did, H5T_NATIVE_REAL, value, dims, ier, mem_sid, sid)
type is(integer(int32))
  call h5dread_f(did, H5T_NATIVE_INTEGER, value, dims, ier, mem_sid, sid)
class default
  error stop 'h5fortran:reader: incorrect type'
end select
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR:reader: could not read ', dname, ' from ', self%filename
  error stop
endif

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_read_1d

module procedure hdf_read_2d
integer(HSIZE_T) :: dims(rank(value))
integer(HID_T) :: did, sid, mem_sid
integer :: ier

did = 0 !< sentinel
sid = H5S_ALL_F
mem_sid = H5S_ALL_F
dims = shape(value)

if(present(istart) .and. present(iend)) then
  if(present(stride)) then
    !! necessary to use this present check for Intel and GCC
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend, stride)
  else
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend)
  endif
else
  call hdf_shape_check(self, dname, dims)
  call h5dopen_f(self%lid, dname, did, ier)
  if(ier /= 0) then
    write(stderr, *) 'h5fortran:ERROR:reader: could not setup read ', dname, ' from ', self%filename
    error stop
  endif
endif

select type(value)
type is(real(real64))
  call h5dread_f(did, H5T_NATIVE_DOUBLE, value, dims, ier, mem_sid, sid)
type is(real(real32))
  call h5dread_f(did, H5T_NATIVE_REAL, value, dims, ier, mem_sid, sid)
type is(integer(int32))
  call h5dread_f(did, H5T_NATIVE_INTEGER, value, dims, ier, mem_sid, sid)
class default
  error stop 'h5fortran:reader: incorrect type'
end select
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR:reader: could not read ', dname, ' from ', self%filename
  error stop
endif

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_read_2d

module procedure hdf_read_3d
integer(HSIZE_T) :: dims(rank(value))
integer(HID_T) :: did, sid, mem_sid
integer :: ier

did = 0 !< sentinel
sid = H5S_ALL_F
mem_sid = H5S_ALL_F
dims = shape(value)

if(present(istart) .and. present(iend)) then
  if(present(stride)) then
    !! necessary to use this present check for Intel and GCC
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend, stride)
  else
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend)
  endif
else
  call hdf_shape_check(self, dname, dims)
  call h5dopen_f(self%lid, dname, did, ier)
  if(ier /= 0) then
    write(stderr, *) 'h5fortran:ERROR:reader: could not setup read ', dname, ' from ', self%filename
    error stop
  endif
endif

select type(value)
type is(real(real64))
  call h5dread_f(did, H5T_NATIVE_DOUBLE, value, dims, ier, mem_sid, sid)
type is(real(real32))
  call h5dread_f(did, H5T_NATIVE_REAL, value, dims, ier, mem_sid, sid)
type is(integer(int32))
  call h5dread_f(did, H5T_NATIVE_INTEGER, value, dims, ier, mem_sid, sid)
class default
  error stop 'h5fortran:reader: incorrect type'
end select
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR:reader: could not read ', dname, ' from ', self%filename
  error stop
endif

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_read_3d

module procedure hdf_read_4d
integer(HSIZE_T) :: dims(rank(value))
integer(HID_T) :: did, sid, mem_sid
integer :: ier

did = 0 !< sentinel
sid = H5S_ALL_F
mem_sid = H5S_ALL_F
dims = shape(value)

if(present(istart) .and. present(iend)) then
  if(present(stride)) then
    !! necessary to use this present check for Intel and GCC
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend, stride)
  else
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend)
  endif
else
  call hdf_shape_check(self, dname, dims)
  call h5dopen_f(self%lid, dname, did, ier)
  if(ier /= 0) then
    write(stderr, *) 'h5fortran:ERROR:reader: could not setup read ', dname, ' from ', self%filename
    error stop
  endif
endif

select type(value)
type is(real(real64))
  call h5dread_f(did, H5T_NATIVE_DOUBLE, value, dims, ier, mem_sid, sid)
type is(real(real32))
  call h5dread_f(did, H5T_NATIVE_REAL, value, dims, ier, mem_sid, sid)
type is(integer(int32))
  call h5dread_f(did, H5T_NATIVE_INTEGER, value, dims, ier, mem_sid, sid)
class default
  error stop 'h5fortran:reader: incorrect type'
end select
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR:reader: could not read ', dname, ' from ', self%filename
  error stop
endif

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_read_4d

module procedure hdf_read_5d
integer(HSIZE_T) :: dims(rank(value))
integer(HID_T) :: did, sid, mem_sid
integer :: ier

did = 0 !< sentinel
sid = H5S_ALL_F
mem_sid = H5S_ALL_F
dims = shape(value)

if(present(istart) .and. present(iend)) then
  if(present(stride)) then
    !! necessary to use this present check for Intel and GCC
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend, stride)
  else
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend)
  endif
else
  call hdf_shape_check(self, dname, dims)
  call h5dopen_f(self%lid, dname, did, ier)
  if(ier /= 0) then
    write(stderr, *) 'h5fortran:ERROR:reader: could not setup read ', dname, ' from ', self%filename
    error stop
  endif
endif

select type(value)
type is(real(real64))
  call h5dread_f(did, H5T_NATIVE_DOUBLE, value, dims, ier, mem_sid, sid)
type is(real(real32))
  call h5dread_f(did, H5T_NATIVE_REAL, value, dims, ier, mem_sid, sid)
type is(integer(int32))
  call h5dread_f(did, H5T_NATIVE_INTEGER, value, dims, ier, mem_sid, sid)
class default
  error stop 'h5fortran:reader: incorrect type'
end select
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR:reader: could not read ', dname, ' from ', self%filename
  error stop
endif

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_read_5d

module procedure hdf_read_6d
integer(HSIZE_T) :: dims(rank(value))
integer(HID_T) :: did, sid, mem_sid
integer :: ier

did = 0 !< sentinel
sid = H5S_ALL_F
mem_sid = H5S_ALL_F
dims = shape(value)

if(present(istart) .and. present(iend)) then
  if(present(stride)) then
    !! necessary to use this present check for Intel and GCC
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend, stride)
  else
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend)
  endif
else
  call hdf_shape_check(self, dname, dims)
  call h5dopen_f(self%lid, dname, did, ier)
  if(ier /= 0) then
    write(stderr, *) 'h5fortran:ERROR:reader: could not setup read ', dname, ' from ', self%filename
    error stop
  endif
endif

select type(value)
type is(real(real64))
  call h5dread_f(did, H5T_NATIVE_DOUBLE, value, dims, ier, mem_sid, sid)
type is(real(real32))
  call h5dread_f(did, H5T_NATIVE_REAL, value, dims, ier, mem_sid, sid)
type is(integer(int32))
  call h5dread_f(did, H5T_NATIVE_INTEGER, value, dims, ier, mem_sid, sid)
class default
  error stop 'h5fortran:reader: incorrect type'
end select
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR:reader: could not read ', dname, ' from ', self%filename
  error stop
endif

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_read_6d

module procedure hdf_read_7d
integer(HSIZE_T) :: dims(rank(value))
integer(HID_T) :: did, sid, mem_sid
integer :: ier

did = 0 !< sentinel
sid = H5S_ALL_F
mem_sid = H5S_ALL_F
dims = shape(value)

if(present(istart) .and. present(iend)) then
  if(present(stride)) then
    !! necessary to use this present check for Intel and GCC
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend, stride)
  else
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend)
  endif
else
  call hdf_shape_check(self, dname, dims)
  call h5dopen_f(self%lid, dname, did, ier)
  if(ier /= 0) then
    write(stderr, *) 'h5fortran:ERROR:reader: could not setup read ', dname, ' from ', self%filename
    error stop
  endif
endif

select type(value)
type is(real(real64))
  call h5dread_f(did, H5T_NATIVE_DOUBLE, value, dims, ier, mem_sid, sid)
type is(real(real32))
  call h5dread_f(did, H5T_NATIVE_REAL, value, dims, ier, mem_sid, sid)
type is(integer(int32))
  call h5dread_f(did, H5T_NATIVE_INTEGER, value, dims, ier, mem_sid, sid)
class default
  error stop 'h5fortran:reader: incorrect type'
end select
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR:reader: could not read ', dname, ' from ', self%filename
  error stop
endif

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_read_7d

end submodule reader
