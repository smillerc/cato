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

submodule(h5fortran:write) writer
!! This submodule is for writing 0-D..7-D data
use hdf5, only: h5dwrite_f

implicit none(type, external)

contains

module procedure hdf_write_scalar

integer(HID_T)  :: sid, did
integer(HSIZE_T), allocatable :: dims(:)
integer :: ier

if(.not. self%is_open) error stop 'h5fortran:writer: file handle is not open'

allocate(dims(0))

did = 0
sid = 0
!! sentinel

select type(value)
type is(character(*))
  call h5ltmake_dataset_string_f(self%lid, dname, value, ier)
type is(real(real64))
  !! NOTE: 0d does not use chunk_size
  call hdf_create(self, dname, H5T_NATIVE_DOUBLE, dims, sid, did)
  call h5dwrite_f(did, H5T_NATIVE_DOUBLE, value, dims, ier)
type is(real(real32))
  call hdf_create(self, dname, H5T_NATIVE_REAL, dims, sid, did)
  call h5dwrite_f(did, H5T_NATIVE_REAL, value, dims, ier)
type is(integer(int32))
  call hdf_create(self, dname, H5T_NATIVE_INTEGER, dims, sid, did)
  call h5dwrite_f(did, H5T_NATIVE_INTEGER, value, dims, ier)
class default
  error stop 'h5fortran: invalid data type'
end select
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR: could not write ', dname, ' to ', self%filename
  error stop
endif

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_write_scalar

module procedure hdf_write_1d
integer(HID_T) :: did, sid, mem_sid
integer(HSIZE_T) :: dims(rank(value))
integer :: ier

dims = shape(value)
select type(value)
type is(real(real64))
  call hdf_create(self, dname, H5T_NATIVE_DOUBLE, dims, sid, did, chunk_size, istart, iend, stride)
type is(real(real32))
  call hdf_create(self, dname, H5T_NATIVE_REAL, dims, sid, did, chunk_size, istart, iend, stride)
type is(integer(int32))
  call hdf_create(self, dname, H5T_NATIVE_INTEGER, dims, sid, did, chunk_size, istart, iend, stride)
class default
  error stop 'h5fortran:write:invalid data type'
end select

mem_sid = H5S_ALL_F !< default

if(present(istart) .and. present(iend)) then
  if(present(stride)) then
    !! necessary to use this present check for Intel and GCC
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend, stride)
  else
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend)
  endif
endif

select type(value)
type is(real(real64))
  call h5dwrite_f(did, H5T_NATIVE_DOUBLE, value, dims, ier, mem_sid, sid)
type is(real(real32))
  call h5dwrite_f(did, H5T_NATIVE_REAL, value, dims, ier, mem_sid, sid)
type is(integer(int32))
  call h5dwrite_f(did, H5T_NATIVE_INTEGER, value, dims, ier, mem_sid, sid)
end select
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR: could not write ', dname, ' to ', self%filename
  error stop
endif

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_write_1d

module procedure hdf_write_2d
integer(HID_T) :: did, sid, mem_sid
integer(HSIZE_T) :: dims(rank(value))
integer :: ier

dims = shape(value)
select type(value)
type is(real(real64))
  call hdf_create(self, dname, H5T_NATIVE_DOUBLE, dims, sid, did, chunk_size, istart, iend, stride)
type is(real(real32))
  call hdf_create(self, dname, H5T_NATIVE_REAL, dims, sid, did, chunk_size, istart, iend, stride)
type is(integer(int32))
  call hdf_create(self, dname, H5T_NATIVE_INTEGER, dims, sid, did, chunk_size, istart, iend, stride)
class default
  error stop 'h5fortran:write:invalid data type'
end select

mem_sid = H5S_ALL_F !< default

if(present(istart) .and. present(iend)) then
  if(present(stride)) then
    !! necessary to use this present check for Intel and GCC
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend, stride)
  else
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend)
  endif
endif

select type(value)
type is(real(real64))
  call h5dwrite_f(did, H5T_NATIVE_DOUBLE, value, dims, ier, mem_sid, sid)
type is(real(real32))
  call h5dwrite_f(did, H5T_NATIVE_REAL, value, dims, ier, mem_sid, sid)
type is(integer(int32))
  call h5dwrite_f(did, H5T_NATIVE_INTEGER, value, dims, ier, mem_sid, sid)
end select
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR: could not write ', dname, ' to ', self%filename
  error stop
endif

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_write_2d

module procedure hdf_write_3d
integer(HID_T) :: did, sid, mem_sid
integer(HSIZE_T) :: dims(rank(value))
integer :: ier

dims = shape(value)
select type(value)
type is(real(real64))
  call hdf_create(self, dname, H5T_NATIVE_DOUBLE, dims, sid, did, chunk_size, istart, iend, stride)
type is(real(real32))
  call hdf_create(self, dname, H5T_NATIVE_REAL, dims, sid, did, chunk_size, istart, iend, stride)
type is(integer(int32))
  call hdf_create(self, dname, H5T_NATIVE_INTEGER, dims, sid, did, chunk_size, istart, iend, stride)
class default
  error stop 'h5fortran:write:invalid data type'
end select

mem_sid = H5S_ALL_F !< default

if(present(istart) .and. present(iend)) then
  if(present(stride)) then
    !! necessary to use this present check for Intel and GCC
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend, stride)
  else
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend)
  endif
endif

select type(value)
type is(real(real64))
  call h5dwrite_f(did, H5T_NATIVE_DOUBLE, value, dims, ier, mem_sid, sid)
type is(real(real32))
  call h5dwrite_f(did, H5T_NATIVE_REAL, value, dims, ier, mem_sid, sid)
type is(integer(int32))
  call h5dwrite_f(did, H5T_NATIVE_INTEGER, value, dims, ier, mem_sid, sid)
end select
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR: could not write ', dname, ' to ', self%filename
  error stop
endif

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_write_3d

module procedure hdf_write_4d
integer(HID_T) :: did, sid, mem_sid
integer(HSIZE_T) :: dims(rank(value))
integer :: ier

dims = shape(value)
select type(value)
type is(real(real64))
  call hdf_create(self, dname, H5T_NATIVE_DOUBLE, dims, sid, did, chunk_size, istart, iend, stride)
type is(real(real32))
  call hdf_create(self, dname, H5T_NATIVE_REAL, dims, sid, did, chunk_size, istart, iend, stride)
type is(integer(int32))
  call hdf_create(self, dname, H5T_NATIVE_INTEGER, dims, sid, did, chunk_size, istart, iend, stride)
class default
  error stop 'h5fortran:write:invalid data type'
end select

mem_sid = H5S_ALL_F !< default

if(present(istart) .and. present(iend)) then
  if(present(stride)) then
    !! necessary to use this present check for Intel and GCC
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend, stride)
  else
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend)
  endif
endif

select type(value)
type is(real(real64))
  call h5dwrite_f(did, H5T_NATIVE_DOUBLE, value, dims, ier, mem_sid, sid)
type is(real(real32))
  call h5dwrite_f(did, H5T_NATIVE_REAL, value, dims, ier, mem_sid, sid)
type is(integer(int32))
  call h5dwrite_f(did, H5T_NATIVE_INTEGER, value, dims, ier, mem_sid, sid)
end select
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR: could not write ', dname, ' to ', self%filename
  error stop
endif

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_write_4d

module procedure hdf_write_5d
integer(HID_T) :: did, sid, mem_sid
integer(HSIZE_T) :: dims(rank(value))
integer :: ier

dims = shape(value)
select type(value)
type is(real(real64))
  call hdf_create(self, dname, H5T_NATIVE_DOUBLE, dims, sid, did, chunk_size, istart, iend, stride)
type is(real(real32))
  call hdf_create(self, dname, H5T_NATIVE_REAL, dims, sid, did, chunk_size, istart, iend, stride)
type is(integer(int32))
  call hdf_create(self, dname, H5T_NATIVE_INTEGER, dims, sid, did, chunk_size, istart, iend, stride)
class default
  error stop 'h5fortran:write:invalid data type'
end select

mem_sid = H5S_ALL_F !< default

if(present(istart) .and. present(iend)) then
  if(present(stride)) then
    !! necessary to use this present check for Intel and GCC
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend, stride)
  else
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend)
  endif
endif

select type(value)
type is(real(real64))
  call h5dwrite_f(did, H5T_NATIVE_DOUBLE, value, dims, ier, mem_sid, sid)
type is(real(real32))
  call h5dwrite_f(did, H5T_NATIVE_REAL, value, dims, ier, mem_sid, sid)
type is(integer(int32))
  call h5dwrite_f(did, H5T_NATIVE_INTEGER, value, dims, ier, mem_sid, sid)
end select
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR: could not write ', dname, ' to ', self%filename
  error stop
endif

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_write_5d

module procedure hdf_write_6d
integer(HID_T) :: did, sid, mem_sid
integer(HSIZE_T) :: dims(rank(value))
integer :: ier

dims = shape(value)
select type(value)
type is(real(real64))
  call hdf_create(self, dname, H5T_NATIVE_DOUBLE, dims, sid, did, chunk_size, istart, iend, stride)
type is(real(real32))
  call hdf_create(self, dname, H5T_NATIVE_REAL, dims, sid, did, chunk_size, istart, iend, stride)
type is(integer(int32))
  call hdf_create(self, dname, H5T_NATIVE_INTEGER, dims, sid, did, chunk_size, istart, iend, stride)
class default
  error stop 'h5fortran:write:invalid data type'
end select

mem_sid = H5S_ALL_F !< default

if(present(istart) .and. present(iend)) then
  if(present(stride)) then
    !! necessary to use this present check for Intel and GCC
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend, stride)
  else
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend)
  endif
endif

select type(value)
type is(real(real64))
  call h5dwrite_f(did, H5T_NATIVE_DOUBLE, value, dims, ier, mem_sid, sid)
type is(real(real32))
  call h5dwrite_f(did, H5T_NATIVE_REAL, value, dims, ier, mem_sid, sid)
type is(integer(int32))
  call h5dwrite_f(did, H5T_NATIVE_INTEGER, value, dims, ier, mem_sid, sid)
end select
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR: could not write ', dname, ' to ', self%filename
  error stop
endif

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_write_6d

module procedure hdf_write_7d
integer(HID_T) :: did, sid, mem_sid
integer(HSIZE_T) :: dims(rank(value))
integer :: ier

dims = shape(value)
select type(value)
type is(real(real64))
  call hdf_create(self, dname, H5T_NATIVE_DOUBLE, dims, sid, did, chunk_size, istart, iend, stride)
type is(real(real32))
  call hdf_create(self, dname, H5T_NATIVE_REAL, dims, sid, did, chunk_size, istart, iend, stride)
type is(integer(int32))
  call hdf_create(self, dname, H5T_NATIVE_INTEGER, dims, sid, did, chunk_size, istart, iend, stride)
class default
  error stop 'h5fortran:write:invalid data type'
end select

mem_sid = H5S_ALL_F !< default

if(present(istart) .and. present(iend)) then
  if(present(stride)) then
    !! necessary to use this present check for Intel and GCC
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend, stride)
  else
    call hdf_get_slice(self, dname, did, sid, mem_sid, istart, iend)
  endif
endif

select type(value)
type is(real(real64))
  call h5dwrite_f(did, H5T_NATIVE_DOUBLE, value, dims, ier, mem_sid, sid)
type is(real(real32))
  call h5dwrite_f(did, H5T_NATIVE_REAL, value, dims, ier, mem_sid, sid)
type is(integer(int32))
  call h5dwrite_f(did, H5T_NATIVE_INTEGER, value, dims, ier, mem_sid, sid)
end select
if(ier /= 0) then
  write(stderr, *) 'h5fortran:ERROR: could not write ', dname, ' to ', self%filename
  error stop
endif

call hdf_wrapup(did, sid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, dname) .and. .not. present(ierr)) error stop

end procedure hdf_write_7d

end submodule writer
