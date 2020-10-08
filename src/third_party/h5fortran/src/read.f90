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

submodule(h5fortran) read
!! This submodule is for reading HDF5 via submodules
use hdf5, only: h5dget_create_plist_f, &
                h5pget_layout_f, h5pget_chunk_f
use H5LT, only: h5ltpath_valid_f

implicit none(type, external)

contains

module procedure hdf_get_ndims
!! get rank or "ndims"
integer :: ier

if(.not. self%is_open) error stop 'h5fortran:read: file handle is not open'

drank = -1

if(self%exist(dname)) then
  call h5ltget_dataset_ndims_f(self%lid, dname, drank, ier)
else
  write(stderr, *) 'ERROR:get_ndims: '//dname//' does not exist in '//self%filename
endif

end procedure hdf_get_ndims

module procedure hdf_get_shape
!! must get dims before info, as "dims" must be allocated or segfault occurs.
integer(SIZE_T) :: dsize
integer :: dtype, drank, ier

if(.not. self%is_open) error stop 'h5fortran:get_shape: file handle is not open'

ier = -1

if(self%exist(dname)) then
  call h5ltget_dataset_ndims_f(self%lid, dname, drank, ier)
else
  write(stderr, *) 'ERROR:get_shape: '//dname//' does not exist in '//self%filename
endif

if(ier == 0) then
  allocate(dims(drank))
  call h5ltget_dataset_info_f(self%lid, dname, dims, dtype, dsize, ier)
endif

if(present(ierr)) ierr = ier
if(ier /= 0) then
  if(present(ierr)) return
  error stop
endif

end procedure hdf_get_shape

module procedure hdf_get_chunk

integer :: ierr, drank
integer(HID_T) :: pid, did

if(.not. self%is_open) error stop 'h5fortran:read: file handle is not open'

chunk_size = -1
if(.not. self%exist(dname)) then
  write(stderr, *) 'ERROR:get_chunk: '//dname//' does not exist in '//self%filename
  ierr = -1
  return
endif

if(.not. self%is_chunked(dname)) return

call h5ltget_dataset_ndims_f(self%lid, dname, drank, ierr)
if(check(ierr, 'ERROR:get_chunk: get rank '//dname//' '//self%filename)) return
call h5dopen_f(self%lid, dname, did, ierr)
if(check(ierr, 'ERROR:get_chunk: open dataset '//dname//' '//self%filename)) return
call h5dget_create_plist_f(did, pid, ierr)
if(check(ierr, 'ERROR:get_chunk: get property list ID '//dname//' '//self%filename)) return

call h5pget_chunk_f(pid, drank, chunk_size, ierr)
if(ierr /= drank) then
  write(stderr, *) 'ERROR:get_chunk read '//dname//' '//self%filename
  return
endif

call h5dclose_f(did, ierr)
if(check(ierr, 'ERROR:get_chunk: close dataset: '//dname//' '//self%filename)) return

end procedure hdf_get_chunk

module procedure hdf_get_layout

integer(HID_T) :: pid, did
integer :: ierr

if(.not. self%is_open) error stop 'h5fortran:read: file handle is not open'

layout = -1

if(.not. self%exist(dname)) then
  write(stderr, *) 'ERROR:get_layout: '//dname//' does not exist in '//self%filename
  return
endif

call h5dopen_f(self%lid, dname, did, ierr)
if(check(ierr, 'ERROR:get_layout: open dataset '//dname//' '//self%filename)) return
call h5dget_create_plist_f(did, pid, ierr)
if(check(ierr, 'ERROR:get_layout: get property list ID '//dname//' '//self%filename)) return
call h5pget_layout_f(pid, layout, ierr)
if(check(ierr, 'ERROR:get_layout read '//dname//' '//self%filename)) return
call h5dclose_f(did, ierr)
if(check(ierr, 'ERROR:get_layout: close dataset: '//dname//' '//self%filename)) return

end procedure hdf_get_layout

module procedure hdf_check_exist

integer :: ierr

exists = .false.

if(.not. self%is_open) error stop 'h5fortran:exist: file handle is not open'

call h5ltpath_valid_f(self%lid, dname, .true., exists, ierr)
!! h5lexists_f can false error with groups--just use h5ltpath_valid

if(ierr /= 0) then
  write(stderr, *) 'ERROR:h5fortran:check_exist: could not determine status of '//dname//' in '//self%filename
  return
endif

end procedure hdf_check_exist

end submodule read
