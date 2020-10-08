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

submodule(h5fortran) write
!! This submodule is for writing HDF5 data via child submodules
use hdf5, only: &
  h5screate_f, H5S_SCALAR_F, &
  h5dcreate_f, &
  h5pset_chunk_f, h5pset_deflate_f, h5pset_shuffle_f, h5pset_fletcher32_f, h5pcreate_f, H5P_DATASET_CREATE_F, h5pclose_f, &
  h5gopen_f

use H5LT, only: h5ltmake_dataset_string_f, h5ltpath_valid_f

implicit none(type, external)

contains

module procedure hdf_create

logical :: exists
integer :: ierr
integer(HID_T) :: pid, space_id, ds_id
integer(HSIZE_T) :: vdims(size(dims))

if(.not. self%is_open) error stop 'h5fortran:write: file handle is not open'

call h5ltpath_valid_f(self%lid, dname, .true., exists, ierr)
!! h5lexists_f can false error with groups--just use h5ltpath_valid

!! stricter than self%exists() since we're creating and/or writing variable
if(ierr /= 0) then
  write(stderr, *) 'ERROR:h5fortran:create: variable path invalid: ', dname, ' in ', self%filename
  error stop 6
endif

!> allow user to specify int4 or int8 dims
select type(dims)
type is(integer(int32))
  vdims = int(dims, int64)
type is(integer(hsize_t))
  vdims = dims
class default
  write(stderr, *) 'ERROR:h5fortran:create: wrong type for dims: ', dname, self%filename
  error stop 5
end select

if(self%debug) print *, 'h5fortran:TRACE:create:exists: '//dname, exists

if(exists) then
  if(.not. present(istart)) call hdf_shape_check(self, dname, vdims)
  !! FIXME: read and write slice shape not checked; but should check in future versions
  !> open dataset
  call h5dopen_f(self%lid, dname, ds_id, ierr)
  if(ierr /= 0) then
    write(stderr, *) 'ERROR:h5fortran:create: could not open ', dname, ' in ', self%filename
    error stop
  endif

  if(present(did)) did = ds_id
  if(present(sid)) then
    call h5dget_space_f(ds_id, sid, ierr)
    if(ierr /= 0) error stop 'h5fortran:create could not get dataset'
  endif
  return
endif

if(self%debug) print *, 'h5fortran:TRACE1: '//dname

!> Only new datasets go past this point

call self%write_group(dname, ierr)
if(ierr /= 0) then
  write(stderr, *) 'ERROR:h5fortran:create: could not create group for ', dname, ' in ', self%filename
  error stop
endif

if(size(vdims) >= 2) then
  if(self%debug) print *, 'h5fortran:TRACE:create: deflate: '//dname
  call set_deflate(self, vdims, pid, ierr, chunk_size)
  if(ierr /= 0) error stop 'ERROR:h5fortran:create: problem setting deflate on'
else
  pid = 0
endif

if(size(vdims) == 0) then
  call h5screate_f(H5S_SCALAR_F, space_id, ierr)
else
  call h5screate_simple_f(size(vdims), vdims, space_id, ierr)
endif
if(check(ierr, self%filename, dname)) error stop

if(pid == 0) then
  call h5dcreate_f(self%lid, dname, dtype, space_id, ds_id, ierr)
else
  call h5dcreate_f(self%lid, dname, dtype, space_id, ds_id, ierr, pid)
  if(check(ierr, self%filename, dname)) error stop
  call h5pclose_f(pid, ierr)
  if(check(ierr, self%filename, dname)) error stop
endif
if(check(ierr, self%filename, dname)) error stop

if(.not.(present(did) .and. present(sid))) then
  if(self%debug) print *, 'h5fortran:TRACE:create: closing dataset ', dname
  call hdf_wrapup(ds_id, space_id, ierr)
endif
if(check(ierr, self%filename, dname)) error stop

if(present(sid)) sid = space_id
if(present(did)) did = ds_id

end procedure hdf_create

subroutine set_deflate(self, dims, pid, ierr, chunk_size)
  class(hdf5_file), intent(inout) :: self
  integer(HSIZE_T), intent(in) :: dims(:)
  integer(HID_T), intent(out) :: pid
  integer, intent(out) :: ierr
  integer, intent(in), optional :: chunk_size(:)

  integer(HSIZE_T) :: cs(size(dims))

  pid = 0
  if(self%comp_lvl < 1 .or. self%comp_lvl > 9) return

  if(present(chunk_size)) then
    cs = chunk_size
    where(cs > dims) cs = dims
    if(self%debug) print *, 'TRACE: user request chunk_size ', cs
  else
  !! guess chunk size, keeping in mind 1 Megabyte recommended maximum chunk size
    call guess_chunk_size(dims, cs)
  endif

  if(any(cs < 1)) return

  if(self%debug) print *, 'DEBUG:set_deflate: dims: ', dims, 'chunk size: ', cs

  call h5pcreate_f(H5P_DATASET_CREATE_F, pid, ierr)
  if(check(ierr, self%filename)) return

  call h5pset_chunk_f(pid, size(dims), cs, ierr)
  if(check(ierr, self%filename)) return

  call h5pset_shuffle_f(pid, ierr)
  if(check(ierr, self%filename)) return

  call h5pset_fletcher32_f(pid, ierr)
  if(check(ierr, self%filename)) return

  call h5pset_deflate_f(pid, self%comp_lvl, ierr)
  if(check(ierr, self%filename)) return

  if(self%debug) print *, 'TRACE:set_deflate done'

end subroutine set_deflate

subroutine guess_chunk_size(dims, chunk_size)
!! based on https://github.com/h5py/h5py/blob/master/h5py/_hl/filters.py
!! refer to https://support.hdfgroup.org/HDF5/Tutor/layout.html
  integer(HSIZE_T), intent(in) :: dims(:)
  integer(HSIZE_T), intent(out) :: chunk_size(:)

  integer(hsize_t), parameter :: &
    CHUNK_BASE = 16000, &    !< Multiplier by which chunks are adjusted
    CHUNK_MIN = 8000, &      !< lower limit: 8 kbyte
    CHUNK_MAX = 1000000, &   !< upper limit: 1 Mbyte
    TYPESIZE = 8             !< bytes, assume real64 for simplicity

  integer(hsize_t) :: dset_size, target_size, chunk_bytes, i, j, ndims

  if(product(dims) * TYPESIZE < CHUNK_MIN) then
    chunk_size = 0
    return
  endif

  ndims = size(chunk_size)
  chunk_size = dims

  dset_size = product(chunk_size) * TYPESIZE
  target_size = int(CHUNK_BASE * (2**log10(real(dset_size) / 1e6)), hsize_t)
  if(target_size > CHUNK_MAX) target_size = CHUNK_MAX

! print *,'target_size [bytes]: ',target_size

  i = 0
  do
  !! Repeatedly loop over the axes, dividing them by 2.
  !! Stop when:
  !!   1a. We're smaller than the target chunk size, OR
  !!   1b. We're within 50% of the target chunk size, AND
  !!    2. The chunk is smaller than the maximum chunk size

    chunk_bytes = product(chunk_size) * TYPESIZE

    if((chunk_bytes < target_size .or. 2 * (abs(chunk_bytes - target_size) / target_size) < 1) .and. &
       chunk_bytes < CHUNK_MAX) exit

    if(product(chunk_size) == 1) exit
  !! Element size larger than CHUNK_MAX
    j = int(modulo(i, ndims), hsize_t) + 1
    if(j < 1 .or. j > ndims) error stop 'auto index bounds error'
    chunk_size(j) = ceiling(real(chunk_size(j)) / 2.0)
    i = i + 1
  end do

end subroutine guess_chunk_size

module procedure hdf_open_group

integer :: ier

if(.not. self%is_open) error stop 'h5fortran:write: file handle is not open'

call h5gopen_f(self%lid, gname, self%gid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename, gname)) then
  if(present(ierr)) return
  error stop
endif

self%glid = self%lid
self%lid = self%gid

end procedure hdf_open_group

module procedure hdf_close_group

integer :: ier

if(.not. self%is_open) error stop 'h5fortran:write: file handle is not open'

call h5gclose_f(self%gid, ier)

if(present(ierr)) ierr = ier
if(check(ier, self%filename)) then
  if(present(ierr)) return
  error stop
endif

self%lid = self%glid

end procedure hdf_close_group

end submodule write
