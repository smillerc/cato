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

submodule(h5fortran:read) reader_lt

implicit none(type, external)

contains

module procedure h5exist

type(hdf5_file) :: h

call h%initialize(filename, status='old', action='r')
h5exist = h%exist(dname)
call h%finalize()

end procedure h5exist

module procedure lt0read
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='old', action='r')
if(ier == 0) call h%read(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt0read

module procedure lt1read
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='old', action='r')
if(ier == 0) call h%read(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt1read

module procedure lt2read
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='old', action='r')
if(ier == 0) call h%read(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt2read

module procedure lt3read
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='old', action='r')
if(ier == 0) call h%read(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt3read

module procedure lt4read
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='old', action='r')
if(ier == 0) call h%read(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt4read

module procedure lt5read
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='old', action='r')
if(ier == 0) call h%read(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt5read

module procedure lt6read
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='old', action='r')
if(ier == 0) call h%read(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt6read

module procedure lt7read
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='old', action='r')
if(ier == 0) call h%read(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt7read

end submodule reader_lt
