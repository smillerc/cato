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

submodule(h5fortran:write) writer_lt

implicit none(type, external)

contains

module procedure lt0write
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='unknown')
if(ier == 0) call h%write(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt0write

module procedure lt1write
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='unknown')
if(ier == 0) call h%write(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt1write

module procedure lt2write
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='unknown')
if(ier == 0) call h%write(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt2write

module procedure lt3write
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='unknown')
if(ier == 0) call h%write(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt3write

module procedure lt4write
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='unknown')
if(ier == 0) call h%write(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt4write

module procedure lt5write
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='unknown')
if(ier == 0) call h%write(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt5write

module procedure lt6write
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='unknown')
if(ier == 0) call h%write(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt6write

module procedure lt7write
type(hdf5_file) :: h
integer :: ier

call h%initialize(filename, ier, status='unknown')
if(ier == 0) call h%write(dname, value, ier)
if(ier == 0) call h%finalize(ier)

if(present(ierr)) ierr = ier
if(check(ier, filename, dname) .and. .not. present(ierr)) error stop

end procedure lt7write

end submodule writer_lt
