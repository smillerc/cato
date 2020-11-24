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

submodule(h5fortran) pathlib
!! vendored from Michael Hirsch's Fortran pathlib

implicit none

contains

module procedure get_tempdir

character(1024) :: argv
integer :: L, i

call get_environment_variable("TMP", argv, L, i)
if(L == 0 .or. i /= 0) call get_environment_variable("TEMP", argv, L, i)
if(L == 0 .or. i /= 0) call get_environment_variable("TMPDIR", argv, L, i)
if(L == 0 .or. i /= 0) argv = "/tmp"

get_tempdir = trim(argv)

end procedure get_tempdir

module procedure is_absolute_path
!! heuristic to determine if is absolute path

is_absolute_path = .false.

if(.false.) then
  if(lge(path(1:1), 'A') .and. lle(path(1:1), 'z') .and. path(2:2) == ':') is_absolute_path = .true.
else
  if(path(1:1) == '/') is_absolute_path = .true.
endif

end procedure is_absolute_path

module procedure unlink
!! deletes file in Fortran standard manner.
integer :: i, u

open(newunit=u, file=filename, iostat=i)
close(u, status='delete', iostat=i)

inquire(file=filename, exist=unlink)

end procedure unlink

end submodule pathlib
