! The MIT License
! ===============

! Copyright (c) 2016 Stefano Zaghi

! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.

!< VTK_Fortran, pure Fortran (2003+) library to parse and emitt VTK files.
module vtk_fortran
  use penf
  use vtk_fortran_pvtk_file, only: pvtk_file
  use vtk_fortran_vtk_file, only: vtk_file
  use vtk_fortran_vtm_file, only: vtm_file

  implicit none
  private
  public :: pvtk_file
  public :: vtk_file
  public :: vtm_file
  public :: write_xml_volatile

contains
  function write_xml_volatile(xml_volatile, filename) result(error)
    !< Write the volatile file into a real file.
    !< This is what a master process should do into a parallel scenario where it being the only process allowed to access to
    !< filesystem: slave processes create XML volatile file econded into a characters string and master process collects and writes
    !< them by means of `write_xml_volatile`.
    character(*), intent(in) :: xml_volatile !< XML volatile file.
    character(*), intent(in) :: filename     !< XML file name.
    integer(I4P)             :: error        !< Status error.
    integer(I4P)             :: xml_unit     !< XML file unit.

    open(newunit=xml_unit, &
         file=trim(adjustl(filename)), &
         form='UNFORMATTED', &
         access='STREAM', &
         action='WRITE', &
         status='REPLACE', &
         iostat=error)
    write(unit=xml_unit, iostat=error) xml_volatile
  endfunction write_xml_volatile
endmodule vtk_fortran
