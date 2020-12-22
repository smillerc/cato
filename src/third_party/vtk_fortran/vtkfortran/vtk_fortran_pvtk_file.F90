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

!< Parallel (partioned) VTK file class.
module vtk_fortran_pvtk_file
  use befor64
  use penf
  use vtk_fortran_vtk_file_xml_writer_abstract
  use vtk_fortran_vtk_file_xml_writer_ascii_local

  implicit none
  private
  public :: pvtk_file

  type :: pvtk_file
    !< VTK parallel (partioned) file class.
    private
    class(xml_writer_abstract), allocatable, public :: xml_writer !< XML writer.
  contains
    procedure, pass(self) :: initialize !< Initialize file.
    procedure, pass(self) :: finalize   !< Finalize file.
  endtype pvtk_file
contains
  function initialize(self, filename, mesh_topology, mesh_kind, nx1, nx2, ny1, ny2, nz1, nz2) result(error)
    !< Initialize file (writer).
    !<
    !< @note This function must be the first to be called.
    !<
    !<### Supported topologies are:
    !<
    !<- PRectilinearGrid;
    !<- PStructuredGrid;
    !<- PUnstructuredGrid.
    !<
    !<### Example of usage
    !<
    !<```fortran
    !< type(pvtk_file) :: pvtk
    !< integer(I4P)    :: nx1, nx2, ny1, ny2, nz1, nz2
    !< ...
    !< error = pvtk%initialize('XML_RECT_BINARY.pvtr','PRectilinearGrid',nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2)
    !< ...
    !<```
    !< @note The file extension is necessary in the file name. The XML standard has different extensions for each
    !< different topologies (e.g. *pvtr* for rectilinear topology). See the VTK-standard file for more information.
    class(pvtk_file), intent(inout)         :: self          !< VTK file.
    character(*), intent(in)            :: filename      !< File name.
    character(*), intent(in)            :: mesh_topology !< Mesh topology.
    character(*), intent(in)            :: mesh_kind     !< Kind of mesh data: Float64, Float32, ecc.
    integer(I4P), intent(in), optional :: nx1           !< Initial node of x axis.
    integer(I4P), intent(in), optional :: nx2           !< Final node of x axis.
    integer(I4P), intent(in), optional :: ny1           !< Initial node of y axis.
    integer(I4P), intent(in), optional :: ny2           !< Final node of y axis.
    integer(I4P), intent(in), optional :: nz1           !< Initial node of z axis.
    integer(I4P), intent(in), optional :: nz2           !< Final node of z axis.
    integer(I4P)                            :: error         !< Error status.

    if(.not. is_initialized) call penf_init
    if(.not. is_b64_initialized) call b64_init
    if(allocated(self%xml_writer)) deallocate(self%xml_writer)
    allocate(xml_writer_ascii_local :: self%xml_writer)
    error = self%xml_writer%initialize(format='ascii', filename=filename, mesh_topology=mesh_topology, &
                                       nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2, nz1=nz1, nz2=nz2, mesh_kind=mesh_kind)
  endfunction initialize

  function finalize(self) result(error)
    !< Finalize file (writer).
    class(pvtk_file), intent(inout) :: self  !< VTK file.
    integer(I4P)                    :: error !< Error status.

    error = 1
    if(allocated(self%xml_writer)) error = self%xml_writer%finalize()
  endfunction finalize
endmodule vtk_fortran_pvtk_file
