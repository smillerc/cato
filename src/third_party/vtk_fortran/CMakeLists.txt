# project(vtkfortran Fortran C)

# Portable environment
add_library(penf 
            penf/penf_b_size.F90
            penf/penf_global_parameters_variables.F90
            penf/penf_stringify.F90
            penf/penf.F90)

# Encoding library
add_library(befor64 
            befor64/befor64_pack_data_m.F90
            befor64/befor64.F90)
target_link_libraries(befor64 PRIVATE penf)

# ANSI color environment
add_library(face face/face.F90)

# String manipulation
add_library(stringifor
            stringifor/stringifor_string_t.F90
            stringifor/stringifor.F90)
target_link_libraries(stringifor PRIVATE penf befor64 face)

# XML library
add_library(foxy 
            foxy/foxy_xml_file.F90
            foxy/foxy_xml_tag.F90
            foxy/foxy.F90)

target_link_libraries(foxy PRIVATE stringifor)


# VTK output
add_library(vtkfortran
            vtkfortran/vtk_fortran_dataarray_encoder.F90
            vtkfortran/vtk_fortran_parameters.F90
            vtkfortran/vtk_fortran_pvtk_file.F90
            vtkfortran/vtk_fortran_vtk_file_xml_writer_abstract.F90
            vtkfortran/vtk_fortran_vtk_file_xml_writer_appended.F90
            vtkfortran/vtk_fortran_vtk_file_xml_writer_ascii_local.F90
            vtkfortran/vtk_fortran_vtk_file_xml_writer_binary_local.F90
            vtkfortran/vtk_fortran_vtk_file.F90
            vtkfortran/vtk_fortran_vtm_file.F90
            vtkfortran/vtk_fortran.F90)

target_link_libraries(vtkfortran PRIVATE penf foxy)