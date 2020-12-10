add_subdirectory(ini_parser)
add_subdirectory(finterp)
add_subdirectory(functional)
add_subdirectory(h5fortran)

if(USE_OPENCL)
  add_subdirectory(focal)
endif()
