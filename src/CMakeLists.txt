add_library(constants
            lib/math_constants.f90
            )

add_library(geometry
            lib/geometry/quad_cell.f90
            lib/mach_cone_geometry.f90
            )

target_link_libraries(geometry PUBLIC constants)

add_library(fvleg
            lib/reconstruction.f90
            lib/mach_cone_geometry.f90
            lib/slope_limiter.f90
            )

include_directories(${HDF5_INCLUDE_DIRS})
target_include_directories(fvleg PUBLIC ${HDF5_INCLUDE_DIRS})

target_compile_options(fvleg PUBLIC ${Coarray_COMPILE_OPTIONS})
target_link_libraries(fvleg
                      PUBLIC 
                      constants
                      geometry
                      ${Coarray_LIBRARIES}
                            #  functional
                             ${HDF5_LIBRARIES}
                             ${HDF5_HL_LIBRARIES})

add_executable(fvleg_2d.x main.f90)

set_target_properties (fvleg_2d.x PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

target_compile_options(fvleg_2d.x PUBLIC ${Coarray_COMPILE_OPTIONS})
target_link_libraries(fvleg_2d.x
                      fvleg
                    #   functional
                    #   ${HDF5_LIBRARIES}
                    #   ${HDF5_HL_LIBRARIES}
                      ${Coarray_LIBRARIES}
                      )
