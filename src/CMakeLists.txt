add_library(ini_parser
            third_party/ini_parser/cfgio_mod.f90
            third_party/ini_parser/string_conv_mod.f90
            )

add_library(constants
            lib/math_constants.f90
            )

add_library(io
            lib/io/mod_input.f90
            )
target_link_libraries(io PUBLIC ini_parser)


add_library(utilities
            lib/limiters/slope_limiter.f90
            )

add_library(geometry
            lib/domain/grid_type.f90
            lib/domain/regular_2d_grid.f90
            lib/domain/grid_factory.f90
            lib/geometry/vector_2d.f90
            lib/geometry/quad_cell.f90
            lib/geometry/mach_cone_geometry.f90
            )

target_link_libraries(geometry PUBLIC constants io)

add_library(fvleg
            lib/domain/conserved_vars.f90
            lib/evolution_operator/abstract_evo_operator.f90
            lib/evolution_operator/local_evo_operator.f90
            lib/reconstruction/abstract_reconstruction.f90
            lib/reconstruction/first_order_reconstruction.f90
            lib/reconstruction/second_order_reconstruction.f90
            lib/reconstruction/reconstruction_factory.f90
            lib/finite_volume_scheme/fv_scheme.f90
            lib/finite_volume_scheme/fvleg.f90
            lib/integrand.f90
            lib/integrators/strategy.f90
            lib/integrators/surrogate.f90
            lib/integrators/runge_kutta_2nd_order.f90
            )

# include_directories(${HDF5_INCLUDE_DIRS})
# target_include_directories(fvleg PUBLIC ${HDF5_INCLUDE_DIRS})

# target_compile_options(fvleg PUBLIC ${Coarray_COMPILE_OPTIONS})
target_link_libraries(fvleg
                      PUBLIC
                      constants
                      geometry
                      utilities
                      io
                      #  functional
                      ${Coarray_LIBRARIES}
                      ${HDF5_LIBRARIES}
                      ${HDF5_HL_LIBRARIES})

add_executable(fvleg_2d.x main.f90)
target_link_libraries(fvleg_2d.x
                      io
                      fvleg
                      #   functional
                      #   ${HDF5_LIBRARIES}
                      #   ${HDF5_HL_LIBRARIES}
                      ${Coarray_LIBRARIES}
                      )
# target_compile_options(fvleg_2d.x PUBLIC ${Coarray_COMPILE_OPTIONS})



# set_target_properties (fvleg_2d.x PROPERTIES
#   Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
