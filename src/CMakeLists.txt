add_library(globals lib/hermetic.f90 lib/globals.f90 lib/math_constants.f90
                    lib/floating_point_utilities.f90 lib/mod_timing.f90)

add_library(ini_parser third_party/ini_parser/cfgio_mod.f90
                       third_party/ini_parser/string_conv_mod.f90)

add_library(finterp third_party/finterp/linear_interpolation_module.f90)

add_library(functional third_party/functional/mod_functional.f90
                       third_party/functional/mod_interfaces.f90)

add_library(io lib/io/mod_hdf5.f90 lib/io/mod_input.f90)

add_library(eos lib/equation_of_state.f90)

include_directories(${HDF5_INCLUDE_DIRS})
target_include_directories(io PUBLIC ${HDF5_INCLUDE_DIRS})

target_link_libraries(io PUBLIC globals ini_parser ${HDF5_LIBRARIES}
                                ${HDF5_HL_LIBRARIES})
target_link_libraries(eos PUBLIC io)

add_library(utilities lib/limiters/slope_limiter.f90)
target_link_libraries(utilities PUBLIC globals)

add_library(
  geometry
  lib/geometry/cone.f90
  lib/domain/grid_type.f90
  lib/domain/regular_2d_grid.f90
  lib/domain/grid_factory.f90
  lib/geometry/vector_2d.f90
  lib/geometry/vector.f90
  lib/geometry/quad_cell.f90)

target_link_libraries(geometry PUBLIC globals io eos)

add_library(
  cato
  lib/finite_volume_scheme/mod_fluid.f90
  lib/io/contour_writer.f90
  lib/boundary_conditions/boundary_conditions.f90
  lib/boundary_conditions/zero_gradient_bc.f90
  lib/boundary_conditions/periodic_bc.f90
  lib/boundary_conditions/pressure_input_bc.f90
  lib/boundary_conditions/reflection_bc.f90
  lib/boundary_conditions/symmetry_bc.f90
  lib/boundary_conditions/bc_factory.f90
  lib/source_terms/source_factory.f90
  lib/source_terms/mod_source.f90
  lib/source_terms/energy_source.f90
  lib/domain/fluxes.f90
  lib/evolution_operator/abstract_evo_operator.f90
  lib/evolution_operator/local_evo_operator.f90
  lib/evolution_operator/evo_operator_factory.f90
  lib/reconstruction/abstract_reconstruction.f90
  lib/reconstruction/first_order_reconstruction.f90
  lib/reconstruction/second_order_reconstruction.f90
  lib/reconstruction/reconstruction_factory.f90
  lib/finite_volume_scheme/fv_scheme.f90
  lib/integrand.f90
  lib/integrators/strategy.f90
  lib/integrators/surrogate.f90
  lib/integrators/time_integrator_factory.f90
  lib/integrators/heun_2nd_order.f90
  lib/integrators/ssp_runge_kutta_3rd_order.f90
  lib/integrators/ssp_runge_kutta_2nd_order.f90
  lib/integrators/ralston_2nd_order.f90)

# target_compile_options(cato PUBLIC ${Coarray_COMPILE_OPTIONS})
target_link_libraries(
  cato
  PUBLIC eos
         globals
         geometry
         utilities
         io
         finterp
         functional # functional # ${Coarray_LIBRARIES}
         ${HDF5_LIBRARIES}
         ${HDF5_HL_LIBRARIES})

add_executable(cato.x main.f90)
target_link_libraries(
  cato.x
  cato
  functional
  finterp
  ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES}
  ${Coarray_LIBRARIES})

# include_directories(${HDF5_INCLUDE_DIRS})

# add_library( cato lib/globals.f90 lib/math_constants.f90
# third_party/ini_parser/cfgio_mod.f90
# third_party/ini_parser/string_conv_mod.f90
# third_party/finterp/linear_interpolation_module.f90
# third_party/functional/mod_functional.f90
# third_party/functional/mod_interfaces.f90 lib/io/mod_hdf5.f90
# lib/io/mod_input.f90 lib/equation_of_state.f90 lib/limiters/slope_limiter.f90
# lib/domain/grid_type.f90 lib/domain/regular_2d_grid.f90
# lib/domain/grid_factory.f90 lib/geometry/vector_2d.f90 lib/geometry/vector.f90
# lib/geometry/mach_cone_geometry.f90 lib/geometry/quad_cell.f90
# lib/io/contour_writer.f90 lib/boundary_conditions/boundary_conditions.f90
# lib/boundary_conditions/periodic_bc.f90 lib/boundary_conditions/bc_factory.f90
# lib/domain/fluxes.f90 lib/evolution_operator/abstract_evo_operator.f90
# lib/evolution_operator/local_evo_operator.f90
# lib/reconstruction/abstract_reconstruction.f90
# lib/reconstruction/first_order_reconstruction.f90
# lib/reconstruction/second_order_reconstruction.f90
# lib/reconstruction/reconstruction_factory.f90
# lib/finite_volume_scheme/fv_scheme.f90 lib/finite_volume_scheme/cato.f90
# lib/integrand.f90 lib/integrators/strategy.f90 lib/integrators/surrogate.f90
# lib/integrators/time_integrator_factory.f90
# lib/integrators/runge_kutta_2nd_order.f90)

# set_property(TARGET cato PROPERTY INTERPROCEDURAL_OPTIMIZATION FALSE)
# set_target_properties(cato PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)

# target_link_libraries(cato ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

# add_executable(cato.x main.f90) set_property(TARGET cato.x PROPERTY
# INTERPROCEDURAL_OPTIMIZATION FALSE) set_target_properties(cato.x PROPERTIES
# INTERPROCEDURAL_OPTIMIZATION FALSE) # message(STATUS
# "INTERPROCEDURAL_OPTIMIZATION: " # ${INTERPROCEDURAL_OPTIMIZATION})
# target_link_libraries(cato.x cato)
