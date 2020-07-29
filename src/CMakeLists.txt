add_library(globals lib/hermetic.f90 lib/units.f90 lib/globals.f90 lib/non_dimensional.f90
                    lib/math_constants.f90 lib/floating_point_utilities.f90)

add_library(ini_parser third_party/ini_parser/cfgio_mod.f90
                       third_party/ini_parser/string_conv_mod.f90)

add_library(finterp third_party/finterp/linear_interpolation_module.f90)

add_library(functional third_party/functional/mod_functional.f90
                       third_party/functional/mod_interfaces.f90)

add_library(io lib/io/mod_hdf5.f90 lib/io/mod_input.f90)

add_library(eos lib/equation_of_state.f90)

include_directories(${HDF5_INCLUDE_DIRS})
target_include_directories(io PUBLIC ${HDF5_INCLUDE_DIRS})

target_link_libraries(io PUBLIC globals ini_parser ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
target_link_libraries(eos PUBLIC io)

add_library(
  utilities lib/limiters/slope_limiter.f90 lib/limiters/flux_limiter.f90
  # lib/limiters/gnoffo_aux_limiter.f90 lib/limiters/post_limiter.f90
)
target_link_libraries(utilities PUBLIC globals)

add_library(
  geometry
  lib/grid/grid_type.f90
  lib/grid/regular_2d_grid.f90
  lib/grid/grid_factory.f90
  lib/geometry/geometry.f90
  lib/geometry/intersections.f90
  lib/geometry/vector_2d.f90
  lib/geometry/vector.f90
  lib/geometry/quad_cell.f90)

target_link_libraries(geometry PUBLIC globals io eos functional)

### Pre-process: .fypp -> .f90 via Fypp

# Create a list of the files to be preprocessed
set(fyppFiles
    lib/limiters/e_mlp_distinguisher.fypp
    lib/flux_solvers/m_ausmpw_plus_solver.fypp
    lib/flux_solvers/ausmpw_plus_solver.fypp
    lib/flux_solvers/slau_solver.fypp)

# Pre-process
foreach(infileName IN LISTS fyppFiles)

    # Generate output file name
    string(REGEX REPLACE ".fypp\$" ".f90" outfileName "${infileName}")
    # Generate input file name
    set(infile "${CMAKE_CURRENT_SOURCE_DIR}/${infileName}")

    # Create the full path for the new file
    get_filename_component(f90_file ${outfileName} NAME)
    set(outfile "${CMAKE_BINARY_DIR}/generated/${f90_file}")

    # Custom command to do the processing
    add_custom_command(
        OUTPUT "${outfile}"
        COMMAND fypp "${infile}" "${outfile}"
        MAIN_DEPENDENCY "${infile}"
        VERBATIM)

    # Finally add output file to a list
    set(outFiles ${outFiles} "${outfile}")

endforeach(infileName)

add_library(
  cato
  ${outFiles}
  lib/mach_cone/mach_cone_utilities.f90
  lib/mach_cone/mach_cone_collection.f90
  lib/mod_timing.f90
  lib/io/contour_writer.f90
  lib/boundary_conditions/boundary_conditions.f90
  lib/boundary_conditions/zero_gradient_bc.f90
  lib/boundary_conditions/inlet_outlet.f90
  lib/boundary_conditions/periodic_bc.f90
  lib/boundary_conditions/pressure_input_bc.f90
  lib/boundary_conditions/symmetry_bc.f90
  lib/boundary_conditions/bc_factory.f90
  lib/source_terms/source_factory.f90
  lib/source_terms/mod_source.f90
  lib/source_terms/energy_source.f90
  lib/fluid/fluxes.f90
  lib/fluid/flux_array.f90
  lib/fluid/mod_fluid.f90
  lib/flux_solvers/flux_solver.f90
  lib/flux_solvers/fvleg_solver.f90
  lib/flux_solvers/ausm_plus_solver.f90
  lib/evolution_operator/local_evo_operator.f90
  lib/spatial_reconstruction/abstract_reconstruction.f90
  lib/spatial_reconstruction/gradients.f90
  lib/spatial_reconstruction/piecewise_constant.f90
  lib/spatial_reconstruction/piecewise_linear.f90
  lib/spatial_reconstruction/edge_interpolator.f90
  lib/spatial_reconstruction/edge_interpolator_factory.f90
  lib/spatial_reconstruction/tvd_2nd_order.f90
  lib/spatial_reconstruction/tvd_3rd_order.f90
  lib/spatial_reconstruction/tvd_5th_order.f90
  lib/spatial_reconstruction/mlp_baseline.f90
  lib/spatial_reconstruction/mlp_3rd_order.f90
  lib/spatial_reconstruction/mlp_5th_order.f90
  lib/spatial_reconstruction/reconstruction_factory.f90
  lib/master_puppeteer/master_puppeteer.f90)

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
         # OpenMP::OpenMP_Fortran
         ${HDF5_LIBRARIES}
         ${HDF5_HL_LIBRARIES})

add_executable(cato.x main.f90)
target_link_libraries(
  cato.x cato functional finterp ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES}
  # ${Coarray_LIBRARIES}
)

# include_directories(${HDF5_INCLUDE_DIRS})

# add_library( cato lib/globals.f90 lib/math_constants.f90 third_party/ini_parser/cfgio_mod.f90
# third_party/ini_parser/string_conv_mod.f90 third_party/finterp/linear_interpolation_module.f90
# third_party/functional/mod_functional.f90 third_party/functional/mod_interfaces.f90
# lib/io/mod_hdf5.f90 lib/io/mod_input.f90 lib/equation_of_state.f90 lib/limiters/slope_limiter.f90
# lib/domain/grid_type.f90 lib/domain/regular_2d_grid.f90 lib/domain/grid_factory.f90
# lib/geometry/vector_2d.f90 lib/geometry/vector.f90 lib/geometry/mach_cone_geometry.f90
# lib/geometry/quad_cell.f90 lib/io/contour_writer.f90
# lib/boundary_conditions/boundary_conditions.f90 lib/boundary_conditions/periodic_bc.f90
# lib/boundary_conditions/bc_factory.f90 lib/domain/fluxes.f90
# lib/evolution_operator/abstract_evo_operator.f90 lib/evolution_operator/local_evo_operator.f90
# lib/reconstruction/abstract_reconstruction.f90 lib/reconstruction/first_order_reconstruction.f90
# lib/reconstruction/second_order_reconstruction.f90 lib/reconstruction/reconstruction_factory.f90
# lib/master_puppeteer/master.f90 lib/master_puppeteer/cato.f90 lib/integrand.f90
# lib/integrators/strategy.f90 lib/integrators/surrogate.f90
# lib/integrators/time_integrator_factory.f90 lib/integrators/runge_kutta_2nd_order.f90)

# set_property(TARGET cato PROPERTY INTERPROCEDURAL_OPTIMIZATION FALSE) set_target_properties(cato
# PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE)

# target_link_libraries(cato ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

# add_executable(cato.x main.f90) set_property(TARGET cato.x PROPERTY INTERPROCEDURAL_OPTIMIZATION
# FALSE) set_target_properties(cato.x PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE) #
# message(STATUS "INTERPROCEDURAL_OPTIMIZATION: " # ${INTERPROCEDURAL_OPTIMIZATION})
# target_link_libraries(cato.x cato)
