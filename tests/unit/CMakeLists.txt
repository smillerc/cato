
# The simple.h5 grid is used by multiple unit tests
set(IC_FILE "simple.h5")
# # add_custom_command(
# #    OUTPUT ${IC_FILE}
# #    COMMAND "python make_simple_grid.py"
# #    DEPENDS make_simple_grid.py)

# add_custom_command(OUTPUT "${IC_FILE}"
#    COMMAND "${PYTHON_EXECUTABLE}" "make_simple_grid.py" --argument
#    DEPENDS generate_test_grid)

# add_custom_target(generate_test_grid ALL DEPENDS "${IC_FILE}")

file(COPY ${IC_FILE} DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

add_library(caf_testing caf_testing.f90)

add_subdirectory(test_field)
# add_subdirectory(test_limiters)
add_subdirectory(test_bc)

# add_subdirectory(test_misc) add_subdirectory(test_geometry) add_subdirectory(test_mach_cones)
# add_subdirectory(test_limiters) add_subdirectory(test_fv_schemes)
# add_subdirectory(test_reconstruction) add_subdirectory(test_interpolation)
# add_subdirectory(test_evo_operators)
