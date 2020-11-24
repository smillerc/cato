
# The simple.h5 grid is used by multiple unit tests
set(IC_FILE "simple.h5")

file(COPY ${IC_FILE} DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

set(caf_test_file "${CMAKE_BINARY_DIR}/generated/caf_testing.f90")
set(test_fypp "${CMAKE_CURRENT_SOURCE_DIR}/caf_testing.fypp")
add_custom_command(OUTPUT "${caf_test_file}" COMMAND fypp "${test_fypp}" "${caf_test_file}"
                   MAIN_DEPENDENCY "${test_fypp}" VERBATIM WORKING_DIRECTORY ".")
                   
add_library(caf_testing "${caf_test_file}")

add_subdirectory(test_field)
# add_subdirectory(test_limiters)
add_subdirectory(test_bc)
add_subdirectory(test_geometry)

# add_subdirectory(test_misc) add_subdirectory(test_geometry) add_subdirectory(test_mach_cones)
# add_subdirectory(test_limiters) add_subdirectory(test_fv_schemes)
# add_subdirectory(test_reconstruction) add_subdirectory(test_interpolation)
# add_subdirectory(test_evo_operators)
