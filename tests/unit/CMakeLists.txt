add_library(caf_testing caf_testing.f90)

add_subdirectory(test_field)
# add_subdirectory(test_limiters)
add_subdirectory(test_bc)

# add_subdirectory(test_misc) add_subdirectory(test_geometry) add_subdirectory(test_mach_cones)
# add_subdirectory(test_limiters) add_subdirectory(test_fv_schemes)
# add_subdirectory(test_reconstruction) add_subdirectory(test_interpolation)
# add_subdirectory(test_evo_operators)
