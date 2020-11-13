# Since we test this unit with 1 image and 4 images, make a library that will
# be used by both
add_library(mod_test_bc mod_test_bc.f90)
target_link_libraries(mod_test_bc PRIVATE caf_testing cato ${Coarray_LIBRARIES})

# Add the single image version test
add_caf_test(test_bc_1_image N_IMAGES 1 SOURCES test_bc.f90 LINK_LIBRARIES mod_test_bc)

# Add the 4 image version test
add_caf_test(test_bc_4_image N_IMAGES 4 SOURCES test_bc.f90 LINK_LIBRARIES mod_test_bc)
