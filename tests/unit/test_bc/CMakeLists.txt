
# Test the single image version
add_caf_test(test_bc_1_image N_IMAGES 1 SOURCES test_bc.f90 LINK_LIBRARIES caf_testing cato)

# # Test the 4 image version
# add_caf_test(test_bc_4_image N_IMAGES 4 SOURCES test_bc.f90 LINK_LIBRARIES caf_testing cato)
