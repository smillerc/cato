add_pfunit_ctest(vector_tests TEST_SOURCES test_vector.pf LINK_LIBRARIES
                 geometry)

add_pfunit_ctest(mach_cone_tests TEST_SOURCES test_mach_cone_geometry.pf
                 LINK_LIBRARIES geometry)

add_pfunit_ctest(quad_element_tests TEST_SOURCES test_quad_cell.pf
                 LINK_LIBRARIES geometry)

file(COPY simple.h5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_pfunit_ctest(grid_tests TEST_SOURCES test_grid.pf LINK_LIBRARIES geometry)
