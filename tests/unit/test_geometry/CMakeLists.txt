add_pfunit_ctest(vector_tests TEST_SOURCES test_vector.pf LINK_LIBRARIES geometry)

add_pfunit_ctest(quad_element_tests TEST_SOURCES test_quad_cell.pf LINK_LIBRARIES geometry)

add_pfunit_ctest(intersection_tests TEST_SOURCES test_intersections.pf LINK_LIBRARIES geometry)

# The simple.h5 grid is also needed by the fvleg tests, but this will place it in the same directory
# to be reused
file(COPY simple.h5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_pfunit_ctest(grid_tests TEST_SOURCES test_grid.pf LINK_LIBRARIES geometry)
