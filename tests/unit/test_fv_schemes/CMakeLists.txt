add_pfunit_ctest(fvleg_tests TEST_SOURCES test_fvleg.pf LINK_LIBRARIES fvleg)

file(COPY ../test_geometry/simple.h5 DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
