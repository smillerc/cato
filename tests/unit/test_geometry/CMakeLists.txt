add_pfunit_ctest (mach_cone_tests
  TEST_SOURCES test_mach_cone_geometry.pf
  LINK_LIBRARIES geometry
  )

add_pfunit_ctest (mesh_tests
  TEST_SOURCES test_quad_cell.pf
  LINK_LIBRARIES geometry
  )

add_pfunit_ctest (grid_tests
  TEST_SOURCES test_grid.pf
  LINK_LIBRARIES geometry
  )
