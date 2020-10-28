add_caf_test(
  test_field
  N_IMAGES
  4
  SOURCES
  test_field.f90
  LINK_LIBRARIES
  caf_testing
  cato)

# add_caf_test(test_data_layout N_IMAGES 1 SOURCES test_data_layout.f90 LINK_LIBRARIES
# OpenMP::OpenMP_Fortran)
