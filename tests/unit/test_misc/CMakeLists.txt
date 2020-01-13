add_pfunit_ctest(flux_tests TEST_SOURCES test_fluxes.pf LINK_LIBRARIES cato)
# set_tests_properties(flux_tests PROPERTIES RUN_SERIAL TRUE) set(NAME
# flux_tests PROPERTIES COMMAND ${Coarray_EXECUTABLE} ${Coarray_NUMPROC_FLAG} 1
# $<TARGET_FILE:flux_tests>)
