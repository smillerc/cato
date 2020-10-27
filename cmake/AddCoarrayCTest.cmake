# --------------------------------------------------------------------------------------------------
# ----------------------------
# Add a CTest unit test that uses coarrays
# --------------------------------------------------------------------------------------------------
# ----------------------------
function(add_caf_test test_target)

  # Get the function arguments
  set(oneValueArgs N_IMAGES)
  set(multiValueArgs SOURCES LINK_LIBRARIES)

  cmake_parse_arguments(TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Make the test executable and link any additional libraries
  add_executable(${test_target} ${TEST_SOURCES})
  if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
    target_compile_options(${test_target} PUBLIC -coarray-num-images=${TEST_N_IMAGES})
  endif()

  target_link_libraries(${test_target} ${TEST_LINK_LIBRARIES} ${Coarray_LIBRARIES})

  # This sets it to build_dir/bin, e.g. project/build/bin
  set(test_dir ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

  # Function to add Coarray Fortran tests.
  if(TARGET ${test_target})
    get_target_property(MIN_TEST_IMGS ${test_target} MIN_IMAGES)
  elseif(TARGET build_${test_target})
    get_target_property(MIN_TEST_IMGS build_${test_target} MIN_IMAGES)
  endif()
  if(MIN_TEST_IMGS)
    if(TEST_N_IMAGES LESS MIN_TEST_IMGS)
      message(
        FATAL_ERROR
          "Test ${test_target} requires ${MIN_TEST_IMGS} but was only given ${num_caf_images}")
    endif()
  endif()

  if(((NUM_IMAGES LESS TEST_N_IMAGES) OR (NUM_IMAGES EQUAL 0)))
    message(
      STATUS
        "Test ${test_target} is oversubscribed: ${TEST_N_IMAGES} CAF images requested with ${NUM_IMAGES} system processor available."
    )
    if(openmpi)
      if(MIN_TEST_IMGS)
        set(TEST_N_IMAGES ${MIN_TEST_IMGS})
      elseif(NUM_IMAGES LESS 2)
        set(TEST_N_IMAGES 2)
      endif()
      set(test_parameters --oversubscribe)
      message(
        STATUS
          "Open-MPI back end detected, passing --oversubscribe for oversubscribed test, ${test_target}, with ${TEST_N_IMAGES} ranks/images."
      )
    endif()
  endif()

  # Set the number of processors to use
  set(test_parameters -np ${TEST_N_IMAGES} ${test_parameters})

  # Set compiler specifics
  if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    add_test(NAME ${test_target} COMMAND "bash" cafrun ${test_parameters}
                                         "${test_dir}/${test_target}")
  elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
    add_test(NAME ${test_target} COMMAND "${test_dir}/${test_target}")
  elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Cray")
    add_test(NAME ${test_target} COMMAND aprun ${test_parameters} "${test_dir}/${test_target}")
  elseif()
    message(WARNING "Unsupported compiler: tests might not launch correctly.")
    add_test(NAME ${test_target} COMMAND "${test_dir}/${test_target}")
  endif()

  # Each test must end with a print*, "Test passed." to succeed
  set_property(TEST ${test_target} PROPERTY PASS_REGULAR_EXPRESSION "Success")
  set_property (TEST ${test_target}
  PROPERTY FAIL_REGULAR_EXPRESSION "Test failure"
  )
endfunction(add_caf_test)
