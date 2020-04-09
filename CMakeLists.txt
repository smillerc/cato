cmake_minimum_required(VERSION 3.8)

project(
  CATO
  DESCRIPTION "A modern fortran implementation of the FVLEG scheme"
  HOMEPAGE_URL "https://github.com/smillerc/cato"
  LANGUAGES Fortran C)

set(CATO_VERSION_MAJOR 1)
set(CATO_VERSION_MINOR 1)
set(CATO_PATCH_VERSION 0)

set(CMAKE_VERBOSE_MAKEFILE Off)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED On)

cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0069 NEW) # Link-time optimization
# include(CheckIPOSupported) check_ipo_supported()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# set output paths for modules, archives, and executables
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Determine the host architecture, e.g. Intel Kabylake, and optimize for it
include(OptimizeForArchitecture)
optimizeforarchitecture()

# Determine the number of cores (useful for single machine, shared-memory build)
include(ProcessorCount)
processorcount(N_CORES)

# Enable ctest unit testing
if(NOT ENABLE_TESTING)
  set(ENABLE_TESTING Off)
else()
  if(NOT SERIAL_BUILD)
    # The unit tests for the coarray portions assume 4 cores (or images). Certain indexing and other
    # features will not be the same when the number of images is not what are expected
    if(N_CORES LESS 6)
      message(WARNING "N_CORES < 6, and the Coarray unit tests assume 4 cores " /
                      "(coarray unit tests will not be enabled)")
    else()
      set(N_CORES 6)
    endif()

  endif()
endif()

# Default to a parallel build
if(DISTRIBUTED_MEMORY AND SHARED_MEMORY) 
  message(FATAL_ERROR "Cannot enable DISTRIBUTED_MEMORY and SHARED_MEMORY at the same time")
endif()

if(NOT SERIAL_BUILD)
  set(SERIAL_BUILD Off)
  # Set for coarray compile type
  if(DISTRIBUTED_MEMORY)
    message(STATUS "Configuring for parallel distributed-memory build")
    set(DISTRIBUTED_MEMORY On)
    set(SHARED_MEMORY Off)
  else()
    message(STATUS "Configuring for parallel shared-memory build")
    set(DISTRIBUTED_MEMORY Off)
    set(SHARED_MEMORY On)
  endif()

else()
  # SERIAL_BUILD
  message(STATUS "Configuring for serial single-core build")
  set(N_CORES 1)
  set(DISTRIBUTED_MEMORY Off)
  set(SHARED_MEMORY On)
endif()

# Packages
find_package(Coarray REQUIRED)
find_package(
  HDF5
  COMPONENTS Fortran_HL
  REQUIRED)

# Set compiler flags
include(SetFortranFlags)

# Profile the code?
if(NOT ENABLE_PROFILING)
  set(ENABLE_PROFILING Off)
endif()


if(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  message(STATUS "Compile Flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE} ")
elseif(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  message(STATUS "Compile Flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
endif()

cmake_host_system_information(RESULT HOST_NAME QUERY HOSTNAME)
cmake_host_system_information(RESULT N_PHYSICAL_CORES QUERY NUMBER_OF_PHYSICAL_CORES)
cmake_host_system_information(RESULT N_LOGICAL_CORES QUERY NUMBER_OF_LOGICAL_CORES)

# Insert the git version information into the version.h header file See https://goo.gl/697j8v (short
# stackoverflow link)
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_local_changes(GIT_LOCAL_CHANGES)

# Include build info in the binary
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/include/version.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/generated/version.h" @ONLY)

# Include the newly generated version.h file
include_directories(${CMAKE_BINARY_DIR}/generated)

add_subdirectory(src)

if(ENABLE_TESTING)
  find_package(PFUNIT REQUIRED)
  message(STATUS "PFUNIT_VERSION: " ${PFUNIT_VERSION})
  message(STATUS "PFUNIT_LIBRARIES: " ${PFUNIT_LIBRARIES})
  message(STATUS "PFUNIT_INCLUDE_DIRS: " ${PFUNIT_INCLUDE_DIRS})
  enable_testing()
  add_subdirectory(tests/unit)
endif()
