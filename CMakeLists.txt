cmake_minimum_required(VERSION 3.8)

project(
  CATO
  DESCRIPTION "A modern fortran implementation of the FVLEG scheme"
  HOMEPAGE_URL "https://github.com/smillerc/cato"
  LANGUAGES Fortran C)

set(CATO_VERSION_MAJOR 1)
set(CATO_VERSION_MINOR 0)
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

# Packages
find_package(Coarray REQUIRED)
set(CMAKE_REQUIRED_FLAGS ${Coarray_COMPILE_OPTIONS})
set(CMAKE_REQUIRED_LIBRARIES ${Coarray_LIBRARY})

find_package(
  HDF5
  COMPONENTS Fortran_HL
  REQUIRED)

# Determine the host architecture, e.g. Intel Kabylake, and optimize for it
include(OptimizeForArchitecture)
optimizeforarchitecture()

# Default to a parallel build
if(NOT SERIAL_BUILD)
  set(SERIAL_BUILD Off)

  # Set for coarray compile type
  if(DISTRIBUTED_MEMORY)
    set(DISTRIBUTED_MEMORY On)
    set(SHARED_MEMORY Off)
  else()
    set(DISTRIBUTED_MEMORY Off)
    set(SHARED_MEMORY On)
  endif()

else()
  set(NUM_IMAGES 1)
  set(DISTRIBUTED_MEMORY Off)
  set(SHARED_MEMORY On)
endif()

# Set compiler flags
include(SetFortranFlags)

# Profile the code?
if(NOT ENABLE_PROFILING)
  set(ENABLE_PROFILING Off)
endif()

# Enable ctest unit testing
if(NOT ENABLE_TESTING)
  set(ENABLE_TESTING Off)
else()
  if(NOT SERIAL_BUILD)
    set(NUM_IMAGES 4)
  endif()
endif()

message(STATUS "Coarray_COMPILE_OPTIONS:" ${Coarray_COMPILE_OPTIONS})
message(STATUS "Coarray_NUMPROC_FLAG:" ${Coarray_NUMPROC_FLAG})
message(STATUS "Coarray_LIBRARIES:" ${Coarray_LIBRARIES})

if(CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  message(STATUS "Compile Flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE} ")
elseif(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  message(STATUS "Compile Flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
endif()

# Host information There is much more info available with this query, but for
# now only get the host name
cmake_host_system_information(RESULT HOST_NAME QUERY HOSTNAME)

# Insert the git version information into the version.h header file See
# https://goo.gl/697j8v (short stackoverflow link)
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_local_changes(GIT_LOCAL_CHANGES)

# This takes variables like HOST_NAME, CMAKE_Fortran_COMPILER_VERSION, etc, and
# includes it in the draco build
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/include/version.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/generated/version.h" @ONLY)

# Include the newly generated version.h file
include_directories(${CMAKE_BINARY_DIR}/generated)

add_subdirectory(src)

if(ENABLE_TESTING)
  find_package(PFUNIT REQUIRED)
  enable_testing()
  add_subdirectory(tests/unit)
endif()
