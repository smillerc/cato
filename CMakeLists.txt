cmake_minimum_required(VERSION 3.8)

project(
  fvleg_2d
  DESCRIPTION "A modern fortran implementation of the FVLEG scheme"
  HOMEPAGE_URL "https://github.com/smillerc/fvleg_2d"
  LANGUAGES Fortran C)

set(fvleg_2d_VERSION_MAJOR 0)
set(fvleg_2d_VERSION_MINOR 1)
set(fvleg_2d_PATCH_VERSION 0)

set(CMAKE_VERBOSE_MAKEFILE "Off")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_policy(SET CMP0074 NEW)
# cmake_policy(SET CMP0069 NEW) # Link-time optimization
# include(CheckIPOSupported) check_ipo_supported()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# set output paths for modules, archives, and executables
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Packages
find_package(
  HDF5
  COMPONENTS Fortran_HL
  REQUIRED)

find_package(PFUNIT REQUIRED)
enable_testing()

message(STATUS "PFUNIT_VERSION: " ${PFUNIT_VERSION})
message(STATUS "PFUNIT_LIBRARIES: " ${PFUNIT_LIBRARIES})
message(STATUS "PFUNIT_INCLUDE_DIRS: " ${PFUNIT_INCLUDE_DIRS})

# find_package(OpenMP) # find_package(OpenCL REQUIRED)

# if build type not specified, default to release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}")

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
  message(STATUS "Debug mode, so IPO is set to OFF")
endif(CMAKE_BUILD_TYPE STREQUAL Debug)
message(STATUS "CMAKE_INTERPROCEDURAL_OPTIMIZATION: "
               ${CMAKE_INTERPROCEDURAL_OPTIMIZATION})

# compiler flags for gfortran
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  find_package(Coarray REQUIRED)
  # list(APPEND CMAKE_Fortran_FLAGS ${Coarray_COMPILE_OPTIONS})

  list(APPEND CMAKE_Fortran_FLAGS "-cpp -std=f2018 -ffree-line-length-none")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -g -fcoarray=single\
 -Wall -Wextra -Wpedantic\
 -Wconversion \
 -fimplicit-none -fbacktrace \
 -fcheck=all") # -pg -ffpe-trap=zero,overflow,underflow -finit-real=nan")

  set(CMAKE_Fortran_FLAGS_RELEASE
      "-O3 -floop-parallelize-all -ftree-parallelize-loops=6 -march=native -fcoarray=single"
  )
endif()

# compiler flags for ifort
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  # list( APPEND CMAKE_Fortran_FLAGS "-fpp -std18 -fp-model precise -fp-model
  # except -check bounds -diag-disable 5268 -diag-disable 8770 -check
  # noarg_temp_created" )
  list(
    APPEND
      CMAKE_Fortran_FLAGS
      "-fpp -std18 -check bounds -diag-disable 5268 -diag-disable 8770 -check noarg_temp_created"
  )

  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -debug all -traceback -fpe-all=0"
  )# -coarray=singl
  # e
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2 -parallel") # -coarray
endif()

message(STATUS "CMAKE_REQUIRED_FLAGS" ${CMAKE_REQUIRED_FLAGS})
message(STATUS "CMAKE_REQUIRED_LIBRARIES" ${CMAKE_REQUIRED_LIBRARIES})
message(STATUS "Coarray_COMPILE_OPTIONS:" ${Coarray_COMPILE_OPTIONS})
message(STATUS "Compile Flags: " ${CMAKE_Fortran_FLAGS})
message(STATUS "Compile Flags Debug: " ${CMAKE_Fortran_FLAGS_DEBUG})
message(STATUS "Compile Flags Release: " ${CMAKE_Fortran_FLAGS_RELEASE})

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

# add_library(clfortran src/lib/opencl/mod_clfortran.f90)
# target_link_libraries(clfortran OpenCL::OpenCL)

add_subdirectory(src)
add_subdirectory(tests/unit)
