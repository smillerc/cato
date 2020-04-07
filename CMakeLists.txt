cmake_minimum_required(VERSION 3.8)

project(
  CATO
  DESCRIPTION "A modern fortran implementation of the FVLEG scheme"
  HOMEPAGE_URL "https://github.com/smillerc/cato"
  LANGUAGES Fortran C)

set(CATO_VERSION_MAJOR 1)
set(CATO_VERSION_MINOR 1)
set(CATO_PATCH_VERSION 0)

set(CMAKE_VERBOSE_MAKEFILE "Off")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_policy(SET CMP0074 NEW)
# cmake_policy(SET CMP0069 NEW) # Link-time optimization
# include(CheckIPOSupported)
# check_ipo_supported()

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

# if build type not specified, default to release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}")

# Profile the code?
if(NOT ENABLE_PROFILING)
  set(ENABLE_PROFILING Off)
endif()

# Set some debug options
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
  message(STATUS "Debug mode, so IPO is set to OFF")
else()

endif(CMAKE_BUILD_TYPE STREQUAL Debug)

include(OptimizeForArchitecture)
optimizeforarchitecture()

message(STATUS "CMAKE_Fortran_COMPILER_VERSION: "
               ${CMAKE_Fortran_COMPILER_VERSION})

# compiler flags for gfortran
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)

  list(APPEND CMAKE_Fortran_FLAGS "-cpp -std=f2018 -ffree-line-length-none")
  list(APPEND CMAKE_Fortran_FLAGS ${Coarray_COMPILE_OPTIONS})
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -g -fcoarray=lib\
 -Wall -Wextra -Wpedantic\
 -Wconversion \
 -fimplicit-none -fbacktrace \
 -fcheck=all -ffpe-trap=zero,overflow,invalid,underflow -finit-real=nan")

  set(CMAKE_Fortran_FLAGS_RELEASE
      "-O3 -floop-parallelize-all -ftree-parallelize-loops=6 -march=native -fcoarray=lib"
  )

  if(ENABLE_PROFILING)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pg -g")
  endif(ENABLE_PROFILING)

endif()

# compiler flags for ifort
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 19.0.3)
    list(
      APPEND
        CMAKE_Fortran_FLAGS
        "-fpp -stand f15 -fp-model precise -fp-model except -diag-disable 5268 -diag-disable 8770"
    )
  else()
    list(
      APPEND
        CMAKE_Fortran_FLAGS
        "-fpp -stand f18 -fp-model precise -fp-model except -diag-disable 5268 -diag-disable 8770"
    )
  endif()

  if(ENABLE_PROFILING)
    set(CMAKE_Fortran_FLAGS
        "${CMAKE_Fortran_FLAGS} -p -qopt-report-phase=all -g -qopt-report-annotate=html -qopt-report-annotate-position=both -qopt-report=5"
    )
  endif(ENABLE_PROFILING)

  # check_ipo_supported()

  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -g -warn all -debug all -traceback -fpe-all=0 -check bounds ")
  set(CMAKE_Fortran_FLAGS_RELEASE
      " -O3 -xHost -ipo -parallel -mtune=${TARGET_ARCHITECTURE}")
endif()

message(STATUS "CMAKE_REQUIRED_FLAGS" ${CMAKE_REQUIRED_FLAGS})
message(STATUS "CMAKE_REQUIRED_LIBRARIES" ${CMAKE_REQUIRED_LIBRARIES})
message(STATUS "Coarray_COMPILE_OPTIONS:" ${Coarray_COMPILE_OPTIONS})
message(STATUS "Coarray_NUMPROC_FLAG:" ${Coarray_NUMPROC_FLAG})
message(STATUS "Coarray_LIBRARIES:" ${Coarray_LIBRARIES})
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

add_subdirectory(src)

if(NOT ENABLE_TESTING)
  set(ENABLE_TESTING Off)
endif()

if(ENABLE_TESTING)
  # find_package(OpenMP)  # Sometimes pFunit need OpenMP, but iris does not
  find_package(PFUNIT REQUIRED)
  message(STATUS "PFUNIT_VERSION: " ${PFUNIT_VERSION})
  message(STATUS "PFUNIT_LIBRARIES: " ${PFUNIT_LIBRARIES})
  message(STATUS "PFUNIT_INCLUDE_DIRS: " ${PFUNIT_INCLUDE_DIRS})
  enable_testing()
  add_subdirectory(tests/unit)
endif()
