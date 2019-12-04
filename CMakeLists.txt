cmake_minimum_required(VERSION 3.12)

project(
  fvleg_2d
  VERSION 0.0.1
  DESCRIPTION "A modern fortran implementation of the FVLEG scheme"
  HOMEPAGE_URL "https://github.com/smillerc/fvleg_2d"
  LANGUAGES Fortran C)

set(CMAKE_VERBOSE_MAKEFILE "Off")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_policy(SET CMP0074 NEW)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# set output paths for modules, archives, and executables
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Packages
find_package(HDF5 COMPONENTS Fortran_HL REQUIRED)

find_package(PFUNIT REQUIRED)
enable_testing()

message(STATUS "PFUNIT_VERSION: " ${PFUNIT_VERSION})
message(STATUS "PFUNIT_LIBRARIES: " ${PFUNIT_LIBRARIES})
message(STATUS "PFUNIT_INCLUDE_DIRS: " ${PFUNIT_INCLUDE_DIRS})


# find_package(OpenMP)
# find_package(OpenCL REQUIRED)

# if build type not specified, default to release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE}")

# compiler flags for gfortran
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  find_package(Coarray REQUIRED)
  list(APPEND CMAKE_Fortran_FLAGS ${Coarray_COMPILE_OPTIONS})

  list(APPEND CMAKE_Fortran_FLAGS "-cpp -std=f2018 -ffree-line-length-none -fcoarray=single")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g\
 -Wall -Wextra \
 -Wconversion \
 -fimplicit-none -fbacktrace \
 -fcheck=all \
 -ffpe-trap=zero,overflow,underflow -finit-real=nan")

  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -floop-parallelize-all -march=native ")
endif()

# compiler flags for ifort
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  list(APPEND CMAKE_Fortran_FLAGS "-fpp -coarray=single")
  set(CMAKE_Fortran_FLAGS_DEBUG " -O0 -g -traceback")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O2")
endif()

message(STATUS "CMAKE_REQUIRED_FLAGS" ${CMAKE_REQUIRED_FLAGS})
message(STATUS "CMAKE_REQUIRED_LIBRARIES" ${CMAKE_REQUIRED_LIBRARIES})
message(STATUS "Coarray_COMPILE_OPTIONS:" ${Coarray_COMPILE_OPTIONS})
message(STATUS "Compile Flags: " ${CMAKE_Fortran_FLAGS})
message(STATUS "Compile Flags Debug: " ${CMAKE_Fortran_FLAGS_DEBUG})
message(STATUS "Compile Flags Release: " ${CMAKE_Fortran_FLAGS_RELEASE})


# add_library(clfortran src/lib/opencl/mod_clfortran.f90)
# target_link_libraries(clfortran OpenCL::OpenCL)

# add_subdirectory(src/third_party)
add_subdirectory(src)
add_subdirectory(tests/unit)
