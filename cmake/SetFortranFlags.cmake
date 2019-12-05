# ##############################################################################
# Determine and set the Fortran compiler flags we want
# ##############################################################################
# https://github.com/SethMMorton/cmake_fortran_template

# ##############################################################################
# Make sure that the default build type is RELEASE if not specified.
# ##############################################################################
include(${CMAKE_MODULE_PATH}/SetCompileFlag.cmake)

# Make sure the build type is uppercase
string(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

message(STATUS "Build type: ${BT}")

if(BT STREQUAL "RELEASE")
  set(CMAKE_BUILD_TYPE
      RELEASE
      CACHE STRING
            "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
            FORCE)
elseif(BT STREQUAL "DEBUG")
  set(CMAKE_BUILD_TYPE
      DEBUG
      CACHE STRING
            "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
            FORCE)
elseif(BT STREQUAL "TESTING")
  set(CMAKE_BUILD_TYPE
      TESTING
      CACHE STRING
            "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
            FORCE)
elseif(NOT BT)
  set(CMAKE_BUILD_TYPE
      RELEASE
      CACHE STRING
            "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
            FORCE)
  message(STATUS "CMAKE_BUILD_TYPE not given, defaulting to RELEASE")
else()
  message(
    FATAL_ERROR
      "CMAKE_BUILD_TYPE not valid, choices are DEBUG, RELEASE, or TESTING")
endif(BT STREQUAL "RELEASE")

# gfortran
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  set(CMAKE_Fortran_FLAGS "-cpp -fconvert=big-endian -ffree-line-length-none -fno-stack-arrays")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Wall -Wextra -Wpedantic -fimplicit-none -fbacktrace -fcheck=all -fcheck=bounds -ffpe-trap=zero,overflow,underflow -finit-real=nan ")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -floop-parallelize-all -funroll-loops -finline-functions")
endif()

# ifort
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)
  set(CMAKE_Fortran_FLAGS "-fpp -convert big_endian -heap-arrays -fp-model precise -fp-model except -check bounds -std08 -diag-disable 5268 -diag-disable 8770")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback -warn all -debug extended")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -axCORE-AVX2 -unroll -inline")
endif()

# There is some bug where -march=native doesn't work on Mac
if(APPLE)
  set(GNUNATIVE "-mtune=native")
else()
  set(GNUNATIVE "-march=native")
endif()
