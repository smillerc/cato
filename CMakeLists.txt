# MIT License

# Copyright (c) 2019 Sam Miller

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
# associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
# NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
# OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
cmake_minimum_required(VERSION 3.8)

project(
  CATO
  VERSION 2.0.0
  DESCRIPTION "A modern fortran implementation of the FVLEG scheme"
  HOMEPAGE_URL "https://github.com/smillerc/cato"
  LANGUAGES Fortran C)

  set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0069 NEW) # Link-time optimization
# include(CheckIPOSupported) check_ipo_supported()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# set output paths for modules, archives, and executables
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(Python)

option(USE_OPENMP_THREADS "Enable OpenMP parallelization" OFF)
option(USE_OPENMP_SIMD "Enable OpenMP SIMD vectorization" OFF)
if(USE_OPENMP_THREADS OR USE_OPENMP_SIMD)
  find_package(OpenMP REQUIRED Fortran)
endif()

# Procedures can't be pure if there is an omp parallel block, but they
# can be if it's only SIMD, so we add a definition here to turn pure
# on for certain functions if threading is off
if(USE_OPENMP_THREADS)
  add_compile_definitions(__OPENMP_THREADS__)
endif()

option(USE_OPENCL "Enable OpenCL to run portions on GPUs" OFF)
if(USE_OPENCL)
  find_package(OpenCL)
endif()

option(ENABLE_COARRAY "Build with coarrays enabled" OFF)
option(ENABLE_COARRAY_SINGLE
       "Build with coarrays but only for a single image (useful for debugging)" OFF)
option(BUILD_FOR_SHARED_MEMORY "Build for shared memory only (single node)" OFF)

option(N_IMAGES "Number of coarray images to use" 0)

option(BUILD_FOR_DISTRIBUTED_MEMORY "Build for distributed memory (multiple nodes)" OFF)
option(OUTPUT_OPTIMIZATION_REPORTS "Turn on optimization reports" OFF)
option(USE_ASAN "Enable GCC's address sanitizer to check for leaks" OFF)
option(USE_TSAN "Enable GCC's thread sanitizer to check for race conditions" OFF)
option(ENABLE_TESTING "Enable unit testing" OFF)
option(ENABLE_PROFILING "Enable profile flags" OFF)

# Determine the number of cores (useful for single machine, shared-memory build)
include(ProcessorCount)
processorcount(N_CORES)

if(ENABLE_COARRAY OR ENABLE_COARRAY_SINGLE)
  find_package(Coarray REQUIRED)
endif()

find_package(HDF5 REQUIRED Fortran HL)

# Set compiler flags
include(SetFortranFlags)

if(APPLE)
  set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
  set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
  set(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
  set(CMAKE_C_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
  set(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
  set(CMAKE_Fortran_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
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
  include(AddCoarrayCTest)
  # find_package(PFUNIT REQUIRED) message(STATUS "PFUNIT_VERSION: " ${PFUNIT_VERSION})
  # message(STATUS "PFUNIT_LIBRARIES: " ${PFUNIT_LIBRARIES}) message(STATUS "PFUNIT_INCLUDE_DIRS: "
  # ${PFUNIT_INCLUDE_DIRS}) 
  enable_testing()
  add_subdirectory(tests/unit)
endif()

message("")
message("===================================================")
message("         ${PROJECT_NAME} ${PROJECT_VERSION} build configuration")
message("===================================================")
message("USE_OPENMP_THREADS: ${USE_OPENMP_THREADS}")
message("USE_OPENMP_SIMD: ${USE_OPENMP_SIMD}")
if(OPENMP_FOUND)
  message("OpenMP")
  message("  Version   : ${OpenMP_Fortran_VERSION}")
  message("  SIMD      : ${USE_OPENMP_SIMD}")
  message("  Threading : ${USE_OPENMP_THREADS}")

  message("")
endif()

message("USE_OPENCL: ${USE_OPENCL}")
if(OPENCL_FOUND)
  message("OpenCL")
  message("  Version      : ${OpenCL_VERSION_STRING}")
  message("  Include dirs : ${OpenCL_INCLUDE_DIRS}")
  message("  Libraries    : ${OpenCL_LIBRARIES}")
  message("")
endif()

message("ENABLE_COARRAY: ${ENABLE_COARRAY}")
message("ENABLE_COARRAY_SINGLE: ${ENABLE_COARRAY_SINGLE}")
message("BUILD_FOR_SHARED_MEMORY: ${BUILD_FOR_SHARED_MEMORY}")
message("BUILD_FOR_DISTRIBUTED_MEMORY: ${BUILD_FOR_DISTRIBUTED_MEMORY}")

if(N_IMAGES GREATER 0)
  message("N_COARRAY_IMAGES: ${N_IMAGES}")
endif()

message("OUTPUT_OPTIMIZATION_REPORTS: ${OUTPUT_OPTIMIZATION_REPORTS}")
message("USE_ASAN: ${USE_ASAN}")
message("USE_TSAN: ${USE_TSAN}")
message("ENABLE_TESTING: ${ENABLE_TESTING}")
message("ENABLE_PROFILING: ${ENABLE_PROFILING}")

if(HDF5_FOUND)
  message("HDF5")
  message("  Version       : ${HDF5_VERSION}")
  message("  Definitions   : ${HDF5_DEFINITIONS}")
  message("  Parallel      : ${HDF5_IS_PARALLEL}")
  message("  Includes dirs : ${HDF5_INCLUDE_DIRS}")
  message("  Libraries     : ${HDF5_LIBRARIES}")
  message("")
endif(HDF5_FOUND)

message("Compiler Flags: ${CMAKE_Fortran_FLAGS}") 
message("")
