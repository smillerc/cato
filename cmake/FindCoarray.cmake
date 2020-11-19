# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying file Copyright.txt or
# https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindCoarray
----------

Finds compiler flags or library necessary to support Fortran 2008/2018 coarrays.

This packages primary purposes are:

* for compilers natively supporting Fortran coarrays without needing compiler options, simply indicating Coarray_FOUND  (example: Cray)
* for compilers with built-in Fortran coarray support, enable compiler option (example: Intel Fortran)
* for compilers needing a library such as OpenCoarrays, presenting library (example: GNU)


Result Variables
^^^^^^^^^^^^^^^^

``Coarray_FOUND``
  indicates coarray support found (whether built-in or library)

``Coarray_LIBRARIES``
  coarray library path
``Coarray_COMPILE_OPTIONS``
  coarray compiler options
``Coarray_EXECUTABLE``
  coarray executable e.g. ``cafrun``
``Coarray_MAX_NUMPROCS``
  maximum number of parallel processes
``Coarray_NUMPROC_FLAG``
  use for executing in parallel: ${Coarray_EXECUTABLE} ${Coarray_NUMPROC_FLAG} ${Coarray_MAX_NUMPROCS} ${CMAKE_CURRENT_BINARY_DIR}/myprogram

Cache Variables
^^^^^^^^^^^^^^^^

The following cache variables may also be set:

``Coarray_LIBRARY``
  The coarray libraries, if needed and found
#]=======================================================================]

cmake_policy(VERSION 3.3)

set(options_coarray Intel) # flags needed
set(opencoarray_supported GNU) # future: Flang, etc.

unset(Coarray_COMPILE_OPTIONS)
unset(Coarray_LIBRARY)
unset(Coarray_REQUIRED_VARS)

if(CMAKE_Fortran_COMPILER_ID IN_LIST options_coarray)

  if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)

    if(WIN32)
      set(Coarray_COMPILE_OPTIONS /Qcoarray:shared)
      list(APPEND Coarray_REQUIRED_VARS ${Coarray_COMPILE_OPTIONS})
    elseif(UNIX AND NOT APPLE)
      if(ENABLE_COARRAY_SINGLE)
        set(IFORT_COARRAY "-coarray=single")
        set(Coarray_MAX_NUMPROCS 1)
      elseif(BUILD_FOR_DISTRIBUTED_MEMORY)
        set(IFORT_COARRAY "-coarray=distributed -coarray-config-file=caf_config.txt")
      elseif(BUILD_FOR_SHARED_MEMORY)
        set(IFORT_COARRAY "-coarray=shared")
      else()
        set(IFORT_COARRAY "-coarray")
        set(Coarray_MAX_NUMPROCS ${N_CORES})
      endif()

      if(N_IMAGES GREATER 0 AND NOT ENABLE_COARRAY_SINGLE)
        set(IFORT_COARRAY "${IFORT_COARRAY} -coarray-num-images=${N_IMAGES}")
      endif()

      set(Coarray_COMPILE_OPTIONS ${IFORT_COARRAY})
      set(Coarray_LIBRARY ${IFORT_COARRAY}) # ifort requires it at build AND link
      list(APPEND Coarray_REQUIRED_VARS ${Coarray_LIBRARY})
    endif()
  endif()

elseif(CMAKE_Fortran_COMPILER_ID IN_LIST opencoarray_supported)

  find_package(OpenCoarrays)

  if(OpenCoarrays_FOUND)
    set(Coarray_LIBRARY OpenCoarrays::caf_mpi)

    set(Coarray_EXECUTABLE cafrun)

    include(ProcessorCount)
    processorcount(Nproc)
    set(Coarray_MAX_NUMPROCS ${Nproc})
    set(Coarray_NUMPROC_FLAG -np)

    list(APPEND Coarray_REQUIRED_VARS ${Coarray_LIBRARY})
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    set(Coarray_COMPILE_OPTIONS -fcoarray=single)
    list(APPEND Coarray_REQUIRED_VARS ${Coarray_COMPILE_OPTIONS})
  endif()

endif()

set(CMAKE_REQUIRED_FLAGS ${Coarray_COMPILE_OPTIONS})
set(CMAKE_REQUIRED_LIBRARIES ${Coarray_LIBRARY})
include(CheckFortranSourceCompiles)
check_fortran_source_compiles("real :: x[*]; end" f08coarray SRC_EXT f90)
unset(CMAKE_REQUIRED_FLAGS)
unset(CMAKE_REQUIRED_LIBRARIES)

list(APPEND Coarray_REQUIRED_VARS ${f08coarray})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Coarray REQUIRED_VARS Coarray_REQUIRED_VARS)

if(Coarray_FOUND)
  set(Coarray_LIBRARIES ${Coarray_LIBRARY})
endif()

message(STATUS "Coarray_COMPILE_OPTIONS:" ${Coarray_COMPILE_OPTIONS})
message(STATUS "Coarray_LIBRARIES:" ${Coarray_LIBRARIES})
message(STATUS "Coarray_MAX_NUMPROCS:" ${Coarray_MAX_NUMPROCS})

mark_as_advanced(Coarray_LIBRARY Coarray_REQUIRED_VARS)
