# Determine and set the Fortran compiler flags we want

# Make sure the build type is uppercase
string(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
message(STATUS "Build type: ${BT}")

if(BT STREQUAL "RELEASE")
  set(CMAKE_BUILD_TYPE
      RELEASE
      CACHE STRING "Choose the type of build, options are DEBUG or RELEASE" FORCE)
elseif(BT STREQUAL "DEBUG")
  set(CMAKE_BUILD_TYPE
      DEBUG
      CACHE STRING "Choose the type of build, options are DEBUG or RELEASE" FORCE)
elseif(NOT BT)
  set(CMAKE_BUILD_TYPE
      RELEASE
      CACHE STRING "Choose the type of build, options are DEBUG or RELEASE" FORCE)
  message(STATUS "CMAKE_BUILD_TYPE not given, defaulting to RELEASE")
else()
  message(FATAL_ERROR "CMAKE_BUILD_TYPE not valid, choices are DEBUG or RELEASE")
endif(BT STREQUAL "RELEASE")

# Determine the host architecture, e.g. Intel Kabylake, and optimize for it
include(OptimizeForArchitecture)
optimizeforarchitecture()

# Set __ALIGNBYTES__ for use in OpenMP/Compiler SIMD alignment statements, e.g. __ALIGNBYTES__ = 32
add_compile_definitions(__ALIGNBYTES__=${MEMORY_ALIGN_BYTES})

# gfortran
if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)

  set(CMAKE_Fortran_FLAGS "-cpp -std=f2018 -ffree-line-length-none")

  if(USE_OPENMP_THREADS)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
  endif()

  if(USE_OPENMP_SIMD)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp-simd")
  endif()

  if(ENABLE_COARRAY_SINGLE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcoarray=lib")
  elseif(ENABLE_COARRAY)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcoarray=single")
  endif()

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-O0 -g \
 -Wall -Wextra -Wpedantic -Wconversion \
 -fimplicit-none -fbacktrace \
 -fcheck=all -ffpe-trap=zero,overflow,invalid,underflow -finit-real=nan")

  set(CMAKE_Fortran_FLAGS_RELEASE
      "-O3 -ftree-vectorize -funroll-loops -finline-functions -march=native -mtune=native")

  if(USE_ASAN)
    set(CMAKE_Fortran_FLAGS
        "${CMAKE_Fortran_FLAGS} -g -fsanitize=leak -fsanitize=address -fno-omit-frame-pointer -fopt-info-all"
    )
  endif()

  if(USE_TSAN)
    set(CMAKE_Fortran_FLAGS
        "${CMAKE_Fortran_FLAGS} -fsanitize=thread -fno-omit-frame-pointer -fopt-info-all")
  endif()

  if(OUTPUT_OPTIMIZATION_REPORTS)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopt-info-all=optim.txt")
    message(STATUS "GFortran optimization reports will be in src/optim.txt")
  endif()

endif()

# ifort
if(CMAKE_Fortran_COMPILER_ID STREQUAL Intel)

  add_compile_definitions(__SIMD_ALIGN_OMP__)
  set(IFORT_FLAGS
      "-fpp -inline-max-size=300 -align array${MEMORY_ALIGN_BYTES}byte -fp-model source -assume contiguous_assumed_shape -diag-disable 5268 -diag-disable 7025 -diag-disable 8770 -diag-disable 6477 ${Coarray_COMPILE_OPTIONS}"
  )
  # set(IFORT_FLAGS "-fpp -fp-model precise -fp-model except -diag-disable 5268 -diag-disable 8770
  # ${Coarray_COMPILE_OPTIONS}" )

  if(USE_OPENMP)
    set(IFORT_FLAGS "${IFORT_FLAGS} ${OpenMP_Fortran_FLAGS}")
  endif()

  # Fortran 2018 standards check based on the version
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 19.0.3)
    set(CMAKE_Fortran_FLAGS "-stand f15 ${IFORT_FLAGS}")
  else()
    set(CMAKE_Fortran_FLAGS "-stand f18 ${IFORT_FLAGS}")
  endif()

  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -warn all -debug all -traceback -fpe-all=0 -check bounds")
  set(CMAKE_Fortran_FLAGS_RELEASE " -O3 -xHost -mtune=${TARGET_ARCHITECTURE} -qopt-prefetch")

  if(OUTPUT_OPTIMIZATION_REPORTS)
    set(CMAKE_Fortran_FLAGS
        "${CMAKE_Fortran_FLAGS}  -g -qopt-report-phase=all -qopt-report=5"
    )
    message(
      STATUS
        "Intel Fortran optimization reports will be in the src/CMakeFiles directories with *.optrpt extensions"
    )
  endif()

endif()
