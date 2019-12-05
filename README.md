# fvleg_2d

...a better name is coming soon



## Build/Install
Requirements:
- CMake
- gfortran 8+ or Intel Fortran 2018+
- OpenCoarrays ([https://github.com/sourceryinstitute/OpenCoarrays](https://github.com/sourceryinstitute/OpenCoarrays))
- HDF5 (for I/O)
- pFUnit (for unit testing) ([https://github.com/Goddard-Fortran-Ecosystem/pFUnit](https://github.com/Goddard-Fortran-Ecosystem/pFUnit))

Sample install script
```bash
mkdir build && cd build
CC=gcc-9 FC=gfortran-9 cmake .. -DCMAKE_BUILD_TYPE="Debug"
```

## Physics
This code solves the Euler fluid equations using the finite volume local evolution Galerkin method (FVLEG). The papers that describe this are listed [here](./papers/Readme.md)

### Influences and Sources of Inspiration

- Scientific Software Design - The Object-Oriented Way ([DOI: https://doi.org/10.1017/CBO9780511977381](https://doi.org/10.1017/CBO9780511977381))
-  Modern Fortran - Building Efficient Parallel Applications ([https://www.manning.com/books/modern-fortran](https://www.manning.com/books/modern-fortran))
-  Modern Fortran Explained - Incorporating Fortran 2018 ([Amazon](https://www.amazon.com/Modern-Fortran-Explained-Incorporating-Mathematics/dp/0198811896))
