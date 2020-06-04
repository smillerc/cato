#!/bin/bash
FC=mpiifort
CC=gcc
cd /tmp && git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git && \
      cd pFUnit && \
      git checkout tags/v4.1.7 -b latest && \
      mkdir build && cd build && \
      FC=${FC} CC=${CC} cmake .. \
      -DSKIP_FHAMCREST=YES \
      -DSKIP_ROBUST=YES \
      -DCMAKE_INSTALL_PREFIX=/software/pfunit_intel && \
      make -j && make install
