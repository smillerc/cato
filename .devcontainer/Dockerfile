# A baseline docker image for mpi dev using openmpi

FROM fedora:latest

RUN dnf groupinstall "Development Tools" -y

RUN dnf install -y \
    libtool file gfortran g++ gcc m4 gdb valgrind \
    openmpi-devel cmake wget zsh curl git \
    cgnslib-devel zlib vim ack

RUN cd /tmp \
    && wget https://github.com/CGNS/CGNS/archive/v3.3.1.tar.gz  \
    && tar -xvf v3.3.1.tar.gz \
    && cd CGNS-3.3.1 \
    && mkdir build && cd build \
    && cmake .. \
    -DCGNS_ENABLE_HDF5=ON -DCGNS_ENABLE_FORTRAN=ON -Wno-dev \
    -DHDF5_DIR=/usr -DCGNS_ENABLE_64BIT=ON -DHDF5_NEED_ZLIB=ON \
    -DCMAKE_INSTALL_PREFIX=/software/cgns \
    && make && make install \
    && cd / && rm -rf /tmp/*

ENV PATH /software/cgns/bin:${PATH}
ENV LD_LIBRARY_PATH /software/cgns/lib:$LD_LIBRARY_PATH

ENV PATH /usr/lib64/openmpi/bin:${PATH}
ENV LD_LIBRARY_PATH /usr/lib64/openmpi/lib:$LD_LIBRARY_PATH

# Add opencoarrays
RUN git clone https://github.com/sourceryinstitute/OpenCoarrays.git && \
    cd OpenCoarrays && mkdir build && cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/software/opencoarrays && \
    make && make install

ENV PATH /software/opencoarrays/bin:${PATH}
ENV LD_LIBRARY_PATH /software/opencoarrays/lib:$LD_LIBRARY_PATH

# Add pFunit testing
RUN cd /tmp && git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git && \
    cd pFUnit && mkdir build && cd build && \
    FC=mpifort CC=gcc cmake .. -DCMAKE_INSTALL_PREFIX=/software && \
    make && make install ; exit 0

ENV PFUNIT_DIR /software/PFUNIT-4.0
ENV FARGPARSE_DIR /software/FARGPARSE-0.9
ENV GFTL_DIR /software/GFTL-1.1
ENV GFTL_SHARED_DIR /software/GFTL_SHARED-1.0

RUN rm -rf /tmp/*

# Add a default non-root user to run mpi jobs
ARG USER=mpi
ENV USER ${USER}
RUN adduser ${USER} \
    && echo "${USER}   ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

ENV USER_HOME /home/${USER}
RUN chown -R ${USER}:${USER} ${USER_HOME}

# Install miniconda to /miniconda
RUN cd /tmp && \
    curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b

ENV PATH=/miniconda/bin:${PATH}

RUN pip install pre-commit fprettify fortran-language-server

# Create working directory
ARG WORKDIR=/project
ENV WORKDIR ${WORKDIR}
RUN mkdir ${WORKDIR}
RUN chown -R ${USER}:${USER} ${WORKDIR}

WORKDIR ${WORKDIR}
USER ${USER}

# Test mpi
RUN mkdir /tmp/mpi-tests
WORKDIR /tmp/mpi-tests
COPY mpi-test .
RUN sh test.sh
RUN rm -rf /tmp/mpi-test

# run the installation script
RUN wget https://github.com/robbyrussell/oh-my-zsh/raw/master/tools/install.sh -O - | zsh || true
# terminal colors with xterm
ENV TERM xterm
# set the zsh theme
ENV ZSH_THEME agnoster

WORKDIR ${USER_HOME}

CMD ["/bin/zsh"]
