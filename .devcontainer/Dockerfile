FROM fedora:latest

ARG IFORTSERIAL
COPY intel.cfg /tmp/intel-base.cfg

RUN dnf clean all
RUN dnf install @development-tools -y

RUN dnf install -y \
    kernel-devel kernel-headers valgrind gdb \
    libtool file gfortran g++ gcc m4 \
    openmpi-devel cmake wget zsh curl git \
    cgnslib-devel zlib vim ack

RUN cd /tmp && \
    (cat intel-base.cfg | sed s/INTELSERIAL/$IFORTSERIAL/ > intel-config.cfg) && \
    wget -O intel.tgz http://registrationcenter-download.intel.com/akdlm/irc_nas/tec/16225/parallel_studio_xe_2020_cluster_edition_online.tgz && \
    tar -xvzf intel.tgz && \
    cd /tmp/parallel_studio_xe_* && \
    mkdir -p /software && \
    bash ./install.sh --silent=/tmp/intel-config.cfg && \
    cd / && \
    rm -rf /tmp/* \
    /software/intel/ism \
    /software/intel/documentation* \
    /software/intel/ide* \
    /software/intel/samples*

ENV PATH /software/intel/bin:$PATH
ENV PATH /software/intel/compilers_and_libraries/linux/mpi/intel64/bin:$PATH
ENV LD_LIBRARY_PATH /software/intel/compilers_and_libraries/linux/mpi/intel64/lib:$LD_LIBRARY_PATH
ENV LD_LIBRARY_PATH /software/intel/lib:$LD_LIBRARY_PATH

ENV FC=ifort

# Add pFunit testing (intel)
RUN cd /tmp && git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git && \
    cd pFUnit && mkdir build && cd build && \
    source /software/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64 && \
    FC=mpiifort CC=gcc cmake .. -DCMAKE_INSTALL_PREFIX=/software/pfunit_intel && \
    make && make install ; exit 0

# ENV PFUNIT_DIR /software/pfunit_intel/PFUNIT-*
# ENV FARGPARSE_DIR /software/pfunit_intel/FARGPARSE-*
# ENV GFTL_DIR /software/pfunit_intel/GFTL-*
# ENV GFTL_SHARED_DIR /software/pfunit_intel/GFTL_SHARED-*

RUN rm -rf /tmp/*

# Add opencoarrays
RUN cd /tmp && git clone https://github.com/sourceryinstitute/OpenCoarrays.git && \
    cd OpenCoarrays && mkdir build && cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/software/opencoarrays && \
    make && make install

ENV PATH /software/opencoarrays/bin:${PATH}
ENV LD_LIBRARY_PATH /software/opencoarrays/lib:$LD_LIBRARY_PATH

# Add pFunit testing (gfortran)
RUN cd /tmp && git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git && \
    cd pFUnit && mkdir build && cd build && \
    source /software/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64 && \
    FC=gfortran CC=gcc cmake .. -DCMAKE_INSTALL_PREFIX=/software/pfunit_gnu && \
    make && make install ; exit 0

# ENV PFUNIT_DIR /software/pfunit_gnu/PFUNIT-*
# ENV FARGPARSE_DIR /software/pfunit_gnu/FARGPARSE-*
# ENV GFTL_DIR /software/pfunit_gnu/GFTL-*
# ENV GFTL_SHARED_DIR /software/pfunit_gnu/GFTL_SHARED-*

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

RUN cd /software && git clone https://github.com/spack/spack.git
ENV PATH /software/spack/bin:$PATH
RUN spack compiler find
RUN spack install hdf5 %intel +fortran -mpi +hl

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

# terminal colors with xterm
ENV TERM xterm

WORKDIR ${USER_HOME}

# Set the prompt look
RUN echo "parse_git_branch() {" >> ~/.bashrc
RUN echo "  git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/ (\1)/'" >> ~/.bashrc
RUN echo "}" >> ~/.bashrc

# RUN echo "parse_git_branch() { git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/ (\1)/' }" >> ~/.bashrc
RUN echo 'export PS1="\n\[$(tput sgr0)\]\[\033[38;5;75m\]\h\[$(tput sgr0)\]\[\033[38;5;15m\]: [\w] \$(parse_git_branch) \n\[$(tput sgr0)\]\[\033[38;5;10m\]\\$>\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]"' >> ~/.bashrc

RUN echo "source /software/intel/compilers_and_libraries_2020/linux/bin/compilervars.sh intel64" >> ~/.bashrc
RUN echo "source /software/intel/compilers_and_libraries_2020/linux/mkl/bin/mklvars.sh intel64" >> ~/.bashrc
RUN echo "source /software/spack/share/spack/setup-env.sh" >> ~/.bashrc
CMD ["/bin/bash"]
