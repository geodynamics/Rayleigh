FROM ubuntu:bionic

RUN apt update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt install --yes \
    g++ \
    gfortran \
    git \
    libblas-dev \
    libfftw3-dev \
    liblapack-dev \
    libopenmpi-dev \
    make \
    nano \
    wget

# Export compilers
ENV CC mpicc
ENV CXX mpicxx
ENV FC mpifort
