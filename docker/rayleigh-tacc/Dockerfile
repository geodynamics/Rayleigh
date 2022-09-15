FROM tacc/tacc-ubuntu18-impi19.0.7-common:latest

LABEL maintainer <rene.gassmoeller@mailbox.org>

RUN apt-get update && apt-get upgrade -y && \
    apt-get install --yes \
    git \
    libblas-dev \
    libfftw3-dev \
    liblapack-dev


# Build aspect, replace git checkout command to create image for release
RUN cd /opt && git clone 'https://github.com/geodynamics/Rayleigh.git'
COPY configure /opt/Rayleigh/

ENV FC mpifc
ENV CC mpicc

RUN cd /opt/Rayleigh && ./configure \
    --with-blas='/usr' \
    --with-fftw='/usr' \
    --with-lapack='/usr' \
  && make -j \
  && make install \
  && make clean

ENV RAYLEIGH_DIR /home/rayleigh_user/Rayleigh

ENV PATH="${RAYLEIGH_DIR}/bin:${PATH}"
