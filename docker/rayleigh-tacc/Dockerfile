FROM tacc/tacc-base:ubuntu22.04-impi19.0.9-common

LABEL maintainer <rene.gassmoeller@mailbox.org>

RUN apt-get update || true && \
    DEBIAN_FRONTEND=noninteractive apt-get install --yes \
    git libmkl-dev libfftw3-dev

RUN cd /opt && git clone 'https://github.com/geodynamics/Rayleigh.git'

ENV FC mpifc
ENV CC mpicc
ENV MKLROOT /opt/intel/compilers_and_libraries_2020.0.166/linux/mkl
ENV LD_LIBRARY_PATH="${MKLROOT}/lib/intel64:${LD_LIBRARY_PATH}"

RUN cd /opt/Rayleigh && \
   ./configure --with-fftw=/usr --FFLAGS_OPT="-march=skylake-avx512 -O3 -ffixed-line-length-132" \
     \
  && make -j \
  && make install \
  && make clean

ENV RAYLEIGH_DIR /opt/Rayleigh

ENV PATH="${RAYLEIGH_DIR}/bin:${PATH}"
