FROM tacc/tacc-ubuntu18-impi19.0.7-common:latest

LABEL maintainer <rene.gassmoeller@mailbox.org>

RUN apt-get update && apt-get upgrade -y && \
    apt-get install --yes \
    git wget add-apt-key

RUN cd /tmp && \
    # now get the key: \
    wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB && \
    # now install that key \
    apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB && \
    # now remove the public key file exit the root shell \
    rm GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB

RUN echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list

RUN apt-get update

RUN apt-get install --yes \
    intel-mkl-2020.0-088

# Build aspect, replace git checkout command to create image for release
RUN cd /opt && git clone 'https://github.com/geodynamics/Rayleigh.git'

ENV FC mpifc
ENV CC mpicc
ENV MKLROOT /opt/intel/compilers_and_libraries_2020.0.166/linux/mkl
ENV LD_LIBRARY_PATH="${MKLROOT}/lib/intel64:${LD_LIBRARY_PATH}"

RUN cd /opt/Rayleigh \
   && ./configure --FFLAGS_OPT="-march=skylake-avx512 -O3 -ffixed-line-length-132" \
    -mkl \
  && make -j \
  && make install \
  && make clean

ENV RAYLEIGH_DIR /opt/Rayleigh

ENV PATH="${RAYLEIGH_DIR}/bin:${PATH}"
