FROM geodynamics/rayleigh-buildenv-bionic

RUN apt update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt install --yes \
    sudo

COPY /entrypoint.sh /entrypoint.sh

ENTRYPOINT /entrypoint.sh

# Export compilers
ENV CC mpicc
ENV CXX mpicxx
ENV FC mpifort
