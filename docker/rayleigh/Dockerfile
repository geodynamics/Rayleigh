FROM geodynamics/rayleigh-buildenv-jammy:latest

# expose the port for Jupyter
EXPOSE 8888

# define Rayleigh directory
ENV RAYLEIGH_DIR /opt/Rayleigh

RUN git clone 'https://github.com/geodynamics/Rayleigh.git' ${RAYLEIGH_DIR}

WORKDIR ${RAYLEIGH_DIR}

RUN ./configure \
    -debian-mkl \
  && make \
  && make install \
  && make clean

# Add Rayleigh to the path
ENV PATH="${RAYLEIGH_DIR}/bin:${PATH}"
ENV PYTHONPATH="${RAYLEIGH_DIR}/post_processing${PYTHONPATH:+:$PYTHONPATH}"

# Allow running OpenMPI as root and avoid warning messages
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ENV OMPI_MCA_btl_vader_single_copy_mechanism=none
