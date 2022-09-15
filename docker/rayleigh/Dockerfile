FROM geodynamics/rayleigh-buildenv-jammy:latest

RUN useradd \
  --create-home \
  rayleigh_user

USER rayleigh_user

WORKDIR /home/rayleigh_user

RUN git clone 'https://github.com/geodynamics/Rayleigh.git'

WORKDIR /home/rayleigh_user/Rayleigh

RUN ./configure \
    -debian-mkl \
  && make \
  && make install \
  && make clean

ENV RAYLEIGH_DIR /home/rayleigh_user/Rayleigh

ENV PATH="${RAYLEIGH_DIR}/bin:${PATH}"
