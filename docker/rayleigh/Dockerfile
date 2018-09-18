FROM geodynamics/rayleigh-buildenv-bionic:latest

RUN useradd \
  --create-home \
  rayleigh_user

USER rayleigh_user

WORKDIR /home/rayleigh_user

RUN git clone 'https://github.com/geodynamics/Rayleigh.git'

WORKDIR /home/rayleigh_user/Rayleigh

RUN ./configure \
    --with-blas='/usr' \
    --with-fftw='/usr' \
    --with-lapack='/usr' \
  && make \
  && make install \
  && make clean

ENV RAYLEIGH_DIR /home/rayleigh_user/Rayleigh

ENV PATH="${RAYLEIGH_DIR}/bin:${PATH}"
