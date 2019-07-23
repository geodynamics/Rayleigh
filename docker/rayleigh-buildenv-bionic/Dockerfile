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
    libmpich-dev \
    make \
    nano \
    wget \
    sphinx-common \
    texlive-base \
    texlive-generic-extra \
    texlive-latex-base \
    texlive-latex-recommended \
    texlive-latex-extra \
    lmodern \
    python3-pip \
    latexmk

RUN pip3 install sphinxcontrib-bibtex

# Export compilers
ENV CC mpicc
ENV CXX mpicxx
ENV FC mpifort
