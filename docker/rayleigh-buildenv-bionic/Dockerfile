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
    python3-pip \
    sphinx-common \
    texlive-base \
    texlive-generic-extra \
    texlive-latex-base \
    texlive-latex-recommended \
    texlive-latex-extra \
    lmodern \
    makedepf90 \
    latexmk \
    pandoc \
    python-numpy \
    python3-numpy \
    python-scipy \
    python3-scipy \
    python-funcsigs \
    python3-funcsigs \
    python-vtk6

RUN pip3 install sphinxcontrib-bibtex nbstripout nbsphinx jupyter-client ipykernel

# Export compilers
ENV CC mpicc
ENV CXX mpicxx
ENV FC mpifort
