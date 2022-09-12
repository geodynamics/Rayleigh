FROM ubuntu:focal

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
    wget \
    python3-pip \
    texlive-base \
    texlive-latex-base \
    texlive-latex-recommended \
    texlive-latex-extra \
    lmodern \
    makedepf90 \
    latexmk \
    pandoc \
    python3-numpy \
    python3-scipy \
    python3-funcsigs

RUN pip3 install sphinxcontrib-bibtex==1.0.0 nbstripout nbsphinx jupyter-client ipykernel recommonmark vtk

# Export compilers
ENV CC mpicc
ENV CXX mpicxx
ENV FC mpifort
