FROM ubuntu:jammy

RUN apt-get update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt-get install --yes \
    g++ \
    gfortran \
    git \
    intel-mkl \
    libopenmpi-dev \
    make \
    nano \
    wget \
    python3-pip \
    sphinx-common \
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
    python3-funcsigs \
    python3-sphinxcontrib.bibtex \
    python3-sphinx-book-theme \
    python3-recommonmark \
    python3-nbsphinx \
    python3-vtk7 \
    jupyter \
    && apt-get clean && rm -r /var/lib/apt/lists/*

RUN pip3 install nbstripout

# Export compilers
ENV CC mpicc
ENV CXX mpicxx
ENV FC mpifort
