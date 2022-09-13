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
    python-pip \
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
    python-scipy \
    python-funcsigs

RUN pip3 install sphinx==4.0 sphinxcontrib-bibtex==1.0.0 nbstripout traitlets==4.3.3 nbsphinx==0.8.5 jupyter-client ipykernel recommonmark sphinx-book-theme
RUN pip2 install vtk

# Export compilers
ENV CC mpicc
ENV CXX mpicxx
ENV FC mpifort
