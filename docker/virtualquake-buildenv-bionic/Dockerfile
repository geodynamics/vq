FROM ubuntu:bionic

RUN apt update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt install --yes \
  build-essential \
  cmake \
  doxygen \
  git \
  libgeographic-dev \
  libgomp1 \
  libopenmpi-dev \
  openmpi-bin \
  openmpi-common \
  libhdf5-openmpi-dev \
  libpython2.7-dev \
  ninja-build \
  python2.7 \
  python-h5py \
  python-matplotlib \
  python-mpltoolkits.basemap \
  swig
