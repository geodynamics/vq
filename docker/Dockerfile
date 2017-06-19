FROM ubuntu:16.04

RUN apt-get update && apt-get install -y \
  build-essential \
  cmake \
  doxygen \
  git \
  libgomp1 \
  libopenmpi-dev \
  openmpi-bin \
  openmpi-common \
  libhdf5-openmpi-10 \
  libhdf5-openmpi-dev \
  libpython2.7-dev \
  libpython2.7 \
  python2.7 \
  python-h5py \
  python-matplotlib \
  python-mpltoolkits.basemap \
  swig \
  vim \
  nano

RUN groupadd \
    --system \
    --gid=2000 \
    virtualquake \
  && useradd \
    --create-home \
    --gid virtualquake \
    --no-log-init \
    --system \
    --uid=2000 \
    virtualquake

ADD GeographicLib-1.48.tar.gz /home/virtualquake/
ADD vq-3.1.0.tar.gz /home/virtualquake/

RUN chown -R virtualquake:virtualquake \
    /home/virtualquake/GeographicLib-1.48 \
  && chown -R virtualquake:virtualquake \
    /home/virtualquake/vq-3.1.0

USER virtualquake

WORKDIR /home/virtualquake/GeographicLib-1.48

RUN mkdir BUILD \
  && cd BUILD \
  && cmake .. \
  && cmake --build . -- all

USER root

RUN cd BUILD \
  && cmake --build . -- install

RUN apt-get install -y \
    libgeographic14

USER virtualquake

WORKDIR /home/virtualquake/vq-3.1.0

RUN mkdir -p build \
  && cd build \
  && cmake .. \
  && cmake --build . -- all 

USER root

RUN cd build \
  && cmake --build . -- install

USER virtualquake

WORKDIR /home/virtualquake

RUN mkdir external_volume

WORKDIR /home/virtualquake/external_volume

CMD [ "/bin/bash" ]
