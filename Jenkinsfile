#!groovy

pipeline {
  agent {
    kubernetes {
      //cloud 'kubernetes'
      label 'mypod'
      containerTemplate {
        name 'ubuntu1604'
        image 'ubuntu:16.04'
        ttyEnabled true
        command 'cat'
      }
    }
  }

  environment {
    DEBIAN_FRONTEND = 'noninteractive'
    DEBCONF_NONINTERACTIVE_SEEN = 'true'
  }

  options {
    timeout(time: 4, unit: 'HOURS') 
  }

  container('ubuntu1604') {
    stages {
      stage ("Prepare Environment") {
        steps {
          sh 'apt update'
          sh '''
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
          '''
        }
      }

      stage('Build') {
        steps {
          sh 'mkdir build'
          sh '''
            cd build
            cmake \
              -G "Ninja" \
              -D GeographicLib_INCLUDE_DIRS="/usr/include/GeographicLib" \
              -D GeographicLib_LIBRARY_DIRS="/usr/lib/x86_64-linux-gnu/" \
              ..
           '''
           sh '''
             cd build
             ninja
           '''
        }
      }

      stage('Test') {
        steps {
          sh '''
            cd build
            ctest
          '''
        }
      }
    } 
  }
}
