#!groovy

pipeline {
  agent {
    kubernetes {
      //cloud 'kubernetes'
      label 'mypod'
      containerTemplate {
        name 'ubuntu1804'
        image 'ubuntu:18.04'
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

  stages {
    stage ("Prepare Environment") {
      steps {
        container('ubuntu1804') {
          sh 'apt update'
          sh '''
            apt install --yes \
              build-essential \
              cmake \
              doxygen \
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
    }

    stage('Build') {
      steps {
        container('ubuntu1804') {
          sh 'mkdir build'
          sh 'cd build'
          sh '''
            cmake \
              -G "Ninja" \
              -D GeographicLib_INCLUDE_DIRS="/usr/include/GeographicLib" \
              -D GeographicLib_LIBRARY_DIRS="/usr/lib/x86_64-linux-gnu/" \
              ..
           '''
           sh 'ninja'
        }
      }
    }

    stage('Test') {
      steps {
        container('ubuntu1804') {
          sh 'cd build'
          sh 'ctest'
        }
      } 
    }
  }
}
