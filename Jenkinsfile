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

  stages {
    stage ("Prepare Environment") {
      steps {
        container('ubuntu1604') {
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
          sh '''sed -i '/[(79),(227),(239)]/s/${MPIEXEC}/${MPIEXEC} --allow-run-as-root/' examples/CMakeLists.txt'''
        }
      }
    }

    stage('Build') {
      steps {
        container('ubuntu1604') {
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
    }

    stage('Test') {
      steps {
        container('ubuntu1604') {
          sh '''
            cd build
            ctest --no-compress-output -T Test
          '''
        }
      }
      post {
        always {
          container('ubuntu1604') {
            xunit testTimeMargin: '3000',
              thresholdMode: 1,
              thresholds: [failed(), skipped()],
              tools: [CTest(pattern: 'build/Testing/**/*.xml')]
          }
        }
      }
    }
  }
}
