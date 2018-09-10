#!groovy

pipeline {
  agent {
    kubernetes {
      //cloud 'kubernetes'
      label 'mypod'
      containerTemplate {
        name 'buildenv-xenial'
        image 'geodynamics/virtualquake-buildenv-xenial:latest'
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
        container('buildenv-xenial') {
          sh '''sed -i '/[(79),(227),(239)]/s/${MPIEXEC}/${MPIEXEC} --allow-run-as-root/' examples/CMakeLists.txt'''
        }
      }
    }

    stage('Build') {
      steps {
        container('buildenv-xenial') {
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
        container('buildenv-xenial') {
          sh '''
            cd build
            ctest --no-compress-output -T Test
          '''
        }
      }
      post {
        always {
          container('buildenv-xenial') {
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
