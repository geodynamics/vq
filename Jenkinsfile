#!groovy

pipeline {
  agent none

  options {
    timeout(time: 4, unit: 'HOURS')
  }

  stages {
    stage('Multiple Build Environments') {
      parallel {
        stage('Xenial') {
          agent {
            kubernetes {
              label 'mypod'
              yaml """
apiVersion: v1
kind: Pod
spec:
  containers:
    - name: 'buildenv-xenial'
      image: 'geodynamics/virtualquake-buildenv-xenial:latest'
      command:
      - 'cat'
      tty: true
"""
            }
          }
          stages {
            stage('Prepare Environment (Xenial)') {
              steps {
                container('buildenv-xenial') {
                  sh '''sed -i '/[(79),(227),(239)]/s/${MPIEXEC}/${MPIEXEC} --allow-run-as-root/' examples/CMakeLists.txt'''
                }
              }
            }

            stage('Build (Xenial)') {
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

            stage('Test (Xenial)') {
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

        stage('Bionic') {
          agent {
            kubernetes {
              label 'mypod'
              yaml """
apiVersion: v1
kind: Pod
spec:
  containers:
    - name: 'buildenv-bionic'
      image: 'geodynamics/virtualquake-buildenv-bionic:latest'
      command:
      - 'cat'
      tty: true
"""
            }
          }
          stages {
            stage('Prepare Environment (Bionic)') {
              steps {
                container('buildenv-bionic') {
                  sh '''sed -i '/[(79),(227),(239)]/s/${MPIEXEC}/${MPIEXEC} --allow-run-as-root/' examples/CMakeLists.txt'''
                }
              }
            }

            stage('Build (Bionic)') {
              steps {
                container('buildenv-bionic') {
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

            stage('Test (Bionic)') {
              steps {
                container('buildenv-bionic') {
                  sh '''
                    cd build
                    ctest --no-compress-output -T Test
                  '''
                }
              }
              post {
                always {
                  container('buildenv-bionic') {
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
      }
    }
  }
}
