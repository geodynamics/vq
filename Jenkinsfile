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
            docker {
              image 'geodynamics/virtualquake-buildenv-xenial:latest'
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
            docker {
              image 'geodynamics/virtualquake-buildenv-bionic:latest'
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
