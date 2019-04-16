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
                sh '''sed -i '/[(79),(227),(239)]/s/${MPIEXEC}/${MPIEXEC} --allow-run-as-root/' examples/CMakeLists.txt'''
              }
            }

            stage('Build (Xenial)') {
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

            stage('Test (Xenial)') {
              steps {
                sh '''
                  cd build
                  ctest --no-compress-output -T Test
                '''
              }
              post {
                always {
                  xunit testTimeMargin: '3000',
                    thresholdMode: 1,
                    thresholds: [failed(), skipped()],
                    tools: [CTest(pattern: 'build/Testing/**/*.xml')]
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
                sh '''sed -i '/[(79),(227),(239)]/s/${MPIEXEC}/${MPIEXEC} --allow-run-as-root/' examples/CMakeLists.txt'''
              }
            }

            stage('Build (Bionic)') {
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

            stage('Test (Bionic)') {
              steps {
                sh '''
                  cd build
                  ctest --no-compress-output -T Test
                '''
              }
              post {
                always {
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
