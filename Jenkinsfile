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
              alwaysPull true
            }
          }

          stages {
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

          post { always { cleanWs() } }
        }

        stage('Bionic') {
          agent {
            docker {
              image 'geodynamics/virtualquake-buildenv-bionic:latest'
              alwaysPull true
            }
          }

          stages {
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

          post { always { cleanWs() } }
        }
      }
    }
  }
}
