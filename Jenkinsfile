#!groovy

pipeline {
  agent {
    dockerfile {
      dir 'docker/rayleigh-buildenv-bionic'
    }
  }

  options {
    timeout(time: 2, unit: 'HOURS')
  }

  stages {
    stage('Build documentation') {
      options {
        timeout(time: 15, unit: 'MINUTES')
      }
      steps {
        // First make sure notebooks do not contain output
        sh 'make clear_ipynb && git diff --exit-code --name-only'

        // Now build the new documentation
        sh '''
          cd doc
          make html
          make latexpdf
        '''
      }
    }

    stage('Build') {
      options {
        timeout(time: 15, unit: 'MINUTES')
      }
      steps {
        sh '''
          ./configure \
            --with-blas='/usr' \
            --with-fftw='/usr' \
            --with-lapack='/usr'
        '''

        sh 'make'
        sh 'make install'
      }
    }

    stage('Test') {
      options {
        timeout(time: 90, unit: 'MINUTES')
      }
      steps {
        // Benchmark regression test
        sh '''
          cd tests/c2001_case0

          # This export avoids a warning about
          # a discovered, but unconnected infiniband network.
          mpirun -np 4 ../../bin/rayleigh.dbg

          cd ..
          git diff > changes.diff
        '''

          archiveArtifacts artifacts: 'tests/changes.diff', fingerprint: true
          sh 'git diff --exit-code --name-only'

        // Generic input test
        sh './tests/generic_input/run_test.sh'
      }
    }
  }

  post { always { cleanWs() } }
}
