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
        sh 'cp input_examples/benchmark_diagnostics_input main_input'

        // This model expects 4 MPI processes, but MPI does not work
        // inside the container at the moment, so instead run in serial
        // also reduce runtime of the model for fast testing
        sh '''
          sed \
            --in-place \
            -e 's/max_iterations = 40000/max_iterations = 400/' \
            main_input
        '''

        sh '''
          # This export avoids a warning about
          # a discovered, but unconnected infiniband network.
          mpirun -np 4 ./bin/rayleigh.dbg
        '''
      }
    }
  }

  post { always { cleanWs() } }
}
