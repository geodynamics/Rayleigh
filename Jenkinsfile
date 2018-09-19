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
        // This model expects 4 MPI processes, but MPI does not work
        // inside the container at the moment, so instead run in serial
        // also reduce runtime of the model for fast testing
        sh '''
          sed \
            --in-place \
            -e 's/max_iterations = 40000/max_iterations = 400/' \
            main_input
        '''

<<<<<<< HEAD
        sh '''
          cd tests/benchmark_diagnostics_input

          # This export avoids a warning about
          # a discovered, but unconnected infiniband network.
          mpirun -np 4 ../../bin/rayleigh.dbg
          git diff > changes.diff
        '''
=======
            # This export avoids a warning about
            # a discovered, but unconnected infiniband network.
            export OMPI_MCA_btl=self,tcp
            ../../bin/rayleigh.dbg
            cd ..
            git diff > changes.diff
          '''
>>>>>>> Simplify test case

          archiveArtifacts artifacts: 'tests/changes.diff', fingerprint: true
          sh 'git diff --exit-code --name-only'
        }
      }
    }
  }

  post { always { cleanWs() } }
}
