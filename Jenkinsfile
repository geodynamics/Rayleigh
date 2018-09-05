#!groovy

pipeline {
  agent {
    kubernetes {
      //cloud 'kubernetes'
      label 'mypod'
      containerTemplate {
        name 'rayleigh'
        image 'gassmoeller/rayleigh:base'
        ttyEnabled true
        command 'cat'
      }
    }
  }

  options {
    timeout(time: 2, unit: 'HOURS')
  }

  stages {
    stage ("Check Permissions") {
      when {
        allOf {
          not {branch 'master'}
          not {changeRequest authorEmail: "feathern@colorado.edu"}
          not {changeRequest authorEmail: "rene.gassmoeller@mailbox.org"}
        }
      }
      steps {
        container('rayleigh') {
          sh '''
            wget -q -O - https://api.github.com/repos/geodynamics/Rayleigh/issues/${CHANGE_ID}/labels | grep 'ready to test' || \
            { echo "This commit will only be tested when it has the label 'ready to test'"; exit 1; }
          '''
        }
      }
    }

    stage('Build') {
      options {
        timeout(time: 15, unit: 'MINUTES')
      }
      steps {
        container('rayleigh') {
          sh '''
            ./configure --with-blas=/usr --with-fftw=/usr --with-lapack=/usr
          '''
          sh '''
            make
            make install
          '''
        }
      }
    }

    stage('Test') {
      options {
        timeout(time: 90, unit: 'MINUTES')
      }
      steps {
        container('rayleigh') {
          sh '''
            # This export avoids a warning about
            # a discovered, but unconnected infiniband network.
            export OMPI_MCA_btl=self,tcp
            cp input_examples/benchmark_diagnostics_input main_input
            # This model expects 4 MPI processes, but MPI does not work
            # inside the container at the moment, so instead run in serial
            # also reduce runtime of the model for fast testing
            # mpirun -np 4 bin/rayleigh.dbg
            sed -i -e 's/nprow = 2/nprow = 1/' \
                   -e 's/npcol = 2/npcol = 1/' \
                   -e 's/max_iterations = 40000/max_iterations = 400/' \
                   main_input
            ./bin/rayleigh.dbg
          '''
        }
      }
    }
  }
}
