# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
from spack import *

class Rayleigh(MakefilePackage):
    __doc__ = 'Rayleigh is a 3-D convection code designed for the study of dynamo behavior in spherical geometry.'
    homepage = 'https://github.com/geodynamics/Rayleigh'
    url = 'https://github.com/geodynamics/Rayleigh/archive/refs/tags/v1.0.1.tar.gz'
    git = 'https://github.com/geodynamics/Rayleigh.git'
    maintainers = [
     'tukss']
    version('main', branch='main')
    version('1.1.0', sha256='93fbbdbde6088807638e4dcbd4d622203fd4753c1831bab2cb8eaeca5cba45c3')
    version('1.0.1', sha256='9c9e3b0b180f32a889f158e2ea2967f4ac2bb2124f5d264f230efb8c8f19ea36')
    version('1.0.0', sha256='4f2e8249256adff8c4b0bc377ceacf8a6441dcee54b7d36c784f05a454dc6dcf')
    version('0.9.1', sha256='ab96445fc61822fe2d2cba8729a85b36de6b541febf5759de6d614847844573f')
    version('0.9.0', sha256='63a80d1619cb639f3cb01ab82a441b77d736eee94469c47c50ab740fa81c08f4')
    depends_on('mpi')
    depends_on('fftw-api@3')
    depends_on('lapack')

    def edit(self, spec, prefix):
        if spec.satisfies('^cray-mpich'):
            # The Cray wrapper takes care of linking MPI correctly for all compilers.
            env['FC'] = 'ftn'
        else:
            env['FC'] = spec['mpi'].mpifc
        args = ['--prefix={}'.format(prefix)]

        if spec.satisfies('^mkl'):
            args.append('--with-mkl={}'.format(spec['mkl'].headers.directories[0][:-len('/include')]))
        else:
            if not spec.satisfies('^cray-fftw'):
                args.append('--with-fftw={}'.format(spec['fftw'].prefix))

            if spec.satisfies('^openblas'):
                args.append('--openblas')
            elif not spec.satisfies('^cray-libsci'):
                args.append('--with-lapack={}'.format(spec['lapack'].prefix))
                args.append('--with-blas={}'.format(spec['blas'].prefix))
            else:
                # Cray options are handled through modules and the compiler wrapper.
                args.append('--LIBFLAGS=')
                args.append('--INCLUDE=')

        args.append('--FFLAGS_DBG=-O0 -g')
        args.append('--FFLAGS_OPT=-O3 {}'.format(spec.target.optimization_flags(spec.compiler.name, str(spec.compiler.version))))

        # /dev/null as input to prevent interactive questions in configure
        configure(*args, **{'input': os.devnull})
