#
#  Copyright (C) 2018 by the authors of the RAYLEIGH code.
#
#  This file is part of RAYLEIGH.
#
#  RAYLEIGH is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3, or (at your option)
#  any later version.
#
#  RAYLEIGH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with RAYLEIGH; see the file LICENSE.  If not see
#  <http://www.gnu.org/licenses/>.
#

from __future__ import print_function
import numpy as np
import os
import glob

maxq = 4000

def get_lut(quantities):
    """return the lookup table based on the quantity codes"""
    nq = len(quantities)
    lut = np.zeros(maxq) + maxq
    for i,q in enumerate(quantities):
        if ((0 <= q) and ( q <= maxq-1)): # quantity must be in [0, maxq-1]
            lut[q] = i
    return lut.astype('int32')

class Spherical_3D:
    """Rayleigh Spherical_3D Structure
    ----------------------------------
    self.basefilename                             : base filename
    self.nr                                       : number of radial points
    self.ntheta                                   : number of theta points
    self.nphi                                     : number of phi points sampled
    self.r                                        : radial coordinates
    self.theta                                    : co-latitudinal coordinates
    self.vals[0:nphi-1,0:ntheta=1,0:nr-1]         : 3-D array of values
    """
 
    def __init__(self,filename,path='Spherical_3D/'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if not Spherical_3D)
        """
        self.basefilename = os.path.split(filename)[-1].split("_")[0]

        grid_file = path+self.basefilename+"_grid"
        fd = open(grid_file,'rb')
        bs = check_endian(fd,314,'int32')
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        ntheta = swapread(fd,dtype='int32',count=1,swap=bs)
        nphi = swapread(fd,dtype='int32',count=1,swap=bs)
        assert(nphi == 2*ntheta)

        self.nr = nr
        self.ntheta = ntheta
        self.nphi = nphi

        rs = swapread(fd,dtype='float64',count=self.nr,swap=bs)
        thetas = swapread(fd,dtype='float64',count=self.ntheta,swap=bs)
        fd.close()

        fd = open(path+filename, 'rb')
        self.r = rs
        self.theta = thetas
        self.vals = np.reshape(swapread(fd,dtype='float64',count=nphi*ntheta*nr,swap=bs),(nphi,ntheta,nr), order = 'F')
        fd.close()


class Spherical_3D_multi:
    """Rayleigh Spherical_3D Structure
    ----------------------------------
    self.basefilename                             : base filename
    self.nr                                       : number of radial points
    self.ntheta                                   : number of theta points
    self.nphi                                     : number of phi points sampled
    self.rs                                       : radial coordinates
    self.thetas                                   : co-latitudinal coordinates
    self.vals[qindex][0:nphi-1,0:ntheta=1,0:nr-1] : dictionary of values, indexed by qindex
    """
 
    def __init__(self,filename,path='Spherical_3D/'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if not Spherical_3D)
        """
        self.basefilename = os.path.split(filename)[-1].split("_")[0]

        grid_file = path+self.basefilename+"_grid"
        fd = open(grid_file,'rb')
        bs = check_endian(fd,314,'int32')
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        ntheta = swapread(fd,dtype='int32',count=1,swap=bs)
        nphi = swapread(fd,dtype='int32',count=1,swap=bs)
        assert(nphi == 2*ntheta)

        self.nr = nr
        self.ntheta = ntheta
        self.nphi = nphi

        rs = swapread(fd,dtype='float64',count=self.nr,swap=bs)
        thetas = swapread(fd,dtype='float64',count=self.ntheta,swap=bs)
        fd.close()

        self.rs = rs
        self.thetas = thetas

        var_files = glob.glob1(path, self.basefilename+"_*")
        self.vals = {}
        for var_file in var_files:
          index = var_file.split("_")[-1]
          if index == "grid": continue
          fd = open(path+var_file, 'rb')
          self.vals[index] = swapread(fd,dtype='float64',count=nphi*ntheta*nr,swap=bs)



class RayleighTiming:

    def __init__(self,filename,byteswap=True):
        """filename  : The reference state file to read.
        """       
        fd = open(filename,'rb')
        # We read an integer to assess which endian the file was written in...
        #bs = check_endian(fd,314,'int32')
        bs = byteswap
        self.ncol    = swapread(fd,dtype='int32',count=1,swap=bs)
        self.nrow    = swapread(fd,dtype='int32',count=1,swap=bs)
        self.ntimers = swapread(fd,dtype='int32',count=1,swap=bs)
        self.nr      = swapread(fd,dtype='int32',count=1,swap=bs)
        self.lmax    = swapread(fd,dtype='int32',count=1,swap=bs)
        self.niter   = swapread(fd,dtype='int32',count=1,swap=bs)
        self.np = self.nrow*self.ncol
        self.col_rank = np.reshape(swapread(fd,dtype='int32',count=self.np,swap=bs),(self.np), order = 'F')
        self.row_rank = np.reshape(swapread(fd,dtype='int32',count=self.np,swap=bs),(self.np), order = 'F')
        tcount = self.np*self.ntimers
        self.times = np.reshape(swapread(fd,dtype='float64',count=tcount,swap=bs),
                        (self.ntimers,self.np), order = 'F')

        self.names = ['Main Loop', 'Legendre Transform', 'FFT',
                      'Implicit Solve', 'Row Transpose', 'Column Transpose',
                      'Hybrid Space (Return)', 'Hybrid Space (Forward)',
                      'Physical Space', 'Post Solve', 'D_by_Dphi', 'Nonlinear Terms',
                      'Sin(theta) Division', 'CFL Calculation', 'NULL', 
                      'Linear Coefficients/Implicit Matrix Computation',
                      'Run Initialization', 'Checkpointing (Read)', 'Checkpointing (Write)',
                      'Total Runtime'] 

class RayleighProfile:
    """Rayleigh Reference State Structure
    ----------------------------------
    self.nr         : number of radial points
    self.nq         : number of quantities in the 2-D structure file
    self.radius      : radial coordinates
    self.vals        : vals[0:nr-1,0:nq-1]

    """

    def __init__(self,filename='none'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
        """
        if (filename != 'none'):
            
            fd = open(filename,'rb')
            # We read an integer to assess which endian the file was written in...
            bs = check_endian(fd,314,'int32')
            
            nr = swapread(fd,dtype='int32',count=1,swap=bs)
            n2 = swapread(fd,dtype='int32',count=1,swap=bs)
            nq = n2-1
            tmp = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr,1), order = 'F')
            self.radius      = tmp[:,0]
            tmp2 = np.reshape(swapread(fd,dtype='float64',count=nq*nr,swap=bs),(nr,nq), order = 'F')
            self.nr = nr
            self.nq = nq
            self.vals = tmp2[:,:]
        else:
            self.nr = 0
            self.nq = 0
            self.vals = []
            # we initialize the object and set its attributes later

        fd.close()

class RayleighArray:
    """Rayleigh 2-D Array Structure
    ----------------------------------
    self.nx         : number of x
    self.ny         : number of y-values in the 2-D structure file
    self.vals        : vals[0:nr-1,0:nq-1]

    """

    def __init__(self,filename ='none'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
        """

        if (filename == 'none'):
            self.nx = 0
            self.ny = 0
            self.vals = []
        else:
            print('Opening: ', filename)
            fd = open(filename,'rb')
            # We read an integer to assess which endian the file was written in...
            bs = check_endian(fd,314,'int32')
            
            nx = swapread(fd,dtype='int32',count=1,swap=bs)
            ny = swapread(fd,dtype='int32',count=1,swap=bs)

            tmp2 = np.reshape(swapread(fd,dtype='float64',count=nx*ny,swap=bs),(nx,ny), order = 'F')
            self.nx = nx
            self.ny = ny
            self.vals = tmp2[:,:]

            fd.close()
    def set_vals(self,vals):
        dims = vals.shape
        self.nx = dims[0]
        self.ny = dims[1]
        self.vals = vals
    def write(self,arrfile,byteswap = False):
        fd = open(arrfile,'wb')
        dims = np.ndarray((3),dtype='int32')
        dims[0] = 314
        dims[1] = self.nx
        dims[2] = self.ny
        swapwrite(dims,fd,swap=byteswap,array=True,verbose=True)
        swapwrite(self.vals,fd,swap=byteswap,array=True,verbose=True)
        fd.close()

class ReferenceState:
    """Rayleigh Reference State Structure
    ----------------------------------
    self.n_r         : number of radial points
    self.radius      : radial coordinates
    self.density     : density
    self.dlnrho      : logarithmic derivative of density
    self.d2lnrho     : d_by_dr of dlnrho
    self.pressure    : pressure
    self.temperature : temperature
    self.dlnt        : logarithmic derivative of temperature
    self.dsdr        : entropy gradient (radial)
    self.entropy     : entropy
    self.gravity     : gravity
    self.heating     : volumetric heating (Q) (only after Jan 2019)
    """

    def __init__(self,filename='none',path='./'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
        """
        if (filename == 'none'):
            the_file = path+'reference'
        else:
            the_file = path+filename
        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        heating_written = True # First assume heating was written (as it will be from here on out
        try:
            tmp = np.reshape(swapread(fd,dtype='float64',count=11*nr,swap=bs),(nr,11), order = 'F')
        except: # Heating was not written (different-size binary 'reference')
            fd.close() # close and reopen the file to start from the beginning
            # Read in two ints first to make sure the float arrays read in start at the right place!
            fd = open(the_file,'rb')
            dummy = swapread(fd, dtype='int32',count=1,swap=bs)
            dummy = swapread(fd, dtype='int32',count=1,swap=bs)
            tmp = np.reshape(swapread(fd,dtype='float64',count=10*nr,swap=bs),(nr,10), order = 'F')
            heating_written = False
        self.nr = nr
        self.radius      = tmp[:,0]
        self.density     = tmp[:,1]
        self.dlnrho      = tmp[:,2]
        self.d2lnrho     = tmp[:,3]
        self.pressure    = tmp[:,4]
        self.temperature = tmp[:,5]
        self.dlnt        = tmp[:,6]
        self.dsdr        = tmp[:,7]
        self.entropy     = tmp[:,8]
        self.gravity     = tmp[:,9]
        if heating_written:
            self.heating     = tmp[:,10]
        self.ref = tmp
        if heating_written:
            self.names = ['radius', 'density', 'dlnrho', 'd2lnrho', 'pressure', 'temperature',
        'dlnt', 'dsdr','entropy','gravity', 'heating']
        else:
            self.names = ['radius', 'density', 'dlnrho', 'd2lnrho', 'pressure', 'temperature',
        'dlnt', 'dsdr','entropy','gravity']
        fd.close()

class TransportCoeffs:
    """Rayleigh Transport Coefficients Structure
    ----------------------------------
    self.n_r         : number of radial points
    self.radius      : radial coordinates
    self.nu          : momentum diffusivity (kinematic viscosity)
    self.dlnu        : logarithmic derivative of the viscosity
    self.kappa       : temperature diffusivity (thermometric conductivity)
    self.eta :       : magnetic diffusivity 
    self.dlneta      : logarithmic derivative of magnetic diffusivity
    """

    def __init__(self,filename='none',path='./'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
        """
        if (filename == 'none'):
            the_file = path+'transport'
        else:
            the_file = path+filename
        fd = open(the_file, 'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd, 314, 'int32')
        
        nr = swapread(fd, dtype='int32', count=1, swap=bs)
        mag_flag = swapread(fd, dtype='int32', count=1, swap=bs)
        if (mag_flag == 0):
            tmp = np.reshape(swapread(fd, dtype='float64', count=5*nr, swap=bs),(nr,5), order = 'F')
        elif (mag_flag == 1):
            tmp = np.reshape(swapread(fd, dtype='float64', count=7*nr, swap=bs),(nr,7), order = 'F')

        self.nr = nr
        self.radius      = tmp[:, 0]
        self.nu          = tmp[:, 1]
        self.dlnu        = tmp[:, 2]
        self.kappa       = tmp[:, 3]
        self.dlnkappa    = tmp[:, 4]
        if (mag_flag == 1):
            self.eta         = tmp[:, 5]
            self.dlneta      = tmp[:, 6]

        self.transport = tmp
        if (mag_flag == 1):
            self.names = ['nr', 'radius', 'nu', 'dlnu', 'kappa', 'dlnkappa', 'eta', 'dlneta']
        elif (mag_flag == 0):
            self.names = ['nr', 'radius', 'nu', 'dlnu', 'kappa', 'dlnkappa']

        fd.close()

class G_Avgs:
    """Rayleigh GlobalAverage Structure
    ----------------------------------
    self.niter                  : number of time steps
    self.nq                     : number of diagnostic quantities output
    self.qv[0:nq-1]             : quantity codes for the diagnostics output
    self.vals[0:niter-1,0:nq-1] : The globally averaged diagnostics 
    self.iters[0:niter-1]       : The time step numbers stored in this output file
    self.time[0:niter-1]        : The simulation time corresponding to each time step
    self.version                : The version code for this particular output (internal use)
    self.lut                    : Lookup table for the different diagnostics output
    """

    def __init__(self,filename='none',path='G_Avgs/'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
        """
        if (filename == 'none'):
            the_file = path+'00000001'
        else:
            the_file = path+filename
        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        version = swapread(fd,dtype='int32',count=1,swap=bs)
        nrec = swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)
        self.niter = nrec
        self.nq = nq
        self.qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        self.vals  = np.zeros((nrec,nq),dtype='float64')
        self.iters = np.zeros(nrec,dtype='int32')
        self.time  = np.zeros(nrec,dtype='float64')
        self.version = version
        for i in range(nrec):
            tmp = np.reshape(swapread(fd,dtype='float64',count=nq,swap=bs),(nq), order = 'F')
            self.vals[i,:] = tmp
            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)

        self.lut = get_lut(self.qv)
        fd.close()

class Shell_Avgs:
    """Rayleigh Shell Average Structure
    ----------------------------------
    self.niter                         : number of time steps
    self.nq                            : number of diagnostic quantities output
    self.nr                            : number of radial points
    self.qv[0:nq-1]                    : quantity codes for the diagnostics output
    self.radius[0:nr-1]                : radial grid

    For version 1:
    self.vals[0:nr-1,0:nq-1,0:niter-1] : The spherically averaged diagnostics
                                             

    For version 2:
    self.vals[0:n-1,0:3,0:nq-1,0:niter-1] : The spherically averaged diagnostics
                                             0-3 refers to moments (index 0 is mean, index 3 is kurtosis)    
    self.iters[0:niter-1]              : The time step numbers stored in this output file
    self.time[0:niter-1]               : The simulation time corresponding to each time step
    self.version                       : The version code for this particular output (internal use)
    self.lut                           : Lookup table for the different diagnostics output
    """
    def __init__(self,filename='none',path='Shell_Avgs/',ntheta=0):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
        """
        if (filename == 'none'):
            the_file = path+'00000001'
        else:
            the_file = path+filename
        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        version = swapread(fd,dtype='int32',count=1,swap=bs)
        nrec = swapread(fd,dtype='int32',count=1,swap=bs)
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)

        self.version = version
        self.niter = nrec
        self.nq = nq
        self.nr = nr
        self.qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        self.radius = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')
        if (self.version == 1):
            self.vals  = np.zeros((nr,nq,nrec),dtype='float64')
        if (self.version > 1):
            self.vals  = np.zeros((nr,4,nq,nrec),dtype='float64')
            #print 'version is: ', self.version
        self.iters = np.zeros(nrec,dtype='int32')
        self.time  = np.zeros(nrec,dtype='float64')

        for i in range(nrec):
            if (self.version == 1):
                tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr,swap=bs),(nr,nq), order = 'F')
                self.vals[:,:,i] = tmp
            if (self.version > 1):
                tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr*4,swap=bs),(nr,4,nq), order = 'F')
                self.vals[:,:,:,i] = tmp
            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)

        
        if ((self.version > 1) and (self.version <=3) ):
            nphi = 2*ntheta
            self.vals[:,2:4,:,:] = 0.0
            if (nphi > 0):
                cfactor = -1.0-1.0/nphi**2+2.0/nphi
                print('This ShellAverage file is version 2, but ntheta was provided.')
                print('The 2nd moment has been corrected.  3rd and 4th moments are set to zero')
                for i in range(nr):
                    self.vals[i,1,:,:] = self.vals[i,1,:,:]+cfactor*self.vals[i,0,:,:]**2
            else:
                print('This ShellAverage file is version 2, and ntheta was not provided.')
                print('The 2nd, 3rd and 4th moments are set to zero')   
                self.vals[:,1,:,:] = 0.0            

        self.lut = get_lut(self.qv)
        fd.close()

class AZ_Avgs:
    """Rayleigh AZ_Avgs Structure
    ----------------------------------
    self.niter                                    : number of time steps
    self.nq                                       : number of diagnostic quantities output
    self.nr                                       : number of radial points
    self.ntheta                                   : number of theta points
    self.qv[0:nq-1]                               : quantity codes for the diagnostics output
    self.radius[0:nr-1]                           : radial grid
    self.costheta[0:ntheta-1]                     : cos(theta grid)
    self.sintheta[0:ntheta-1]                     : sin(theta grid)
    self.vals[0:ntheta-1,0:nr-1,0:nq-1,0:niter-1] : The phi-averaged diagnostics 
    self.iters[0:niter-1]                         : The time step numbers stored in this output file
    self.time[0:niter-1]                          : The simulation time corresponding to each time step
    self.version                                  : The version code for this particular output (internal use)
    self.lut                                      : Lookup table for the different diagnostics output
    """


    def __init__(self,filename='none',path='AZ_Avgs/'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
        """
        if (filename == 'none'):
            the_file = path+'00000001'
        else:
            the_file = path+filename
        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        version = swapread(fd,dtype='int32',count=1,swap=bs)
        nrec = swapread(fd,dtype='int32',count=1,swap=bs)
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        ntheta = swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)

        self.version = version
        self.niter = nrec
        self.nq = nq
        self.nr = nr
        self.ntheta = ntheta


        self.qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        self.radius = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')
        self.costheta = np.reshape(swapread(fd,dtype='float64',count=ntheta,swap=bs),(ntheta), order = 'F')
        self.sintheta = (1.0-self.costheta**2)**0.5
        self.vals  = np.zeros((ntheta,nr,nq,nrec),dtype='float64')
        self.iters = np.zeros(nrec,dtype='int32')
        self.time  = np.zeros(nrec,dtype='float64')

        for i in range(nrec):
            tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr*ntheta,swap=bs),(ntheta,nr,nq), order = 'F')
            self.vals[:,:,:,i] = tmp
            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)

        self.lut = get_lut(self.qv)
        fd.close()

class Point_Probes:
    """Rayleigh Point Probes Structure
    ----------------------------------
    self.niter                                    : number of time steps
    self.nq                                       : number of diagnostic quantities output
    self.nr                                       : number of radial points
    self.ntheta                                   : number of theta points
    self.nphi                                     : number of phi points sampled
    self.qv[0:nq-1]                               : quantity codes for the diagnostics output
    self.radius[0:nr-1]                           : radial grid
    self.costheta[0:ntheta-1]                     : cos(theta grid)
    self.sintheta[0:ntheta-1]                     : sin(theta grid)
    self.phi[0:nphi-1]                            : phi values (radians)
    self.phi_indices[0:nphi-1]                    : phi indices (from 1 to nphi)
    self.vals[0:nphi-1,0:ntheta-1,0:nr-1,0:nq-1,0:niter-1] : The meridional slices 
    self.iters[0:niter-1]                         : The time step numbers stored in this output file
    self.time[0:niter-1]                          : The simulation time corresponding to each time step
    self.version                                  : The version code for this particular output (internal use)
    self.lut                                      : Lookup table for the different diagnostics output
    """


    def __init__(self,filename='none',path='Point_Probes/'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
        """
        if (filename == 'none'):
            the_file = path+'00000001'
        else:
            the_file = path+filename
        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        version = swapread(fd,dtype='int32',count=1,swap=bs)
        nrec = swapread(fd,dtype='int32',count=1,swap=bs)
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        ntheta = swapread(fd,dtype='int32',count=1,swap=bs)
        nphi = swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)


        self.version = version
        self.niter = nrec
        self.nq = nq
        self.nr = nr
        self.ntheta = ntheta
        self.nphi = nphi

        #print nr,ntheta,nphi,nq
        #print 'nrec is: ', nrec
        hsize = (nr+ntheta+nphi)*12 + nq*4 + 8 + 16+4
        #print 'hsize is: ', hsize
        recsize = nq*nphi*ntheta*nr*8 + 12

        #print 'expected filesize (bytes): ', recsize*nrec+hsize
        #print 'single rec size (bytes): ', recsize

        self.qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')

        self.radius = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')
        self.rad_inds = np.reshape(swapread(fd,dtype='int32',count=nr,swap=bs),(nr), order = 'F')

        self.costheta = np.reshape(swapread(fd,dtype='float64',count=ntheta,swap=bs),(ntheta), order = 'F')
        self.theta_inds = np.reshape(swapread(fd,dtype='int32',count=ntheta,swap=bs),(ntheta), order = 'F')

        self.phi = np.reshape(swapread(fd,dtype='float64',count=nphi,swap=bs),(nphi), order = 'F')
        self.phi_inds = np.reshape(swapread(fd,dtype='int32',count=nphi,swap=bs),(nphi), order = 'F')


        # convert from Fortran 1-based to Python 0-based indexing
        self.rad_inds   = self.rad_inds - 1
        self.theta_inds = self.theta_inds - 1
        self.phi_inds   = self.phi_inds - 1

        #print 'rad inds: ', self.rad_inds
        #print 'theta inds: ', self.theta_inds
        #print 'phi_inds: ', self.phi_inds
        #print 'qvals : ', self.qv
        #print ''
        #print 'radius: ', self.radius
        #print 'ctheta: ', self.costheta
        self.sintheta = (1.0-self.costheta**2)**0.5
        self.vals  = np.zeros((nphi,ntheta,nr,nq,nrec),dtype='float64')
        self.iters = np.zeros(nrec,dtype='int32')
        self.time  = np.zeros(nrec,dtype='float64')

        for i in range(nrec):
            #print i
            tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr*ntheta*nphi,swap=bs),(nphi,ntheta,nr,nq), order = 'F')
            self.vals[:,:,:,:,i] = tmp
            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)

        #print 'iters: ', self.iters
        #print 'times: ', self.time

        self.lut = get_lut(self.qv)
        fd.close()


class Meridional_Slices:
    """Rayleigh Meridional Slice Structure
    ----------------------------------
    self.niter                                    : number of time steps
    self.nq                                       : number of diagnostic quantities output
    self.nr                                       : number of radial points
    self.ntheta                                   : number of theta points
    self.nphi                                     : number of phi points sampled
    self.qv[0:nq-1]                               : quantity codes for the diagnostics output
    self.radius[0:nr-1]                           : radial grid
    self.costheta[0:ntheta-1]                     : cos(theta grid)
    self.sintheta[0:ntheta-1]                     : sin(theta grid)
    self.phi[0:nphi-1]                            : phi values (radians)
    self.phi_indices[0:nphi-1]                    : phi indices (from 1 to nphi)
    self.vals[0:nphi-1,0:ntheta-1,0:nr-1,0:nq-1,0:niter-1] : The meridional slices 
    self.iters[0:niter-1]                         : The time step numbers stored in this output file
    self.time[0:niter-1]                          : The simulation time corresponding to each time step
    self.version                                  : The version code for this particular output (internal use)
    self.lut                                      : Lookup table for the different diagnostics output
    """


    def __init__(self,filename='none',path='Meridional_Slices/'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
        """
        if (filename == 'none'):
            the_file = path+'00000001'
        else:
            the_file = path+filename
        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        version = swapread(fd,dtype='int32',count=1,swap=bs)
        nrec = swapread(fd,dtype='int32',count=1,swap=bs)
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        ntheta = swapread(fd,dtype='int32',count=1,swap=bs)
        nphi = swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)

        self.version = version
        self.niter = nrec
        self.nq = nq
        self.nr = nr
        self.ntheta = ntheta
        self.nphi = nphi


        self.qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        self.radius = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')
        self.costheta = np.reshape(swapread(fd,dtype='float64',count=ntheta,swap=bs),(ntheta), order = 'F')
        self.phi_inds = np.reshape(swapread(fd,dtype='int32',count=nphi,swap=bs),(nphi), order = 'F')
        self.phi = np.zeros(nphi,dtype='float64')
      
        # convert from Fortran 1-based to Python 0-based indexing
        self.phi_inds = self.phi_inds - 1

        dphi = (2*np.pi)/(ntheta*2)
        for i in range(nphi):
            self.phi[i] = self.phi_inds[i]*dphi

        self.sintheta = (1.0-self.costheta**2)**0.5
        self.vals  = np.zeros((nphi,ntheta,nr,nq,nrec),dtype='float64')
        self.iters = np.zeros(nrec,dtype='int32')
        self.time  = np.zeros(nrec,dtype='float64')

        for i in range(nrec):
            tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr*ntheta*nphi,swap=bs),(nphi,ntheta,nr,nq), order = 'F')
            self.vals[:,:,:,:,i] = tmp
            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)

        self.lut = get_lut(self.qv)
        fd.close()



class Equatorial_Slices:
    """Rayleigh Equatorial Slice Structure
    ----------------------------------
    self.niter                                    : number of time steps
    self.nq                                       : number of diagnostic quantities output
    self.nr                                       : number of radial points
    self.nphi                                     : number of phi points
    self.qv[0:nq-1]                               : quantity codes for the diagnostics output
    self.radius[0:nr-1]                           : radial grid
    self.vals[0:phi-1,0:nr-1,0:nq-1,0:niter-1]    : The equatorial_slices
    self.phi[0:nphi-1]                            : phi values (in radians)
    self.iters[0:niter-1]                         : The time step numbers stored in this output file
    self.time[0:niter-1]                          : The simulation time corresponding to each time step
    self.version                                  : The version code for this particular output (internal use)
    self.lut                                      : Lookup table for the different diagnostics output
    """


    def __init__(self,filename='none',path='Equatorial_Slices/'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
        """
        if (filename == 'none'):
            the_file = path+'00000001'
        else:
            the_file = path+filename
        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        version = swapread(fd,dtype='int32',count=1,swap=bs)
        nrec = swapread(fd,dtype='int32',count=1,swap=bs)
        nphi = swapread(fd,dtype='int32',count=1,swap=bs)
        nr= swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)

        self.version = version
        self.niter = nrec
        self.nq = nq
        self.nr = nr
        self.nphi= nphi


        self.qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        self.radius = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')
        self.vals  = np.zeros((nphi,nr,nq,nrec),dtype='float64')
        self.iters = np.zeros(nrec,dtype='int32')
        self.time  = np.zeros(nrec,dtype='float64')

        self.phi = np.zeros(nphi,dtype='float64')
        dphi = 2.0*np.pi/(nphi)
        for i in range(nphi):
            self.phi[i] = i*dphi

        for i in range(nrec):
            tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr*nphi,swap=bs),(nphi,nr,nq), order = 'F')
            self.vals[:,:,:,i] = tmp
            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)

        self.lut = get_lut(self.qv)
        fd.close()


class Shell_Slices:
    """Rayleigh Shell Slice Structure
    ----------------------------------
    self.niter                                    : number of time steps
    self.nq                                       : number of diagnostic quantities output
    self.nr                                       : number of shell slices output
    self.ntheta                                   : number of theta points
    self.nphi                                     : number of phi points
    self.qv[0:nq-1]                               : quantity codes for the diagnostics output
    self.radius[0:nr-1]                           : radii of the shell slices output
    self.inds[0:nr-1]                             : radial indices of the shell slices output
    self.costheta[0:ntheta-1]                     : cos(theta grid)
    self.sintheta[0:ntheta-1]                     : sin(theta grid)
    self.vals[0:nphi-1,0:ntheta-1,0:nr-1,0:nq-1,0:niter-1] 
                                                  : The shell slices 
    self.iters[0:niter-1]                         : The time step numbers stored in this output file
    self.time[0:niter-1]                          : The simulation time corresponding to each time step
    self.version                                  : The version code for this particular output (internal use)
    self.lut                                      : Lookup table for the different diagnostics output
    """

    def print_info(self, print_costheta = False):
        """ Prints all metadata associated with the shell-slice object."""
        print( 'version  : ', self.version)
        print( 'niter    : ', self.niter)
        print( 'nq       : ', self.nq)
        print( 'nr       : ', self.nr)
        print( 'ntheta   : ', self.ntheta)
        print( 'nphi     : ', self.nphi)
        print( '.......................')
        print( 'radius   : ', self.radius)
        print( '.......................')
        print( 'inds     : ', self.inds)
        print( '.......................')
        print( 'iters    : ', self.iters)
        print( '.......................')
        print( 'time     : ', self.time)
        print( '.......................')
        print( 'qv       : ', self.qv)
        if (print_costheta):
            print('.......................')
            print('costheta : ', self.costheta)

    def __init__(self,filename='none',path='Shell_Slices/',slice_spec = [], rec0 = False):
        """filename   : The reference state file to read.
           path       : The directory where the file is located (if full path not in filename
           slice_spec : Optional list of [time index, quantity code, radial index].  If 
                        specified, only a single shell is read.  time indexing and radial 
                        indexing start at 0
           rec0      : Set to true to read the first timestep's data only.
        """
        if (filename == 'none'):
            the_file = path+'00000001'
        else:
            the_file = path+filename

        #slice_spec is [time, qcode, radindex] ; time and rad_index start at 0

        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        version = swapread(fd,dtype='int32',count=1,swap=bs)
        nrec = swapread(fd,dtype='int32',count=1,swap=bs)

        ntheta = swapread(fd,dtype='int32',count=1,swap=bs)
        nphi = 2*ntheta
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)


        self.nq = nq
        self.nr = nr
        self.ntheta = ntheta
        self.nphi = nphi


        qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')

        self.lut = get_lut(qv)




        radius = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')
        inds = np.reshape(swapread(fd,dtype='int32',count=nr,swap=bs),(nr), order = 'F')
        self.costheta = np.reshape(swapread(fd,dtype='float64',count=ntheta,swap=bs),(ntheta), order = 'F')
        self.sintheta = (1.0-self.costheta**2)**0.5

        # convert from Fortran 1-based to Python 0-based indexing
        inds = inds - 1

        if (len(slice_spec) == 3):

            self.iters = np.zeros(1,dtype='int32')
            self.qv    = np.zeros(1,dtype='int32')
            self.inds  = np.zeros(1,dtype='int32')
            self.time  = np.zeros(1,dtype='float64')
            self.iters = np.zeros(1,dtype='float64')
            self.radius = np.zeros(1,dtype='float64')
            self.version = version

            self.niter = 1
            self.nq = 1
            self.nr = 1
            self.vals  = np.zeros((nphi,ntheta,1,1,1),dtype='float64')


            error = False
            tspec = slice_spec[0]
            qspec = slice_spec[1]
            rspec = slice_spec[2]

            if (tspec > (nrec-1)):
                print(" ")
                print( "---------------------------------------------------------")
                print( " Error: specified time index out of range.")
                print( " Number of records in this file : ", nrec)
                print( " Specified time index           : ", tspec)
                print( " Valid time indices range from 0 through "+str(nrec-1)+".")
                print( "---------------------------------------------------------")
                print( " " )
                error = True

            if (rspec > (nr-1)):
                error = True
                print(" ")
                print("---------------------------------------------------------")
                print(" Error: specified radial index out of range.")
                print(" Number of radii in this file     : ", nr)
                print(" Specified radial index           : ", rspec)
                print(" Valid radial indices range from 0 through "+str(nr-1)+".") 
                print( "---------------------------------------------------------")
                print( " ")
            qind = -1
            for i in range(nq):
                if (qv[i] == qspec):
                    qind = i
            if (qind == -1):
                print(" ")
                print("---------------------------------------------------------")
                print(" Error: Quantity code not found")
                print(" Specified quantity code: ", qind)
                print(" Valid quantity codes: ")
                print(" ", qv)
                print("---------------------------------------------------------")
                print(" ")
                error = True

            if (error):
                print(" Returning zero shell-slice structure.")
                fd.close()
                return

            self.lut[:] = maxq
            self.lut[qspec] = 0
            slice_size  = ntheta*nphi*8
            qsize       = nr*slice_size
            rec_size    = nq*qsize+12  # 8-byte timestamp and 4-byte time index.
            seek_offset = rec_size*tspec+qsize*qind+slice_size*rspec
            seek_bytes  = seek_offset 
            fd.seek(seek_bytes,1)


            self.radius[0] = radius[rspec]
            self.inds[0]   = inds[rspec]
            self.qv[0] = qspec
            tmp = np.reshape(swapread(fd,dtype='float64',count=ntheta*nphi,swap=bs),(nphi,ntheta), order = 'F')
            self.vals[:,:,0,0,0] = tmp
            self.time[0]  = swapread(fd, dtype='float64',count=1,swap=bs)
            self.iters[0] = swapread(fd, dtype='int32'  ,count=1,swap=bs)
        else:
            if (rec0):
                #If true, read only the first record of the file.
                nrec = 1  
            self.radius = radius
            self.inds   = inds
            self.qv     = qv
            self.niter = nrec
            self.vals  = np.zeros((nphi,ntheta,nr,nq,nrec),dtype='float64')
            self.iters = np.zeros(nrec,dtype='int32')
            self.time  = np.zeros(nrec,dtype='float64')
            self.version = version
            for i in range(nrec):
                tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr*ntheta*nphi,swap=bs),(nphi,ntheta,nr,nq), order = 'F')
                self.vals[:,:,:,:,i] = tmp
                self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
                self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)

        fd.close()

class SPH_Modes:
    """Rayleigh Shell Spectrum Structure
    ----------------------------------
    self.niter                                    : number of time steps
    self.nq                                       : number of diagnostic quantities output
    self.nr                                       : number of shell slices output
    self.nell                                     : number of ell values
    self.qv[0:nq-1]                               : quantity codes for the diagnostics output
    self.radius[0:nr-1]                           : radii of the shell slices output
    self.inds[0:nr-1]                             : radial indices of the shell slices output
    self.lvals[0:nell-1]                          : ell-values output
    self.vals[0:lmax,0:nell-1,0:nr-1,0:nq-1,0:niter-1] 
                                                  : The complex spectra of the SPH modes output
                                                  :  (here lmax denotes the maximum l-value output; not the simulation lmax)
    self.iters[0:niter-1]                         : The time step numbers stored in this output file
    self.time[0:niter-1]                          : The simulation time corresponding to each time step
    self.version                                  : The version code for this particular output (internal use)
    self.lut                                      : Lookup table for the different diagnostics output
    """




    def __init__(self,filename='none',path='SPH_Modes/'):
        """
           filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
        """
        if (filename == 'none'):
            the_file = path+'00000001'
        else:
            the_file = path+filename
        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        version = swapread(fd,dtype='int32',count=1,swap=bs)
        nrec = swapread(fd,dtype='int32',count=1,swap=bs)
        nell = swapread(fd,dtype='int32',count=1,swap=bs)
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)

        self.niter = nrec
        self.nq = nq
        self.nr = nr
        self.nell = nell


        self.qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        self.radius = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')
        self.inds = np.reshape(swapread(fd,dtype='int32',count=nr,swap=bs),(nr), order = 'F')
        self.lvals = np.reshape(swapread(fd,dtype='int32',count=nell,swap=bs),(nell), order = 'F')
        lmax = np.max(self.lvals)
        nm = lmax+1

        # convert from Fortran 1-based to Python 0-based indexing
        self.inds = self.inds - 1

        #print self.lvals
        #print lmax, nm
        #print self.inds
        self.vals  = np.zeros((nm,nell,nr,nq,nrec),dtype='complex128')
        
        self.iters = np.zeros(nrec,dtype='int32')
        self.time  = np.zeros(nrec,dtype='float64')
        self.version = version
        for i in range(nrec):
            for qv in range(nq):
                for p in range(2):
                    for rr in range(nr):
                        for lv in range(nell):
                            lval = self.lvals[lv]
                            nm = lval+1
                            tmp = np.reshape(swapread(fd,dtype='float64',count=nm,swap=bs),(nm), order = 'F')

                            if (p == 0):
                                self.vals[0:nm,lv,rr,qv,i].real = tmp
                                #if (lval == 0):
                                #    print 'real: ', tmp
                            else:
                                self.vals[0:nm,lv,rr,qv,i].imag = tmp
                                #if (lval == 0):
                                #    print 'imag: ', tmp

            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)

        if (self.version < 4):
            # The m>0 --power-- is too high by a factor of 2
            # We divide the --complex amplitude-- by sqrt(2)
            sqrttwo = np.sqrt(2)
            for k in range(nrec):
                for q in range(nq):
                    for j in range(nr):
                        for m in range(1,nm):
                            self.vals[m,:,j,q,k] = self.vals[m,:,j,q,k]/sqrttwo


        self.lut = get_lut(self.qv)
        fd.close()





class Shell_Spectra:
    """Rayleigh Shell Spectrum Structure
    ----------------------------------
    self.niter                                    : number of time steps
    self.nq                                       : number of diagnostic quantities output
    self.nr                                       : number of shell slices output
    self.nell                                     : number of ell values
    self.nm                                       : number of m values
    self.lmax                                     : maximum spherical harmonic degree l
    self.mmax                                     : maximum spherical harmonic degree m
    self.qv[0:nq-1]                               : quantity codes for the diagnostics output
    self.radius[0:nr-1]                           : radii of the shell slices output
    self.inds[0:nr-1]                             : radial indices of the shell slices output
    self.vals[0:lmax,0:mmax,0:nr-1,0:nq-1,0:niter-1] 
                                                  : The complex spectra of the shells output 
    self.lpower[0:lmax,0:nr-1,0:nq-1,0:niter-1,3]    : The power as a function of ell, integrated over m
                                                     :  index indicates (0:total,1:m=0, 2:total-m=0 power)
    self.iters[0:niter-1]                         : The time step numbers stored in this output file
    self.time[0:niter-1]                          : The simulation time corresponding to each time step
    self.version                                  : The version code for this particular output (internal use)
    self.lut                                      : Lookup table for the different diagnostics output
    """

    def print_info(self):
        """ Prints all metadata associated with the shell-spectra object."""
        print( 'version  : ', self.version)
        print( 'niter    : ', self.niter)
        print( 'nq       : ', self.nq)
        print( 'nr       : ', self.nr)
        print( 'nell     : ', self.nell)
        print( 'nm       : ', self.nm)
        print( 'lmax     : ', self.lmax)
        print( 'mmax     : ', self.mmax)
        print( '.......................')
        print( 'radius   : ', self.radius)
        print( '.......................')
        print( 'inds     : ', self.inds)
        print( '.......................')
        print( 'iters    : ', self.iters)
        print( '.......................')
        print( 'time     : ', self.time)
        print( '.......................')
        print( 'qv       : ', self.qv)


    def __init__(self,filename='none',path='Shell_Spectra/'):
        """
           filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename
        """
        if (filename == 'none'):
            the_file = path+'00000001'
        else:
            the_file = path+filename
        fd = open(the_file,'rb')
        # We read an integer to assess which endian the file was written in...
        bs = check_endian(fd,314,'int32')
        version = swapread(fd,dtype='int32',count=1,swap=bs)
        nrec = swapread(fd,dtype='int32',count=1,swap=bs)
        lmax = swapread(fd,dtype='int32',count=1,swap=bs)
        nell = lmax+1
        nm = nell   
        mmax = nm-1
        nr = swapread(fd,dtype='int32',count=1,swap=bs)
        nq = swapread(fd,dtype='int32',count=1,swap=bs)

        self.niter = nrec
        self.nq = nq
        self.nr = nr
        self.nell = nell
        self.nm   = nm
        self.lmax = lmax
        self.mmax = mmax

        self.qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        self.radius = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')
        self.inds = np.reshape(swapread(fd,dtype='int32',count=nr,swap=bs),(nr), order = 'F')
        self.vals  = np.zeros((nell,nm,nr,nq,nrec),dtype='complex128')
        
        self.iters = np.zeros(nrec,dtype='int32')
        self.time  = np.zeros(nrec,dtype='float64')
        self.version = version

        # convert from Fortran 1-based to Python 0-based indexing
        self.inds = self.inds - 1

        for i in range(nrec):

            tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr*nell*nm,swap=bs),(nm,nell,nr,nq), order = 'F')
            self.vals[:,:,:,:,i].real = tmp

            tmp2 = np.reshape(swapread(fd,dtype='float64',count=nq*nr*nell*nm,swap=bs),(nm,nell,nr,nq), order = 'F')
            self.vals[:,:,:,:,i].imag = tmp2

            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)

        if (self.version < 4):
            # The m>0 --power-- is too high by a factor of 2
            # We divide the --complex amplitude-- by sqrt(2)
            sqrttwo = np.sqrt(2)
            for k in range(nrec):
                for q in range(nq):
                    for j in range(nr):
                        for m in range(1,nm):
                            self.vals[:,m,j,q,k] = self.vals[:,m,j,q,k]/sqrttwo


        self.lut = get_lut(self.qv)
        fd.close()

        self.lpower  = np.zeros((nell,nr,nq,nrec,3),dtype='float64')
        #!Finally, we create the power
        for k in range(nrec):
            for q in range(nq):
                for j in range(nr):
                    # Load the m=0 power
                    self.lpower[:,j,q,k,1] = self.lpower[:,j,q,k,1]+np.real(self.vals[:,0,j,q,k])**2 +np.imag(self.vals[:,0,j,q,k])**2

                    # m !=0 (convective) power

                    for m in range(1,nm):
                        self.lpower[:,j,q,k,2] = self.lpower[:,j,q,k,2]+np.real(self.vals[:,m,j,q,k])**2 +np.imag(self.vals[:,m,j,q,k])**2


                    self.lpower[:,j,q,k,0] = self.lpower[:,j,q,k,2]+self.lpower[:,j,q,k,1] # total power


class Power_Spectrum():
    """Rayleigh Power Spectrum Structure
    ----------------------------------
    self.niter                                    : number of time steps
    self.nr                                       : number of radii at which power spectra are available
    self.lmax                                     : maximum spherical harmonic degree l
    self.radius[0:nr-1]                           : radii of the shell slices output
    self.inds[0:nr-1]                             : radial indices of the shell slices output
    self.power[0:lmax,0:nr-1,0:niter-1,0:2]       : the velocity power spectrum.  The third
                                                  : index indicates (0:total,1:m=0, 2:total-m=0 power)
    self.mpower[0:lmax,0:nr-1,0:niter-1,0:2]      : the magnetic power spectrum
    self.iters[0:niter-1]                         : The time step numbers stored in this output file
    self.time[0:niter-1]                          : The simulation time corresponding to each time step
    self.magnetic                                 : True if mpower exists
    """
    #Power Spectrum Class - generated using shell spectra files
    def __init__(self,infile, dims=[],power_file = False, magnetic = False, path="Shell_Spectra"):
        self.magnetic = magnetic
        if (power_file):
            self.power_file_init(infile) 
        elif (infile == 'Blank' or infile =='blank'):
            self.blank_init(dims)      
        else:
            self.spectra_file_init(path+'/'+infile)

    def blank_init(self,dims):
        print('blank init')
        self.lmax = dims[0]
        self.nr = dims[1]
        self.niter = dims[2]
        self.power = np.zeros((self.lmax+1,self.nr,self.niter,3),dtype='float64')
    def set_pars(self,iters,time,inds,radius):
        self.iters = np.zeros(self.niter,dtype='int32')
        self.time = np.zeros(self.niter,dtype='float64')

        self.inds = np.zeros(self.nr,dtype='int32')
        self.radius = np.zeros(self.nr,dtype='float64')
    
        self.iters[:]  = iters[:]
        self.time[:]   = time[:]
        self.inds[:]   = inds[:]
        self.radius[:] = radius[:]
    def power_file_init(self,pfile):
        fd = open(pfile,'rb')  
        bs = check_endian(fd,314,'int32')
        lmax  = swapread(fd,dtype='int32',count=1,swap=bs)
        nr    = swapread(fd,dtype='int32',count=1,swap=bs)
        niter = swapread(fd,dtype='int32',count=1,swap=bs)
        magint = swapread(fd,dtype='int32',count=1,swap=bs)
        if (magint == 1):
            self.magnetic = True
        else:
            self.magnetic = False
        self.iters = np.reshape(swapread(fd,dtype='int32',count=niter,swap=bs),(niter), order = 'F')
        self.time = np.reshape(swapread(fd,dtype='float64',count=niter,swap=bs),(niter), order = 'F')
        self.inds = np.reshape(swapread(fd,dtype='int32',count=nr,swap=bs),(nr), order = 'F')
        self.radius = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')
        pcount = (lmax+1)*nr*niter*3
        pdim = (lmax+1,nr,niter,3)
        self.power = np.reshape(swapread(fd,dtype='float64',count=pcount,swap=bs),pdim, order = 'F')
        if (self.magnetic):
            self.mpower = np.reshape(swapread(fd,dtype='float64',count=pcount,swap=bs),pdim, order = 'F')

        self.niter = niter
        self.nr = nr
        self.lmax = lmax
        
        fd.close()
    def write_power(self,ofile):
        fd = open(ofile,'wb') #w = write, b = binary
        dims = np.zeros(5,dtype='int32')
        dims[0] = 314
        dims[1] = self.lmax
        dims[2] = self.nr
        dims[3] = self.niter
        if (self.magnetic):
            dims[4] = 1
        else:
            dims[4] = 0
        dims.tofile(fd)
        self.iters.tofile(fd)
        self.time.tofile(fd)
        self.inds.tofile(fd)
        self.radius.tofile(fd)
        tmp = np.transpose(self.power)
        tmp.tofile(fd)
        if (self.magnetic):
            tmp = np.transpose(self.mpower)
            tmp.tofile(fd)
        fd.close()

    def spectra_file_init(self,sfile):

        a = Shell_Spectra(filename=sfile,path='./')
        lmax = a.lmax
        nr = a.nr
        nt = a.niter

        self.lmax = lmax
        self.nr = nr
        self.niter = nt
        self.radius = a.radius
        self.inds = a.inds
        self.iters = a.iters
        self.time = a.time

        # We use the lookup table to find where vr, vtheta, and vphi are stored
        vr_index = a.lut[1]
        vt_index = a.lut[2]
        vp_index = a.lut[3]
        #the last index indicates 0:full power, 1:m0 power, 2:full-m0
        power = np.zeros((lmax+1,nr,nt,3),dtype='float64')
        
        # Next we grab one radial index and one time instance of each variable
        dims = (a.nell,a.nm,nr,nt)
        vrc = np.reshape(a.vals[:,:,:, vr_index, :],dims)
        vtc = np.reshape(a.vals[:,:,:, vt_index, :],dims)
        vpc = np.reshape(a.vals[:,:,:, vp_index, :],dims)

        #print 'first index: '
        #for i in range(0,lmax+1,10):
        #    print a.vals[i:i+10,0,0,0,0]   


        #print 'second index: '
        #for i in range(0,lmax+1,10):
        #    print a.vals[0,i:i+10,0,0,0]   


        for k in range(nt):
            for j in range(nr):
                # Load the m=0 power
                power[:,j,k,1] = power[:,j,k,1]+np.real(vrc[:,0,j,k])**2 +np.imag(vrc[:,0,j,k])**2
                power[:,j,k,1] = power[:,j,k,1]+np.real(vtc[:,0,j,k])**2 +np.imag(vtc[:,0,j,k])**2
                power[:,j,k,1] = power[:,j,k,1]+np.real(vpc[:,0,j,k])**2 +np.imag(vpc[:,0,j,k])**2

                # m !=0 (convective) power

                for m in range(1,a.nm):
                    power[:,j,k,2] = power[:,j,k,2]+np.real(vrc[:,m,j,k])**2 +np.imag(vrc[:,m,j,k])**2
                    power[:,j,k,2] = power[:,j,k,2]+np.real(vtc[:,m,j,k])**2 +np.imag(vtc[:,m,j,k])**2
                    power[:,j,k,2] = power[:,j,k,2]+np.real(vpc[:,m,j,k])**2 +np.imag(vpc[:,m,j,k])**2

                power[:,j,k,0] = power[:,j,k,2]+power[:,j,k,1] # total power

        self.power = power

        if(self.magnetic):
            #Do the same thing for the magnetic field components
            # We use the lookup table to find where br, vtheta, and bphi are stored
            br_index = a.lut[401]
            bt_index = a.lut[402]
            bp_index = a.lut[403]
            #the last index indicates 0:full power, 1:m0 power, 2:full-m0
            mpower = np.zeros((lmax+1,nr,nt,3),dtype='float64')
            
            # Next we grab one radial index and one time instance of each variable
            dims = (a.nell,a.nm,nr,nt)
            brc = np.reshape(a.vals[:,:,:, br_index, :],dims)
            btc = np.reshape(a.vals[:,:,:, bt_index, :],dims)
            bpc = np.reshape(a.vals[:,:,:, bp_index, :],dims)

            for k in range(nt):
                for j in range(nr):
                    mpower[:,j,k,1] = mpower[:,j,k,1]+np.real(brc[:,0,j,k])**2 +np.imag(brc[:,0,j,k])**2
                    mpower[:,j,k,1] = mpower[:,j,k,1]+np.real(btc[:,0,j,k])**2 +np.imag(btc[:,0,j,k])**2
                    mpower[:,j,k,1] = mpower[:,j,k,1]+np.real(bpc[:,0,j,k])**2 +np.imag(bpc[:,0,j,k])**2
                    for m in range(a.nm):
                        mpower[:,j,k,0] = mpower[:,j,k,0]+np.real(brc[:,m,j,k])**2 +np.imag(brc[:,m,j,k])**2
                        mpower[:,j,k,0] = mpower[:,j,k,0]+np.real(btc[:,m,j,k])**2 +np.imag(btc[:,m,j,k])**2
                        mpower[:,j,k,0] = mpower[:,j,k,0]+np.real(bpc[:,m,j,k])**2 +np.imag(bpc[:,m,j,k])**2
                    mpower[:,j,k,2] = mpower[:,j,k,0]-mpower[:,j,k,1]
            self.mpower = mpower

def swapread(fd,dtype='float64',count=1,swap=False):
        #simple wrapper to numpy.fromfile that allows byteswapping based on Boolean swap
        if (swap):
                val = np.fromfile(fd,dtype=dtype,count=count).byteswap()
        else:
                val = np.fromfile(fd,dtype=dtype,count=count)
        if (len(val) == 1):
                val = val[0]
        return val

def check_endian(fd,sig,sigtype):
    # returns False if first element read from file matches sig
    # True otherwise
    chk = np.fromfile(fd,dtype=sigtype,count=1)
    if (chk == sig):
        return False
    else:
        return True

def build_file_list(istart,iend,path = '.',diter = -1,ndig = 8,special=False):
    files = []
    if (diter < 1):
        # Examine the directory and grab all files that fall between istart and iend
        allfiles = os.listdir(path)
        allfiles.sort()
        for f in allfiles:
            if ( ('special' in f) and special ):
                fint = int(f[0:7])
                if ( (fint >= istart ) and (fint <= iend)  ):
                    files.append(path+'/'+f)
            if ( (not 'special' in f) and not(special) ):
                fint = int(f)
                if ( (fint >= istart ) and (fint <= iend)  ):
                    files.append(path+'/'+f)
    else:
        # Generate filename manually (no ls)
        i = istart
        digmod = "%0"+str(ndig)+"d"
        while (i <= iend):
            fiter = digmod % i           
            if (special):
                fiter=fiter+'_special'
            files.append(path+'/'+fiter)
            i = i+diter
    return files

########################################################
#  These routines allow us to time averages or compile multiple diagnostic files
#  and write out a single file in the same format (for use with viz routines for example).
def Compile_GlobalAverages(file_list,ofile):
    nfiles = len(file_list)
    #   We read the first file, assume that nrec doesn't change
    #   and use the nrecs + nq in the file to create our combined array
    a = GlobalAverage(file_list[0], path = '')
    nfiles = len(file_list)
    niter_estimate = a.niter*nfiles
    nq = a.nq
    combined = np.zeros((niter_estimate,a.nq),dtype='float64')
    time = np.zeros(niter_estimate,dtype='float64')
    iters = np.zeros(niter_estimate,dtype='int32')
    ncount = 0 # total number of iterations read so far
    ind = 0

    # We open the file that we want to store the compiled time traces into and write a header
    fd = open(ofile,'wb') #w = write, b = binary
    dims = np.zeros(4,dtype='int32')
    dims[0] = 314
    dims[1] = a.version
    dims[2] = a.niter   # We will fix this at the end
    dims[3] = a.nq
    dims.tofile(fd)
    a.qv.tofile(fd)

    tmp = np.zeros(a.nq,dtype='float64')
    simtime   = np.zeros(1,dtype='float64')
    iteration = np.zeros(1,dtype='int32')
    icount = np.zeros(1,dtype='int32')
    icount[0] = 0
    for i in range(nfiles):
        the_file = file_list[i]
        a = GlobalAverage(the_file,path='')
        nrec = a.niter
        for j in range(nrec):
            tmp[:] = a.vals[j,:]
            tmp.tofile(fd)
            iteration = a.iters[j]
            simtime = a.time[j]
            simtime.tofile(fd)
            iteration.tofile(fd)
            icount[0] = icount[0]+1
    fd.seek(8)
    icount.tofile(fd)   # insert the proper number of iterations
    fd.close()

def TimeAvg_AZAverages(file_list,ofile):
    nfiles = len(file_list)
    #   We read the first file, assume that nrec doesn't change
    #   and use the nrecs + nq in the file to create our combined array
    a = AzAverage(file_list[0], path = '')
    nfiles = len(file_list)


    nr = a.nr
    ntheta = a.ntheta
    nq = a.nq
    tmp = np.zeros((ntheta,nr,nq),dtype='float64')
    simtime   = np.zeros(1,dtype='float64')
    iteration = np.zeros(1,dtype='int32')
    icount = np.zeros(1,dtype='int32')
    ifinal = np.zeros(1,dtype='int32')
    tfinal = np.zeros(1,dtype='float64')
    icount[0] = 0
    i0 = a.iters[0]
    t0 = a.time[0]
    for i in range(0,nfiles):
        the_file = file_list[i]
        print('Adding '+the_file+' to the average...')
        b = AzAverage(the_file,path='')
        nrec = b.niter
        for j in range(nrec):
            tmp[0:ntheta,0:nr,0:nq] += b.vals[0:ntheta,0:nr,0:nq,j].astype('float64')

            tfinal[0] = b.time[j]
            ifinal[0] = b.iters[j]
            icount[0] = icount[0]+1
    div = np.float(icount[0])
    tmp = tmp/div

    # We open the file that we want to store the compiled time traces into and write a header
    fd = open(ofile,'wb') #w = write, b = binary
    dims = np.zeros(6,dtype='int32')
    dims[0] = 314
    dims[1] = a.version
    dims[2] = 1
    dims[3] = a.nr
    dims[4] = a.ntheta
    dims[5] = a.nq
    dims.tofile(fd)
    a.qv.tofile(fd)
    a.radius.tofile(fd)
    a.costheta.tofile(fd)


    test = np.transpose(tmp)
    test.tofile(fd)
    t0.tofile(fd)
    i0.tofile(fd)
    # The final structure is identical to a normal az_average file save for the fact that final iteration adn final time are saved
    tfinal.tofile(fd)
    ifinal.tofile(fd)
    fd.close()

def TimeAvg_ShellAverages(file_list,ofile):
    nfiles = len(file_list)
    #   We read the first file, assume that nrec doesn't change
    #   and use the nrecs + nq in the file to create our combined array
    #print file_list
    a = ShellAverage(file_list[0], path = '')
    nfiles = len(file_list)


    nr = a.nr
    nq = a.nq
    tmp = np.zeros((nr,nq),dtype='float64')
    simtime   = np.zeros(1,dtype='float64')
    iteration = np.zeros(1,dtype='int32')
    icount = np.zeros(1,dtype='int32')
    ifinal = np.zeros(1,dtype='int32')
    tfinal = np.zeros(1,dtype='float64')
    icount[0] = 0
    i0 = a.iters[0]
    t0 = a.time[0]
    for i in range(0,nfiles):
        the_file = file_list[i]
        b = ShellAverage(the_file,path='')
        nrec = b.niter
        for j in range(nrec):
            tmp[0:nr,0:nq] += b.vals[0:nr,0:nq,j].astype('float64')

            tfinal[0] = b.time[j]
            ifinal[0] = b.iters[j]
            icount[0] = icount[0]+1
    div = np.float(icount[0])
    tmp = tmp/div

    # We open the file that we want to store the compiled time traces into and write a header
    fd = open(ofile,'wb') #w = write, b = binary
    dims = np.zeros(5,dtype='int32')
    dims[0] = 314
    dims[1] = a.version
    dims[2] = 1
    dims[3] = a.nr
    dims[4] = a.nq
    dims.tofile(fd)
    a.qv.tofile(fd)
    a.radius.tofile(fd)


    test = np.transpose(tmp)
    test.tofile(fd)
    t0.tofile(fd)
    i0.tofile(fd)
    # The final structure is identical to a normal az_average file save for the fact that final iteration adn final time are saved
    tfinal.tofile(fd)
    ifinal.tofile(fd)
    fd.close()

def integrate_dr(radius,f):
    n_r = len(radius)
    weight = np.zeros(n_r,dtype='float64')
    fpr = np.zeros(n_r,dtype='float64')
    weight[:] = 1.0
    fpr[:] = 1.0
    dr = 0.5*(radius[0]-radius[1])
    intf = dr*f[0]*fpr[0]*weight[0]
    dr = 0.5*(radius[n_r-2]-radius[n_r-1])*fpr[n_r-1]
    intf = intf+dr*f[n_r-1]*weight[n_r-1]

    for i in range(1,n_r-1):
        dr0 = 0.5*(radius[i]-radius[i+1])
        dr1 = 0.5*(radius[i-1]- radius[i])
        intf = intf+(dr1+dr0)*f[i]*fpr[i]*weight[i]
    return intf

def swapwrite(val,fd,swap=False,verbose=False, array = False):
        #simple wrapper to numpy.tofile that allows byteswapping based on Boolean swap
        #set swap to true to write bytes in different endianness than current machine

        if (swap):
                if (verbose):
                    print("Swapping on write.")
                if (array):
                    if (verbose):
                        print("Swapping entire array of bytes")
                    val2 = val.byteswap().newbyteorder() 

                    tmp = np.transpose(val2)
                    tmp.tofile(fd)        
                else:    
                    val2 = val.newbyteorder()
                    val2.tofile(fd)
        else:
            if (array):
                tmp = np.transpose(val)
                tmp.tofile(fd)
            else:
                val.tofile(fd)

###########################################################################
#  This portion file contains utilities useful for plotting the AZ-Average files


def get_lims(arr,boundstype='minmax',boundsfactor=1,themin=True):
    import numpy as np
    if (themin):
        if (boundstype == 'minmax'):
            val=arr.min()*boundsfactor
        elif (boundstype == 'rms'):
            val = -np.std(arr)*boundsfactor
    else:
        if(boundstype == 'minmax'):
            val=max(abs(arr.min()),arr.max())*boundsfactor
        elif (boundstype == 'rms'):
            val = np.std(arr)*boundsfactor
    return val
def plot_azav(fig,ax,field,radius,costheta,sintheta,r_bcz=0.71,mini=-1,maxi=-1,mycmap='jet',cbar=True, 
    boundsfactor = 1, boundstype = 'minmax', units = '',fontsize = 12, underlay = [0], nlevs = 6):
    import numpy as np
    import pylab as p 
    import matplotlib.pyplot as plt
    from matplotlib import ticker, font_manager
    #Modified version of Antoine Strukarek's routine
    #r = radius/6.9599e10
    r = radius/np.max(radius)
    n_r=len(r)
    n_t=len(costheta)
    rtmp = r.reshape(1,n_r)
    cthtmp = costheta.reshape(n_t,1)
    sthtmp = sintheta.reshape(n_t,1)
    xr = p.dot(cthtmp,rtmp)
    yr = p.dot(sthtmp,rtmp)

    if (mini == -1):
        mini = get_lims(field,boundsfactor=boundsfactor,boundstype=boundstype,themin = True)
    if (maxi == -1):
        maxi = get_lims(field,boundsfactor=boundsfactor,boundstype=boundstype,themin = False)
    if (len(underlay)!=1):
        umini = get_lims(underlay,boundsfactor=boundsfactor,boundstype=boundstype,themin = True)
        umaxi = get_lims(underlay,boundsfactor=boundsfactor,boundstype=boundstype,themin = False)


    #plt.hold(True)
    if (len(underlay) == 1):
        img = ax.pcolormesh(yr,xr,field,cmap=mycmap)
    else:
        img = ax.pcolormesh(yr,xr,underlay,cmap=mycmap)
    #ax.plot(r_bcz*sintheta,r_bcz*costheta,'k--',[0,1],[0,0],'k--')
    ax.axis('equal')
    ax.axis('off')
    #cbar=False
    if (cbar):
        cbar = plt.colorbar(img,orientation='horizontal', shrink=0.5, aspect = 15, ax=ax)
        cbar.set_label(units)
        
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.ax.tick_params(labelsize=fontsize)   #font size for the ticks

        t = cbar.ax.xaxis.label
        t.set_fontsize(fontsize)  # font size for the axis title
    if (len(underlay) == 1):
       img.set_clim((mini,maxi))
    else:
        img.set_clim((umini,umaxi))
    ax.set_xlim((0,1))
    ax.set_ylim((-1,1))
    ax.plot(r[0]*sintheta,r[0]*costheta,'k')
    ax.plot(r[n_r-1]*sintheta,r[n_r-1]*costheta,'k')
    ax.plot([0,0],[-r[n_r-1],r[n_r-1]],'k--')
    ax.plot([0,0],[-r[0],-r[n_r-1]],'k',[0,0],[r[n_r-1],r[0]],'k')

    levs=mini+np.linspace(1,nlevs,nlevs)/float(nlevs)*(maxi-mini)
    ax.contour(yr,xr,field,colors='w',levels=levs)

def streamfunction(vr,vt,r,cost,order=0):
    """------------------------------------------------------------
    This routine takes as input a divergenceless axisymmetric 
    vector field in spherical coordinates and computes from 
    it a streamfunction (a.k.a. a flux flunction).  The grid
    is decribed by r and costheta and can be non-uniform.
   ------------------------------------------------------------
    INPUTS:
   
    Vr, Vtheta = the 2-d vector velocity (or magnetic) field.
                 Dimensions are (N_Theta,N_R)
    r,cost     = the radius and cos(colatitude) of the grid.
                 r is assumed to vary from rmax to rmin and 
                 costheta from  1 to -1 (i.e. 90 degrees
                 to -90 degrees in latitude).
                 Dimensions are r(N_R), costheta(N_Theta)
    order      = If greater than zero, integration begins at the
                 outer shell and the north pole and proceeds
                 inward and southward.  If less than zero,
                 integration begins at the inner shell and 
                 south pole and proceeds upward and northward.
                 If equal to zero, both are done and an average
                 is taken.
   ------------------------------------------------------------
    OUTPUTS:
   
    psi = the streamfunction
   ------------------------------------------------------------
    """
    import numpy

    (n_t,n_r)=vr.shape
    nr = n_r
    dtheta = numpy.zeros(n_t)
    dr     = numpy.zeros(n_r)

    psi = numpy.zeros((n_t,n_r))

    dpsi_dr = numpy.zeros((n_t,n_r))
    dpsi_dt = numpy.zeros((n_t,n_r))

    theta = numpy.arccos(cost)
    sint  = numpy.sqrt(1.0-cost**2)

    for i in range(0,n_t):
        dpsi_dr[i,:] = -r[:]*sint[i]*vt[i,:]
        dpsi_dt[i,:] = r[:]*r[:]*sint[i]*vr[i,:]

    if (order >= 0):
        # double precision accumulation
        dtheta[1:n_t] = theta[1:n_t]-theta[0:n_t-1]
        dr[1:n_r] = r[1:n_r]-r[0:n_r-1]

        dtheta[0]=0 
        dr[0]=0

        for i in range(1,nr):
            psi[1:n_t,i] = psi[1:n_t,i-1] + dpsi_dr[1:n_t,i]*dr[i]
        for i in range(1,n_t):
            psi[i,1:n_r] = psi[i-1,1:n_r] + dpsi_dt[i,1:n_r]*dtheta[i]

    if (order <= 0):
        psi2=numpy.zeros((n_t,n_r))
        
        dtheta[0:n_t-1] = theta[0:n_t-1]-theta[1:n_t]
        dr[0:n_r-1] = r[0:n_r-1]-r[1:n_r]
        
        dtheta[n_t-1]=0 
        dr[n_r-1]=0
        
        for i in range(0,n_r-1,-1):
            psi[0:n_t-1,i] = psi[0:n_t-1,i+1] + dpsi_dr[0:n_t-1,i]*dr[i]
        for i in range(0,n_t-1,-1):
            psi[i,0:n_r-1] = psi[i+1,0:n_r-1] + dpsi_dt[i,0:n_r-1]*dtheta[i]
        
        if (order < 0):
            return psi2
        else:
            psi=0.5*(psi+psi2)
            
    return psi
