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
from collections import OrderedDict

maxq = 30000

def get_lut(quantities):
    """return the lookup table based on the quantity codes"""
    nq = len(quantities)
    lut = np.zeros(maxq) + maxq
    for i,q in enumerate(quantities):
        if ((0 <= q) and ( q <= maxq-1)): # quantity must be in [0, maxq-1]
            lut[q] = i
    return lut.astype('int32')

class main_input:
    """ 
    Class for working with main_input files.
    
    Attributes:
        vals : Dictionary that reflects the main_input file structure.
                   
    Initializiation syntax:
        a = main_input(file)  , where
        file = 'main_input', 'jobinfo.txt' or any similar file containing Fortran namelists.
        
    """
    def __init__(self,file):
        """ 
        Initialization routine.
        
        Input parameters:
            file : A main_input or jobinfo.txt file from which the namelist structure
                   and values may be read.
        """
        self.vals = OrderedDict()
        self.namelists = []
        self.read(file)
        self.namelists = list(self.vals.keys())
                    
    def __repr__(self, verbose = False):
        """
        Provides 'print()' functionality for the main_input class.
        """
        self.write(verbose = verbose)
        return ""

    def unset(self):
        """
        Resets all keys in self.vals to the null string "".
        """
        for nml in self.namelists:
            for var in self.vals[nml].keys():
                self.vals[nml][var] =""  

                
    def set(self,nml="",var="",val="", force = False):
        """
            Sets the value of a key in vals and checks for valid key names.
            
            Parameters:
                nml : the namelist to be modified
                var : the variable within nml to modify
                val : the value to which vals[nml][var] will be set
                
                force (optional) : When set to True, nml and/or var will
                                   be created on the fly if they do not exist.
                                   
            Calling syntax:
                mi.set(nml = 'problemsize', var ='n_r', val = 128)
        
        """
        nml_lower = nml.strip().lower()
        var_lower = var.strip().lower()
        
        if (force):
            if nml_lower not in self.namelists:
                self.vals[nml_lower] = OrderedDict()
                self.namelists = list(self.vals.keys())
                
            var_names = list(self.vals[nml_lower].keys())
            
            if var_lower not in var_names:
                self.vals[nml_lower][var_lower]=""
        
        if nml_lower in self.namelists:
            var_names = list(self.vals[nml_lower].keys())
            if var_lower in var_names:
                self.vals[nml_lower][var_lower]=val
            else:
                print('\nError:  variable '+var+' is not present in '+nml_lower+' namelist\n')
        else:
            print('\nError: '+nml_lower+' is not a valid namelist name.\n')
        
        
    def write(self, verbose = False, file=None, ndecimal=6, namelist=None):
        """
        Displays the contents of the vals dictionary and generates
        a new main_input file if desired.  Floats are expressed in
        scientific notation.
        
        Input parameters:
            verbose (optional)  : When set to True, all values, including null strings
                                  are displayed.  When set to False (default), only
                                  non-empty dictionary entries are displayed.
                                 
            file (optional)     : If a value for file is provided, the contents of
                                  the vals dictionary will be written to file in
                                  Fortran-readable namelist format.
                                 
            ndecimal (optional) : Number of digits to the right of decimal place
                                  to maintain when expressing floats in sci-not.
        
        """

        if type(ndecimal) is not type(6):
            ndecimal = 6
        dstr = str(ndecimal)
        
        lprint = print
        endl=""
        if file is not None:
            fd = open(file, "w")
            lprint = fd.write
            endl="\n"
        
        if (namelist in list(self.vals.keys())):
            namelists = [namelist]
        else:
            namelists = self.namelists
        
        
        for nml in namelists:
            
            nml_line="&"+nml+"_namelist"+endl
            lprint(nml_line)             

            for var in self.vals[nml].keys():
                val = self.vals[nml][var]          
                if type(val) is type(3.14):
                    fstring = "{:."+dstr+"e}"
                    val = fstring.format(val)
                if type(val) is not type('astring'):
                    val = str(val)
                if ((val != "") or verbose):
                    val_line = var+' = '+val+endl
                    lprint(val_line)

            lprint("/"+endl)
            
        if file is not None:
            fd.close()

    def read_file_lines(self,filename):
        """
        
           Helper routine.  Reads a file and returns it as a list of
           strings with one file line per list element.
           
           Input parameters:
               filename :  The file to be read
           
        """
        fd = open(filename, "r")
        flines = []
        while True:                            
            oneline = fd.readline()   
            if len(oneline) == 0:              
                break                         
            else:
                flines.append(oneline)        
        fd.close()
        return flines
    
    def read(self, infile):
        """
        
           Populates the vals dictionary using values obtained from a main_input file. 
        
           Input parameters:
               infile     :  A main_input or jobinfo.txt file from which to extract 
                             namelist information.
                             
           Calling Syntax:
               mi.read('main_input')
               
           Note:
               Only namelist variables specified within infile are assigned
               an associated key/value combination.  If a variable is not defined in
               infile, it's corresponding key in the vals dictionary will remain 
               undefined.  If a value for the main_input variable already exists in 
               vals, it will be overwritten using the associated value from infile. 
        
        """
        
        
        flines = self.read_file_lines(infile)

        #####################################
        # Iterate over the main_input file and 
        # extract namelist and associated 
        # variable names 
        
        nlines = len(flines)
        nread = 0
        
        while (nread < nlines):
            
            nextline = flines[nread].strip().lower() #strip white-space and force lower case
            nread+=1
 
            if (len(nextline) != 0):
                
                # Test for beginning of new namelist
                # (first character of line is &)
                tchar = nextline[0]
                if (tchar == '&'):  
                    
                    # Extract namelist name. 
                    # Update the namelist list as we go.
                    
                    nml_name = nextline[1:].split('namelist')[0][:-1]
                    if nml_name not in self.namelists:
                        self.vals[nml_name] = OrderedDict()
                        self.namelists = list(self.vals.keys())
                        
                    # Read all variable names until we reach the end of the namelist
                    reading_nml = True
                    while (reading_nml):
                        nextline = flines[nread].strip()
                        nread +=1
                        if(len(nextline) != 0):
                            if (nextline[0]=='/'):  # A forward slash marks the end of a namelist
                                reading_nml = False

                            # Extract the variable names and values. Strip out comments.
                            # Some lines in jobinfo.txt are split across multiple lines.
                            # Only those with the '=' sign contain variable names
                            lsep = nextline.split('=')                        
                            if (len(lsep) > 1 and reading_nml and (lsep[0][0] != '!')):
                                var_name = lsep[0].strip().lower()
                                var_val = lsep[1].strip().split('!')[0]
                                self.vals[nml_name][var_name] = var_val

    

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
          self.vals[index] = np.reshape(self.vals[index],(nphi,ntheta,nr), order = 'F')

class RayleighTiming:

    def __init__(self,filename,byteswap=True):
        """filename  : The timing file file to read.
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
                      'Total Runtime',
                      'Transpose 1a2a pre', 'Transpose 1a2a all-to-all', 'Transpose 1a2a post',
                      'Transpose 2a3a pre', 'Transpose 2a3a all-to-all', 'Transpose 2a3a post',
                      'Transpose 3b2b pre', 'Transpose 3b2b all-to-all', 'Transpose 3b2b post',
                      'Transpose 2b1b pre', 'Transpose 2b1b all-to-all', 'Transpose 2b1b post',
                     ]

class RayleighProfile:
    """Rayleigh Profile Structure
    ----------------------------------
    self.nr         : number of radial points
    self.nq         : number of quantities in the 2-D structure file
    self.radius      : radial coordinates
    self.vals        : vals[0:nr-1,0:nq-1]

    """

    def __init__(self,filename='none'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename)
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
           path      : The directory where the file is located (if full path not in filename)
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
    self.nr          : number of radial points
    self.radius      : radial coordinates
    self.density     : density
    self.dlnrho      : logarithmic derivative of density
    self.d2lnrho     : d_by_dr of dlnrho
    self.pressure    : pressure (only before Jul 2019)
    self.temperature : temperature
    self.dlnt        : logarithmic derivative of temperature
    self.dsdr        : entropy gradient (radial)
    self.entropy     : entropy (only before Jul 2019)
    self.gravity     : gravity (only before Jul 2019)
    self.heating     : volumetric heating (Q) (only after Jan 2019)
    """

    def __init__(self,filename='reference',path='./'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename)
        """
        the_file = os.path.join(path, filename)
        fd = open(the_file, 'rb')

        bs = check_endian(fd, 314, 'int32')

        nr = swapread(fd,dtype='int32',count=1,swap=bs)

        # there are 3 "versions" of the reference file
        #   1) has heating and entropy and gravity and pressure
        #   2) has entropy and gravity and pressure, no heating
        #   3) has heating, no entropy/gravity/pressure
        # check which version here:
        try:
            tmp = np.reshape(swapread(fd, dtype='float64', count=11*nr, swap=bs),
                             (nr,11), order='F')
            has_heating = True
            has_other = True
        except:
            try:
                fd.close()
                fd = open(the_file, 'rb')
                dummy = swapread(fd, dtype='int32', count=1, swap=bs)
                dummy = swapread(fd, dtype='int32', count=1, swap=bs)
                tmp = np.reshape(swapread(fd, dtype='float64', count=10*nr, swap=bs),
                                 (nr,10), order='F')
                has_heating = False
                has_other = True
            except:
                fd.close()
                fd = open(the_file, 'rb')
                dummy = swapread(fd, dtype='int32', count=1, swap=bs)
                dummy = swapread(fd, dtype='int32', count=1, swap=bs)
                tmp = np.reshape(swapread(fd, dtype='float64', count=8*nr, swap=bs),
                                 (nr,8), order='F')
                has_heating = True
                has_other = False
        fd.close()

        self.ref     = tmp
        self.nr      = nr
        self.radius  = tmp[:,0]
        self.density = tmp[:,1]
        self.dlnrho  = tmp[:,2]
        self.d2lnrho = tmp[:,3]

        if (has_heating and has_other):
            self.pressure    = tmp[:,4]
            self.temperature = tmp[:,5]
            self.dlnt        = tmp[:,6]
            self.dsdr        = tmp[:,7]
            self.entropy     = tmp[:,8]
            self.gravity     = tmp[:,9]
            self.heating     = tmp[:,10]
            self.names = ['radius', 'density', 'dlnrho', 'd2lnrho', 'pressure',
                          'temperature', 'dlnt', 'dsdr', 'entropy', 'gravity', 'heating']
        elif (has_other):
            self.pressure    = tmp[:,4]
            self.temperature = tmp[:,5]
            self.dlnt        = tmp[:,6]
            self.dsdr        = tmp[:,7]
            self.entropy     = tmp[:,8]
            self.gravity     = tmp[:,9]
            self.names = ['radius', 'density', 'dlnrho', 'd2lnrho', 'pressure',
                          'temperature', 'dlnt', 'dsdr', 'entropy', 'gravity']
        elif (has_heating):
            self.temperature = tmp[:,4]
            self.dlnt        = tmp[:,5]
            self.dsdr        = tmp[:,6]
            self.heating     = tmp[:,7]
            self.names = ['radius', 'density', 'dlnrho', 'd2lnrho', 'temperature',
                          'dlnt', 'dsdr', 'heating']
        else:
            print("Should produce an error on reading if you made it here")

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

    def __init__(self,filename='transport',path='./'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename)
        """
        the_file = os.path.join(path, filename)

        fd = open(the_file, 'rb')

        bs = check_endian(fd, 314, 'int32')

        nr = swapread(fd, dtype='int32', count=1, swap=bs)
        mag_flag = swapread(fd, dtype='int32', count=1, swap=bs)

        if (mag_flag == 0):
            tmp = np.reshape(swapread(fd, dtype='float64', count=5*nr, swap=bs),
                             (nr,5), order='F')
        elif (mag_flag == 1):
            tmp = np.reshape(swapread(fd, dtype='float64', count=7*nr, swap=bs),
                             (nr,7), order='F')
        else:
            print("Should never be here")

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
            self.names = ['nr', 'radius', 'nu', 'dlnu', 'kappa', 'dlnkappa',
                          'eta', 'dlneta']
        elif (mag_flag == 0):
            self.names = ['nr', 'radius', 'nu', 'dlnu', 'kappa', 'dlnkappa']
        else:
            print("Should never be here")

        fd.close()

class PDE_Coefficients:
    """Rayleigh PDE Coefficients Structure
    --------------------------------------
    self.nconst                      : number of functions
    self.nfunc                       : number of coefficients
    self.nr                          : number of radial points
    self.radius                      : radial coordinates
    self.cset[0:nconst-1]            : determines if the constant is set or not
    self.fset[0:nfunc-1]             : determines if the function is set or not
    self.constants[0:nconst-1]       : the constant coefficients
    self.functions[0:nr-1,0:nfunc-1] : the non-constant coefficients
    self.version                     : version number of coefficients format
    self.density, self.rho           : [0:nr-1] alias for func_1
    self.dlnrho                      : [0:nr-1] alias for func_8
    self.d2lnrho                     : [0:nr-1] alias for func_9
    self.temperature, self.T         : [0:nr-1] alias for func_4
    self.dlnT                        : [0:nr-1] alias for func_10
    self.dsdr                        : [0:nr-1] alias for func_14
    self.heating                     : [0:nr-1] alias for func_6 * const_10 / func_1 / func_4
    self.nu                          : [0:nr-1] alias for func_3
    self.dlnu                        : [0:nr-1] alias for func_11
    self.kappa                       : [0:nr-1] alias for func_5
    self.dlnkappa                    : [0:nr-1] alias for func_12
    self.eta                         : [0:nr-1] alias for func_7
    self.dlneta                      : [0:nr-1] alias for func_13

    """

    nconst = 10
    nfunc  = 14

    def __init__(self, filename='equation_coefficients', path='./'):
        """filename  : The pde coefficient file to read.
           path      : The directory where the file is located (if full path not in filename)
        """
        the_file = os.path.join(path, filename)

        fd = open(the_file, 'rb')

        bs = check_endian(fd, 314, 'int32')
        version = swapread(fd, dtype='int32', count=1, swap=bs)
        cset = swapread(fd, dtype='int32', count=self.nconst, swap=bs)
        fset = swapread(fd, dtype='int32', count=self.nfunc,  swap=bs)
        const = swapread(fd, dtype='float64', count=self.nconst, swap=bs)
        nr = swapread(fd, dtype='int32', count=1, swap=bs)
        radius = swapread(fd, dtype='float64', count=nr, swap=bs)
        tmp = swapread(fd, dtype='float64', count=self.nfunc*nr, swap=bs)
        funcs = np.reshape(tmp, (nr,self.nfunc), order='F')

        fd.close()

        self.version = version
        self.cset = cset
        self.fset = fset
        self.constants = const
        self.nr = nr
        self.radius = radius
        self.functions = funcs

        # aliases
        self.density = self.rho = self.functions[:,1-1]
        self.dlnrho  = self.functions[:,8-1]
        self.d2lnrho = self.functions[:,9-1]

        self.temperature = self.T = self.functions[:,4-1]
        self.dlnT        = self.functions[:,10-1]

        self.dsdr = self.functions[:,14-1]

        self.heating = self.functions[:,6-1]*self.constants[10-1]/self.rho/self.T

        self.nu   = self.functions[:,3-1]
        self.dlnu = self.functions[:,11-1]
        self.kappa    = self.functions[:,5-1]
        self.dlnkappa = self.functions[:,12-1]
        self.eta    = self.functions[:,7-1]
        self.dlneta = self.functions[:,13-1]

class GridInfo:
    """Rayleigh Grid Structure
    ----------------------------------
    self.n_r        : number of radial points
    self.n_theta    : number of latitudinal points
    self.n_phi      : number of zonal points
    self.radius     : radial coordinates
    self.rweights   : radial integration weights: to average over radius with weighting r^2,
    	i.e. sum(f(radius(i))*rweights(i)) = 
    			int_(rmin)^(rmax) f(r) r^2 dr / ((1/3)(rmax^3 - rmin^3))
    self.theta      : (co)latitudinal coordinates
    self.costheta   : cos(theta)
    self.sintheta   : sin(theta)
    self.tweights   : latitudinal integration weights: to average over colatitude with 
    	weighting sin(theta), i.e. sum(f(theta(i))*tweights(i)) = 
    	int_0^pi f(theta) sin(theta) dtheta / 2   														
    self.phi        : zonal coordinates
    self.pweights   : zonal integration weight(s) (all 1/n_phi, since uniform grid)
    """

    def __init__(self, filename='none', path='./'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename)
        """
        if (filename == 'none'):
            the_file = path + 'grid_info'
        else:
            the_file = path+filename
        fd = open(the_file, 'rb')
        # We read an integer to assess which endian the file was 
        # written in...
        bs = check_endian(fd, 314, 'int32')
        nr = swapread(fd, dtype='int32', count=1, swap=bs)
        ntheta = swapread(fd, dtype='int32', count=1, swap=bs)
        nphi = swapread(fd, dtype='int32', count=1, swap=bs)
        self.nr = nr
        self.ntheta = ntheta
        self.nphi = nphi
        self.radius = swapread(fd, dtype='float64', count=nr, swap=bs)
        self.rweights = swapread(fd, dtype='float64', count=nr, swap=bs)
        self.theta = swapread(fd, dtype='float64', count=ntheta, swap=bs)
        self.costheta = swapread(fd, dtype='float64', count=ntheta,\
                swap=bs)
        self.sintheta = swapread(fd, dtype='float64', count=ntheta,\
                swap=bs)
        self.tweights = swapread(fd, dtype='float64', count=ntheta,\
                swap=bs)
        self.phi = swapread(fd, dtype='float64', count=nphi,\
                swap=bs)
        self.pweights = swapread(fd, dtype='float64', count=nphi, swap=bs)
        self.names = ['nr', 'ntheta', 'nphi', 'radius', 'rweights',\
                'theta', 'costheta', 'sintheta', 'tweights', 'phi',\
                'dphi']        

        fd.close()
        


class G_Avgs:
    """Rayleigh GlobalAverage Data Structure
    ----------------------------------
    self.niter                  : number of time steps
    self.nq                     : number of diagnostic quantities output
    self.qv[0:nq-1]             : quantity codes for the diagnostics output
    self.vals[0:niter-1,0:nq-1] : The globally averaged diagnostics 
    self.iters[0:niter-1]       : The time step numbers stored in this output file
    self.time[0:niter-1]        : The simulation time corresponding to each time step
    self.version                : The version code for this particular output (internal use)
    self.lut                    : Lookup table for the different diagnostics output

    Initialization Examples:
        (1):  Read in a single G_Avgs file
              a = G_Avgs('00000001',path='./G_Avgs/')
              
        (2):  Concatenate time-series data from multiple G_Avgs files:
              a = G_Avgs(['00000001','00000002'])

        (3):  Concatenate time-series data from multiple G_Avgs files + save data to new G_Avgs file:
              a = G_Avgs(['00000001','00000002'],ofile='my_save_file.dat')
              
        (4):  Concatenate time-series data from multiple G_Avgs files.
              Extract only quantity codes 401 and 402:
              a = G_Avgs(['00000001','00000002'], qcodes = [401,402])

    Additional Notes:
        For concatenated data:
        
            - Output files are written in identical format to a standard G_Avgs file. They may be
              read in using the same syntax described above.
              
    """

    def __init__(self,filename='none',path='G_Avgs/', ofile='none', qcodes=[], nfiles=-1):

        """
           Initializer for the G_Avgs class.
           
           Input parameters:
               filename  : (a) {string} The G_Avgs file to read.
                         : (b) {list of strings} The G_Avgs files whose time-series
                               data you wish to concatenate.
               path      : The directory where the file is located (if full path not in filename)
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)
               ofile     : {optional; string} Filename to save time-averaged data to, if desired.
        """
        # Check to see if we are compiling data from multiple files
        # This is true if (a) ofile is set or (b) filename is a list of files, rather than a single string
        # if (a) is false, but (b) is true, no output is written, but compiled output is returned        
        
        #(a)
        multiple_files = False
        mainfile=filename
        if (ofile != 'none'):
            multiple_files = True
            
        #(b)
        ltype = type([])        
        ftype=type(filename)
        if (ltype == ftype):
            multiple_files = True
        else:
            if (mainfile == 'none'):
                mainfile = path+'00000001'
            else:
                mainfile = path+mainfile             
        
        if (multiple_files):
            mainfile = filename[0]
        
        # Read the main file no matter what    
        self.read_dimensions(mainfile)   
        self.read_data(qcodes=qcodes)
        
        if (multiple_files):
            self.compile_multiple_files(filename, qcodes = qcodes,path=path)

        if (ofile != 'none'):
            self.write(ofile)
            
            
    def write(self,outfile):
        """
            Writing routine for G_Avgs class.
            
            Input parameters:
               outfile  : (a) {string} The G_Avgs file to be written.

            Calling Syntax:
                my_gavgs.write("my_file.dat")
                (writes all data contained in my_gavgs to my_file.dat, in standard G_Avgs format)

        """
        one_rec = np.dtype([ ('vals', np.float64, [self.nq,]), ('times',np.float64), ('iters', np.int32) ])
        fstruct = np.dtype([ ('fdims', np.int32, 4), ('qvals', np.int32, [self.nq,]), ('fdata', one_rec, [self.niter,]) ])
        
        odata = np.zeros((1,),dtype=fstruct)
        odata['fdims'][0,0]=314
        odata[0]['fdims'][1]=self.version
        odata[0]['fdims'][2]=self.niter
        odata[0]['fdims'][3]=self.nq
        odata[0]['qvals'][:]=self.qv[:]
        odata[0]['fdata']['times'][:]=self.time
        odata[0]['fdata']['iters'][:]=self.iters
        odata[0]['fdata']['vals'][:,:]=self.vals
        fd = open(outfile,'wb')
        odata.tofile(fd)
        fd.close()
        

    def compile_multiple_files(self,filelist,qcodes=[],path=''):
        """
           Time-series concatenation routine for the G_Avgs class.
           
           Input parameters:
               filelist  : {list of strings} The G_Avgs files to be concatenated.
               path      : The directory where the files are located (if full path not in filename)
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)
        """
        nfiles = len(filelist)
        new_nrec = self.niter*nfiles
        self.niter = new_nrec
        self.vals = np.zeros((self.niter,self.nq),dtype='float64')
        self.iters = np.zeros(self.niter          ,dtype='int32')
        self.time  = np.zeros(self.niter          ,dtype='float64')        
        
        k = 0
        for i in range(nfiles):
            a = G_Avgs(filelist[i],qcodes=qcodes,path=path)

            #If the iteration count isn't constant throughout a run,
            #we may need to resize the arrays
            if (k+a.niter > self.niter):
                self.niter = k+a.niter
                self.vals.resize((self.niter,self.nq))
                self.time.resize((self.niter))
                self.iters.resize( (self.niter))

            self.time[k:k+a.niter]   = a.time[:]
            self.iters[k:k+a.niter]  = a.iters[:]
            self.vals[k:k+a.niter,:] = a.vals[:]
            k+=a.niter
            
        # Trim off any excess zeros 
        # (in case last final was incomplete)
        self.niter = k
        self.time = self.time[0:self.niter]
        self.iters = self.iters[0:self.niter]
        self.vals = self.vals[0:self.niter,:]
        
        
    def read_data(self,qcodes = []):
        """
           Data-reading routine for the G_Avgs class.
           
           Input parameters:
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)
           
           Notes:
                - This routine does not read header information.  Read_dimensions must be called first
                  so that dimentions and the file descriptor are initialized. 
        """
    
        self.qv    = np.zeros(self.nq             ,dtype='int32')
        self.vals  = np.zeros((self.niter,self.nq),dtype='float64')
        self.iters = np.zeros(self.niter          ,dtype='int32')
        self.time  = np.zeros(self.niter          ,dtype='float64')

        
        one_rec = np.dtype([ ('vals', np.float64, [self.nq,]), ('times', np.float64), ('iters', np.int32)  ])
        fstruct = np.dtype([ ('qvals', np.int32, [self.nq,]), ('fdata', one_rec, [self.niter,]) ])
        
        if (self.byteswap):
                fdata = np.fromfile(self.fd,dtype=fstruct,count=1).byteswap()
        else:
                fdata = np.fromfile(self.fd,dtype=fstruct)
                
        self.time[:] = fdata['fdata']['times'][0,:]
        self.iters[:] = fdata['fdata']['iters'][0,:]
        self.qv[:] = fdata['qvals'][0,:]
   
        self.fd.close()
        
        self.lut = get_lut(self.qv)  # Lookup table
        
        ########################################################
        # We may want to extract a subset of the quantity codes
        if (len(qcodes) == 0):
            self.vals[:,:] = fdata['fdata']['vals'][0,:,:]
        else:
            # number of quantity codes in the file
            qget = np.array(qcodes,dtype='int32')
            self.qv = qget  # Don't update the lookup table yet
            self.nq = len(self.qv)  # number of quantity codes we will extract
            self.vals  = np.zeros((self.niter,self.nq),dtype='float64')
            for q in range(self.nq):
                qcheck = self.lut[qget[q]]
                if (qcheck < maxq):
                    self.vals[:,q] = fdata['fdata']['vals'][0,:,qcheck]
            self.lut = get_lut(self.qv)  # Rebuild the lookup table since qv has changed

            
    def read_dimensions(self,the_file,closefile=False):
        """
        
            Header-reading routine for the G_Avgs class.
            
            This routine initialized the file descriptor and reads the dimensions 
            and Endian information of a G_Avgs file only. It does not read the G_Avgs data itself.

           
            Input parameters:
               the_file    : {string} G_Avgs file whose header is to be read.
               closefile   : (Boolean; optional; default=False) Set to True to close the file 
                             after reading the header information. 
                   

        """
        self.fd = open(the_file,'rb')        
        specs = np.fromfile(self.fd,dtype='int32',count=1)
        bcheck = specs[0]       # If not 314, we need to swap the bytes
        self.byteswap = False
        if (bcheck != 314):
            specs.byteswap()
            self.byteswap = True
            
        self.version = swapread(self.fd,dtype='int32',count=1,swap=self.byteswap)
        self.niter = swapread(self.fd,dtype='int32',count=1,swap=self.byteswap)
        self.nq = swapread(self.fd,dtype='int32',count=1,swap=self.byteswap)
        if (closefile):
            self.fd.close()
   
class Shell_Avgs:
    """Rayleigh Shell-Averaged Data Structure
    ----------------------------------

    self.niter                         : number of time steps
    self.nq                            : number of diagnostic quantities output
    self.nr                            : number of radial points
    self.qv[0:nq-1]                    : quantity codes for the diagnostics output
    self.radius[0:nr-1]                : radial grid

    self.vals[0:n-1,0:3,0:nq-1,0:niter-1] : The spherically averaged diagnostics
                                            0-3 refers to moments taken over spherical shells
                                            (index 0 is mean, index 3 is kurtosis)   
                                            Note - only the mean was output in version 1 outputs. 
                                            
    self.iters[0:niter-1]              : The time step numbers stored in this output file
    self.time[0:niter-1]               : The simulation time corresponding to each time step
    self.version                       : The version code for this particular output (internal use)
    self.lut     

    self.time_averaged                 : If True, data has been time averaged.
                                         In this case, iters and time have dimension 2, and contain the
                                         initial and final times and iterations of the averaged data.

    Initialization Examples:
        (1):  Read in a single Shell_Avgs file
              a = Shell_Avgs('00000001',path='./Shell_Avgs/')
              
        (2):  Time-average data from multiple Shell_Avgs files:
              a = Shell_Avgs(['00000001','00000002'])
              
        (2):  Concatenate time-series data from multiple Shell_Avgs files:
              a = Shell_Avgs(['00000001','00000002'],time_average=False)

        (3):  Time-averaged data from multiple Shell_Avgs files + save data to new Shell_Avgs file:
              a = Shell_Avgs(['00000001','00000002'],ofile='my_save_file.dat')
              
        (4):  Concatenate time-series data from multiple Shell_Avgs files.
              Extract only quantity codes 401 and 402:
              a = Shell_Avgs(['00000001','00000002'], qcodes = [401,402], time_average = False)

    Additional Notes:
        For concatenated or averaged data:
        
            - Output files are written in identical format to a standard Shell_Avgs file. They may be
              read in using the same syntax described above.
              
            - To save concatenated or averaged data to a new file after initialization:
              a.write('filename')  [ where 'a' is the Shell_Avgs object name]
              
    """

    def __init__(self,filename='none',path='Shell_Avgs/', ofile=None, qcodes=[],time_average=True, ntheta=0, nfiles = -1, dt = -1):

        """
           Initializer for the Shell_Avgs class.
           
           Input parameters:
               filename     : (a) {string} The Shell_Avgs file to read.
                            : (b) {list of strings} The Shell_Avgs files whose time-series
                               data you wish to concatenate or time-average.
               path         : The directory where the file is located (if full path not in filename)
               qcodes       : {optional; list of ints} Quantity codes you wish to extract (default is to extract all)
               ofile        : {optional; string}; default = None Filename to save time-averaged data to, if desired.
               time_average : {optional; Boolean; default = True} Time-average data from multiple files if a list is provided.
                              If a list of files is provided and this flag is set to False, data will be concatenated instead. 
               ntheta       : {optional; int; default = 0} Set this value to correct the variance in 
                              version=2 Shell_Avgs (only mean is preserved otherwise for that version).
               nfiles    : optional -- number of files to read relative to last file 
                                       in the list (time-averaging mode only; default is all files)
               dt        : optional -- maximum time to average over, relative to
                                       final iteration of last file in the list
                                       (time-averaging mode only; default is all time)                              
        """
        # Check to see if we are compiling data from multiple files
        # This is true if (a) ofile is set or (b) filename is a list of files, rather than a single string
        # if (a) is false, but (b) is true, no output is written, but compiled output is returned        
        
        #(a)
        multiple_files = False
        mainfile=filename
        if ofile is not None:
            multiple_files = True
            
        #(b)
        ltype = type([])        
        ftype=type(filename)
        if (ltype == ftype):
            multiple_files = True
        else:
            if (mainfile == 'none'):
                mainfile = path+'00000001'
            else:
                mainfile = path+mainfile             
        
        if (multiple_files):
            mainfile = filename[0]
        
        # Read the main file no matter what    
        self.read_dimensions(mainfile)   
        self.read_data(qcodes=qcodes, ntheta=ntheta)
        
        if (multiple_files):
            if (not time_average):
                self.compile_multiple_files(filename, qcodes = qcodes,path=path,ntheta=ntheta)
            else:
                self.time_average_files(filename, qcodes = qcodes,path=path, ntheta=ntheta, nfiles = nfiles, dt = dt)
                
        if ofile is not None:
            self.write(ofile)
            
            
    def write(self,outfile):
        """
            Writing routine for Shell_Avgs class.
            
            Input parameters:
               outfile  : (a) {string} The Shell_Avgs file to be written.

            Calling Syntax:
                my_shellavgs.write("my_file.dat")
                (writes all data contained in my_shellavgs to my_file.dat, in standard Shell_Avgs format)

        """

        numdim = 5
        if (self.version > 5):
            numdim = 6
        
        if (self.version ==1):
            one_rec = np.dtype([('vals', np.float64, [self.nq,self.nr]), ('times',np.float64), ('iters', np.int32)  ])
        else:
            one_rec = np.dtype([('vals', np.float64, [self.nq, 4, self.nr]), ('times',np.float64), ('iters', np.int32)  ])   

        fstruct = np.dtype([('fdims', np.int32, numdim), ('qvals', np.int32,(self.nq)), ('radius',np.float64,(self.nr)), ('fdata', one_rec, [self.niter,])  ])
        
                        
        odata = np.zeros((1,),dtype=fstruct)
        odata['fdims'][0,0]=314
        odata[0]['fdims'][1]=self.version
        odata[0]['fdims'][2]=self.niter
        odata[0]['fdims'][3]=self.nr
        odata[0]['fdims'][4]=self.nq
        if (self.version > 5):
            odata[0]['fdims'][5]=1 #Data has been collated now.  Set npcol to 1.
	                
        odata[0]['qvals'][:]=self.qv[:]
        odata[0]['radius'][:]=self.radius[:]
        odata[0]['fdata']['times'][:]=self.time[0:self.niter]
        odata[0]['fdata']['iters'][:]=self.iters[0:self.niter]
        if (self.version != 1):
            odata[0]['fdata']['vals'][:,:,:,:]=np.transpose(self.vals[:,:,:,:])
        else:
            odata[0]['fdata']['vals'][:,:,:]=np.transpose(self.vals[:,0,:,:])	
            
        fd = open(outfile,'wb')
        odata.tofile(fd)
        if (self.time_averaged):
            self.time[1].tofile(fd)
            self.iters[1].tofile(fd)
        fd.close()
        

    def compile_multiple_files(self,filelist,qcodes=[],path='', ntheta=0):
        """
           Time-series concatenation routine for the Shell_Avgs class.
           
           Input parameters:
               filelist  : {list of strings} The Shell_Avgs files to be concatenated.
               path      : The directory where the files are located (if full path not in filename)
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)
               ntheta    : {optional; int; default = 0} Set this value to correct the variance in 
                           version=2 Shell_Avgs (only mean is preserved otherwise for that version).
           Notes:
                  
                - This routine assumes radial resolution does not change across files in filelist.

        """
            
        nfiles = len(filelist)
        new_nrec = self.niter*nfiles
        self.niter = new_nrec
        self.vals = np.zeros((self.nr,4,self.nq,self.niter),dtype='float64')
        self.iters = np.zeros(self.niter          ,dtype='int32')
        self.time  = np.zeros(self.niter          ,dtype='float64')        
        
        k = 0
        for i in range(nfiles):

            a = Shell_Avgs(filelist[i],qcodes=qcodes,path=path, ntheta=ntheta)

            #If the iteration count isn't constant throughout a run,
            #we may need to resize the arrays
            if (k+a.niter > self.niter):
                self.niter = k+a.niter
                self.time.resize((self.niter))
                self.iters.resize( (self.niter))
                # Note that using numpy's resize routine gets a little tricky here
                # due to the striping of the vals array.  Handle this 'manually'
                # for now
                vals = np.zeros((self.nr,4,self.nq,self.niter),dtype='float64')
                vals[:,:,:,0:k] = self.vals[:,:,:,0:k]
                self.vals = vals

            self.time[k:k+a.niter]   = a.time[:]
            self.iters[k:k+a.niter]  = a.iters[:]
            self.vals[:,:,:,k:k+a.niter] = a.vals[:,:,:,:]
            k+=a.niter
            
        # Trim off any excess zeros 
        # (in case last file was incomplete)
        self.niter = k
        self.time = self.time[0:self.niter]
        self.iters = self.iters[0:self.niter]
        self.vals = self.vals[:,:,:,0:self.niter]


    def time_average_files(self,filelist,qcodes=[],path='',ntheta=0, dt=-1,nfiles=-1):
        """
           Time-series concatenation routine for the Shell_Avgs class.
           
           Input parameters:
               filelist  : {list of strings} The Shell_Avgs files to be time-averaged.
               path      : The directory where the files are located (if full path not in filename)
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)
               ntheta    : {optional; int; default = 0} Set this value to correct the variance in 
                           version=2 Shell_Avgs (only mean is preserved otherwise for that version).               
               nfiles    : optional -- number of files to read relative to last file 
                                       in the list (default is all files)
               dt        : optional -- maximum time to average over, relative to
                                       final iteration of last file in the list
                                       (default is all time)
           Notes:
                  
                - This routine assumes radial resolution does not change across files in filelist.
                - This routine assumes that each file contains AT LEAST 2 timesteps.
        """
        
        numfiles = len(filelist)
        flast = -1
        if (nfiles > 0):
            flast = numfiles-nfiles
        if (flast < 0):
            flast = 0
            

        self.niter = 1
        self.vals = np.zeros((self.nr,4,self.nq,1),dtype='float64')
        self.iters = np.zeros(2          ,dtype='int32')
        self.time  = np.zeros(2          ,dtype='float64')        
        
        #Integrate in time using a midpoint method.
        last_dt = 0.0
        total_time = 0.0

        # Read the initial file (last in the list). 
        # Go ahead and store the last time and iteration in that file.
        a = Shell_Avgs(filelist[numfiles-1],qcodes=qcodes,path=path,ntheta=ntheta)
        self.iters[1] = a.iters[a.niter-1]
        self.time[1] = a.time[a.niter-1]

        for i in range(numfiles-1,flast-1,-1):

            weights = np.zeros(a.niter,dtype='float64')
            
            if (i != 0):
                #Read in the next file for time-step information
                b = Shell_Avgs(filelist[i-1],qcodes=qcodes,path=path,ntheta=ntheta)
            
            weights[a.niter-1] = 0.5*(last_dt+a.time[a.niter-1]-a.time[a.niter-2])
            for j in range(a.niter-2,0,-1):
                weights[j] = 0.5*(a.time[j+1] -a.time[j-1])
                
            if (i != 0):
                last_dt = a.time[0]-b.time[b.niter-1]
                weights[0] = 0.5*(a.time[1]-b.time[b.niter-1])
            else:
                weights[0] = 0.5*(a.time[1]-a.time[0])

            
            for j in range(a.niter):
                time_check = True
                if (dt > 0):
                    dt0 = self.time[1]-a.time[j]
                    if (dt0 > dt):
                        time_check = False
                if (time_check):
                    self.vals[:,:,:,0]+=a.vals[:,:,:,j]*weights[j]
                    total_time += weights[j]
                    self.iters[0] = a.iters[j]
                    self.time[0] = a.time[j]

            if (i != flast):
                a = b

        self.vals = self.vals/total_time
        self.version=-self.version # negative version numbers indicate time-averaging
        self.time_averaged = True
        
    def read_data(self,qcodes = [],ntheta=0):
        """
           Data-reading routine for the Shell_Avgs class.
           
           Input parameters:
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)
               ntheta    : {optional; int; default = 0} Set this value to correct the variance in 
                           version=2 Shell_Avgs (only mean is preserved otherwise for that version).
           
           Notes:
                - This routine does not read header information.  Read_dimensions must be called first
                  so that dimentions and the file descriptor are initialized. 
        """
    
        self.qv    = np.zeros(self.nq             ,dtype='int32')
        self.vals  = np.zeros((self.nr, 4, self.nq, self.niter),dtype='float64')
        self.iters = np.zeros(self.niter+self.time_averaged          ,dtype='int32')
        self.time  = np.zeros(self.niter+self.time_averaged          ,dtype='float64')
   
            
        self.radius = np.zeros(self.nr, dtype='float64')

        # Set up the record structure
        if (self.version < 6):

            if (self.version ==1):
                one_rec = np.dtype([('vals', np.float64, [self.nq,self.nr]), ('times',np.float64), ('iters', np.int32)  ])
            else:
                one_rec = np.dtype([('vals', np.float64, [self.nq, 4, self.nr]), ('times',np.float64), ('iters', np.int32)  ])   
            
        else:
            # Things are a little more complicated following the parallel I/O redo.
            # Store all data values in a 1-D array and rerrange each record at the end
            one_rec = np.dtype([('vals', np.float64, [self.nr*4*self.nq]), ('times',np.float64), ('iters', np.int32)  ])

        fstruct = np.dtype([ ('qvals', np.int32,(self.nq)), ('radius',np.float64,(self.nr)), ('fdata', one_rec, [self.niter,])  ])

        
        if (self.byteswap):
                fdata = np.fromfile(self.fd,dtype=fstruct,count=1).byteswap()
                if (self.time_averaged):
                    self.time[1] = np.fromfile(self.fd,dtype='float64',count=1).byteswap() 
                    self.iters[1] = np.fromfile(self.fd,dtype='int32',count=1).byteswap()    
        else:
                fdata = np.fromfile(self.fd,dtype=fstruct)
                if (self.time_averaged):
                    self.time[1] = np.fromfile(self.fd,dtype='float64',count=1) 
                    self.iters[1] = np.fromfile(self.fd,dtype='int32',count=1)                  

        self.fd.close()
        
        self.time[0:self.niter] = fdata['fdata']['times'][0,:self.niter]
        self.iters[0:self.niter] = fdata['fdata']['iters'][0,:self.niter]
        self.qv[:] = fdata['qvals'][0,:]
        self.radius[:] = fdata['radius'][0,:]
        
        ####################################################3
        # Reading in the values is a little more complicated since the Shell_Avgs
        # data structure has undergone a few different iterations.
        if (self.version >= 6):
            vals = np.zeros((self.nr,4,self.nq,self.niter),dtype='float64')
            for i in range(self.niter):
                rind=0
                kone = 0
                nr_base = self.nr//self.npcol
                nr_mod = self.nr % self.npcol
                for j in range(self.npcol):
                    nrout= nr_base
                    if (j < nr_mod) :
                        nrout=nrout+1
                    ktwo = kone+nrout*4*self.nq
                    tmp = np.reshape(fdata['fdata']['vals'][0,i,kone:ktwo], (nrout,4,self.nq), order = 'F'   )
                    vals[rind:rind+nrout,:,:,i] = tmp[:,:,:]
                    rind=rind+nrout            
                    kone = ktwo
        elif (self.version == 1):
            #vals = numpy.zeros((self.nr,self.nq,self.niter)
            vals0 = np.transpose(fdata['fdata']['vals'][0,:,:,:])
            vals = np.zeros((self.nr,4,self.nq,self.niter),dtype='float64')
            vals[:,0,:,:] = vals0[:,:,:]
        else:
            #vals = numpy.zeros((self.niter,self.nq,4,self.nr))
            vals = np.transpose(fdata['fdata']['vals'][0,:,:,:,:])

        self.lut = get_lut(self.qv)  # Lookup table
        
        ########################################################
        # We may want to extract a subset of the quantity codes
        if (len(qcodes) == 0):
            self.vals = vals
        else:
            # nqfile = self.nq        # number of quantity codes in the file
            qget = np.array(qcodes,dtype='int32')
            self.qv = qget  # Don't update the lookup table yet
            self.nq = len(self.qv)  # number of quantity codes we will extract
            

            self.vals  = np.zeros((self.nr,4, self.nq, self.niter),dtype='float64')            
                
            for q in range(self.nq):
                qcheck = self.lut[qget[q]]
                if (qcheck < maxq):
                    self.vals[:,:,q,:] = vals[:,:,qcheck,:]
            self.lut = get_lut(self.qv)  # Rebuild the lookup table since qv has changed


            # Version 2 had an error with the 2nd through 4th moments.  The mean was correct.
            # The variance can be corrected after the fact if the value of ntheta is provided.
            if (self.version==2):
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

            
    def read_dimensions(self,the_file,closefile=False):
        """
        
            Header-reading routine for the Shell_Avgs class.
            
            This routine initialized the file descriptor and reads the dimensions 
            and Endian information of a Shell_Avgs file only. It does not read the Shell_Avgs data itself.

           
            Input parameters:
               the_file    : {string} Shell_Avgs file whose header is to be read.
               closefile   : (Boolean; optional; default=False) Set to True to close the file 
                             after reading the header information. 
                   

        """

        self.fd = open(the_file,'rb')        
        specs = np.fromfile(self.fd,dtype='int32',count=6)
        bcheck = specs[0]       # If not 314, we need to swap the bytes
        self.byteswap = False
        if (bcheck != 314):
            specs.byteswap()
            self.byteswap = True
            
        self.version = specs[1]
        self.niter   = specs[2]
        self.nr      = specs[3]
        self.nq      = specs[4]
        if (self.version >= 6):
            self.npcol = specs[5]
        else:
            #versions < 6 do not have npcol stored.
            #rewind by 4 bytes
            self.fd.seek(-4,1)
        
        self.time_averaged = (self.version < 0)
        if (closefile):
            self.fd.close()
   



class AZ_Avgs:
    """Rayleigh Azimuthally-Averaged Data Structure
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

    self.time_averaged                 : If True, data has been time averaged.
                                         In this case, iters and time have dimension 2, and contain the
                                         initial and final times and iterations of the averaged data.

    Initialization Examples:
        (1):  Read in a single AZ_Avgs file
              a = AZ_Avgs('00000001',path='./AZ_Avgs/')
              
        (2):  Time-average data from multiple AZ_Avgs files:
              a = AZ_Avgs(['00000001','00000002'])
              
        (2):  Concatenate time-series data from multiple AZ_Avgs files:
              a = AZ_Avgs(['00000001','00000002'],time_average=False)

        (3):  Time-averaged data from multiple AZ_Avgs files + save data to new AZ_Avgs file:
              a = AZ_Avgs(['00000001','00000002'],ofile='my_save_file.dat')
              
        (4):  Concatenate time-series data from multiple AZ_Avgs files.
              Extract only quantity codes 401 and 402:
              a = AZ_Avgs(['00000001','00000002'], qcodes = [401,402], time_average = False)

    Additional Notes:
        For concatenated or averaged data:
        
            - Output files are written in identical format to a standard AZ_Avgs file. They may be
              read in using the same syntax described above.
              
            - To save concatenated or averaged data to a new file after initialization:
              a.write('filename')  [ where 'a' is the AZ_Avgs object name]
              
    """

    def __init__(self,filename='none',path='AZ_Avgs/', ofile=None, qcodes=[],time_average=True, dt = -1, nfiles=-1):

        """
           Initializer for the AZ_Avgs class.
           
           Input parameters:
               filename     : (a) {string} The AZ_Avgs file to read.
                            : (b) {list of strings} The AZ_Avgs files whose time-series
                               data you wish to concatenate or time-average.
               path         : The directory where the file is located (if full path not in filename)
               qcodes       : {optional; list of ints} Quantity codes you wish to extract (default is to extract all)
               ofile        : {optional; string}; default = None Filename to save time-averaged data to, if desired.
               time_average : {optional; Boolean; default = True} Time-average data from multiple files if a list is provided.
                              If a list of files is provided and this flag is set to False, data will be concatenated instead. 
               nfiles    : optional -- number of files to read relative to last file 
                                       in the list (time-averaging mode only; default is all files)
               dt        : optional -- maximum time to average over, relative to
                                       final iteration of last file in the list
                                       (time-averaging mode only; default is all time)                                  
        """
        # Check to see if we are compiling data from multiple files
        # This is true if (a) ofile is set or (b) filename is a list of files, rather than a single string
        # if (a) is false, but (b) is true, no output is written, but compiled output is returned        
        
        #(a)
        multiple_files = False
        mainfile=filename
        if ofile is not None:
            multiple_files = True
            
        #(b)
        ltype = type([])        
        ftype=type(filename)
        if (ltype == ftype):
            multiple_files = True
        else:
            if (mainfile == 'none'):
                mainfile = path+'00000001'
            else:
                mainfile = path+mainfile             
        
        if (multiple_files):
            mainfile = filename[0]
        
        # Read the main file no matter what    
        self.read_dimensions(mainfile)   
        self.read_data(qcodes=qcodes)
        
        if (multiple_files):
            if (not time_average):
                self.compile_multiple_files(filename, qcodes = qcodes,path=path)
            else:
                self.time_average_files(filename, qcodes = qcodes,path=path, dt = dt, nfiles = nfiles)
                
        if ofile is not None:
            self.write(ofile)
            
            
    def write(self,outfile):
        """
            Writing routine for AZ_Avgs class.
            
            Input parameters:
               outfile  : (a) {string} The AZ_Avgs file to be written.

            Calling Syntax:
                my_shellavgs.write("my_file.dat")
                (writes all data contained in my_shellavgs to my_file.dat, in standard AZ_Avgs format)

        """

        numdim = 6

        one_rec = np.dtype([('vals', np.float64, [self.nq, self.nr, self.ntheta]), ('times',np.float64), ('iters', np.int32)  ])   

        fstruct = np.dtype([ ('fdims', np.int32, numdim), ('qvals', np.int32,(self.nq)), ('radius',np.float64,(self.nr)), 
                             ('costheta',np.float64,(self.ntheta)), ('fdata', one_rec, [self.niter,])  ])
        
                        
        odata = np.zeros((1,),dtype=fstruct)
        odata['fdims'][0,0]=314
        odata[0]['fdims'][1]=self.version
        odata[0]['fdims'][2]=self.niter
        odata[0]['fdims'][3]=self.nr
        odata[0]['fdims'][4]=self.ntheta
        odata[0]['fdims'][5]=self.nq

	                
        odata[0]['qvals'][:]=self.qv[:]
        odata[0]['radius'][:]=self.radius[:]
        odata[0]['costheta'][:]=self.costheta[:]
        odata[0]['fdata']['times'][:]=self.time[0:self.niter]
        odata[0]['fdata']['iters'][:]=self.iters[0:self.niter]
        
        odata[0]['fdata']['vals'][:,:,:,:]=np.transpose(self.vals[:,:,:,:])
            
        fd = open(outfile,'wb')
        odata.tofile(fd)
        if (self.time_averaged):
            self.time[1].tofile(fd)
            self.iters[1].tofile(fd)
        fd.close()
        

    def compile_multiple_files(self,filelist,qcodes=[],path=''):
        """
           Time-series concatenation routine for the AZ_Avgs class.
           
           Input parameters:
               filelist  : {list of strings} The AZ_Avgs files to be concatenated.
               path      : The directory where the files are located (if full path not in filename)
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)

           Notes:
                - This routine is incompatibile with version=1 AZ_Avgs due to lack of 
                  moments output in that original version.  All other versions are compatible.
                  
                - This routine assumes radial resolution does not change across files in filelist.

        """

            
        nfiles = len(filelist)
        new_nrec = self.niter*nfiles
        self.niter = new_nrec
        self.vals = np.zeros((self.nr,4,self.nq,self.niter),dtype='float64')
        self.iters = np.zeros(self.niter          ,dtype='int32')
        self.time  = np.zeros(self.niter          ,dtype='float64')        
        
        k = 0
        for i in range(nfiles):

            a = AZ_Avgs(filelist[i],qcodes=qcodes,path=path)

            #If the iteration count isn't constant throughout a run,
            #we may need to resize the arrays
            if (k+a.niter > self.niter):
                self.niter = k+a.niter
                self.time.resize((self.niter))
                self.iters.resize( (self.niter))
                # Note that using numpy's resize routine gets a little tricky here
                # due to the striping of the vals array.  Handle this 'manually'
                # for now
                vals = np.zeros((self.ntheta,self.nr, self.nq,self.niter),dtype='float64')
                vals[:,:,:,0:k] = self.vals[:,:,:,0:k]
                self.vals = vals

            self.time[k:k+a.niter]   = a.time[:]
            self.iters[k:k+a.niter]  = a.iters[:]
            self.vals[:,:,:,k:k+a.niter] = a.vals[:,:,:,:]
            k+=a.niter
            
        # Trim off any excess zeros 
        # (in case last file was incomplete)
        self.niter = k
        self.time = self.time[0:self.niter]
        self.iters = self.iters[0:self.niter]
        self.vals = self.vals[:,:,:,0:self.niter]


    def time_average_files(self,filelist,qcodes=[],path='', dt=-1,nfiles=-1):
        """
           Time-series concatenation routine for the Shell_Avgs class.
           
           Input parameters:
               filelist  : {list of strings} The Shell_Avgs files to be time-averaged.
               path      : The directory where the files are located (if full path not in filename)
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)             
               nfiles    : optional -- number of files to read relative to last file 
                                       in the list (default is all files)
               dt        : optional -- maximum time to average over, relative to
                                       final iteration of last file in the list
                                       (default is all time)
           Notes:
                  
                - This routine assumes radial resolution does not change across files in filelist.
                - This routine assumes that each file contains AT LEAST 2 timesteps.
        """
        
        numfiles = len(filelist)
        flast = -1
        if (nfiles > 0):
            flast = numfiles-nfiles
        if (flast < 0):
            flast = 0
            

        self.niter = 1
        self.vals = np.zeros((self.ntheta, self.nr,self.nq,1),dtype='float64')
        self.iters = np.zeros(2          ,dtype='int32')
        self.time  = np.zeros(2          ,dtype='float64')        
        
        #Integrate in time using a midpoint method.
        last_dt = 0.0
        total_time = 0.0

        # Read the initial file (last in the list). 
        # Go ahead and store the last time and iteration in that file.
        a = AZ_Avgs(filelist[numfiles-1],qcodes=qcodes,path=path)
        self.iters[1] = a.iters[a.niter-1]
        self.time[1] = a.time[a.niter-1]

        for i in range(numfiles-1,flast-1,-1):

            weights = np.zeros(a.niter,dtype='float64')
            
            if (i != 0):
                #Read in the next file for time-step information
                b = AZ_Avgs(filelist[i-1],qcodes=qcodes,path=path)
            
            weights[a.niter-1] = 0.5*(last_dt+a.time[a.niter-1]-a.time[a.niter-2])
            for j in range(a.niter-2,0,-1):
                weights[j] = 0.5*(a.time[j+1] -a.time[j-1])
                
            if (i != 0):
                last_dt = a.time[0]-b.time[b.niter-1]
                if (a.niter > 1):
                    weights[0] = 0.5*(a.time[1]-b.time[b.niter-1])
                else:
                   weights[0] = (a.time[0]-b.time[b.niter-1])
            else:
                weights[0] = 0.5*(a.time[1]-a.time[0])

            
            for j in range(a.niter):
                time_check = True
                if (dt > 0):
                    dt0 = self.time[1]-a.time[j]
                    if (dt0 > dt):
                        time_check = False
                if (time_check):
                    self.vals[:,:,:,0]+=a.vals[:,:,:,j]*weights[j]
                    total_time += weights[j]
                    self.iters[0] = a.iters[j]
                    self.time[0] = a.time[j]

            if (i != flast):
                a = b

        self.vals = self.vals/total_time
        self.version=-self.version # negative version numbers indicate time-averaging
        self.time_averaged = True

        
        
    def read_data(self,qcodes = []):
        """
           Data-reading routine for the AZ_Avgs class.
           
           Input parameters:
               qcodes    : {optional; list of ints} Quantity codes you wish to extract (if not all)
               ntheta    : {optional; int; default = 0} Set this value to correct the variance in 
                           version=2 AZ_Avgs (only mean is preserved otherwise for that version).
           
           Notes:
                - This routine does not read header information.  Read_dimensions must be called first
                  so that dimentions and the file descriptor are initialized. 
        """
    
        self.vals  = np.zeros((self.ntheta, self.nr, self.nq, self.niter),dtype='float64')
        self.iters = np.zeros(self.niter+self.time_averaged          ,dtype='int32')
        self.time  = np.zeros(self.niter+self.time_averaged          ,dtype='float64')

        #############################################
        # Revision:
        nq = self.nq
        ntheta=self.ntheta
        niter=self.niter
        nr = self.nr
        bs = self.byteswap
        fd = self.fd
        nrec = niter

        self.qv = np.reshape(swapread(fd,dtype='int32',count=nq,swap=bs),(nq), order = 'F')
        self.radius = np.reshape(swapread(fd,dtype='float64',count=nr,swap=bs),(nr), order = 'F')
        self.costheta = np.reshape(swapread(fd,dtype='float64',count=ntheta,swap=bs),(ntheta), order = 'F')
        self.sintheta = (1.0-self.costheta**2)**0.5


        for i in range(nrec):
            tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr*ntheta,swap=bs),(ntheta,nr,nq), order = 'F')
            self.vals[:,:,:,i] = tmp
            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)

        if (self.time_averaged):
            self.time[1] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[1] = swapread(fd,dtype='int32',count=1,swap=bs) 
            
        self.fd.close()

        self.lut = get_lut(self.qv)  # Lookup table        
            
    def read_dimensions(self,the_file,closefile=False):
        """
        
            Header-reading routine for the AZ_Avgs class.
            
            This routine initialized the file descriptor and reads the dimensions 
            and Endian information of a AZ_Avgs file only. It does not read the AZ_Avgs data itself.

           
            Input parameters:
               the_file    : {string} AZ_Avgs file whose header is to be read.
               closefile   : (Boolean; optional; default=False) Set to True to close the file 
                             after reading the header information. 
                   

        """
        self.fd = open(the_file,'rb')        
        specs = np.fromfile(self.fd,dtype='int32',count=6)
        bcheck = specs[0]       # If not 314, we need to swap the bytes
        self.byteswap = False
        if (bcheck != 314):
            specs.byteswap()
            self.byteswap = True
            
        self.version = specs[1]
        self.niter   = specs[2]
        self.nr      = specs[3]
        self.ntheta  = specs[4]
        self.nq      = specs[5]

        if (self.niter > 7):
            self.niter = 7
        
        self.time_averaged = (self.version < 0)
        if (closefile):
            self.fd.close()


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
    self.rad_inds[0:self.nr-1]                    : radial indices (from the full simulation radial grid) 
                                                  : corresponding to each point in self.radius
    self.theta_inds[0:self.ntheta-1]              : theta indices (from the full simulation theta grid) 
                                                  : corresponding to each point in self.costheta
    self.phi_inds[0:self.nphi-1]                  : phi indices (from the full simulation phi grid) 
                                                  : corresponding to each point in self.phi
    self.vals[0:nphi-1,0:ntheta-1,0:nr-1,0:nq-1,0:niter-1] : The meridional slices 
    self.iters[0:niter-1]                         : The time step numbers stored in this output file
    self.time[0:niter-1]                          : The simulation time corresponding to each time step
    self.version                                  : The version code for this particular output (internal use)
    self.lut                                      : Lookup table for the different diagnostics output
    
    --Note that the indices (phi_inds,rad_inds, theta_inds) use Python's 0-based array indexing.
    --This means that if rad_inds are 1,2,5, say, then in Rayleigh they correspond to points 2,3,6 
    --on the global grid that runs from 1 through N_R.
    """

    def __init__(self,filename='none',path='Point_Probes/'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename)
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

        # hsize = (nr+ntheta+nphi)*12 + nq*4 + 8 + 16+4
        # recsize = nq*nphi*ntheta*nr*8 + 12

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
    self.phi_inds[0:nphi-1]                       : phi indices (from the full simulation phi grid) 

    self.vals[0:nphi-1,0:ntheta-1,0:nr-1,0:nq-1,0:niter-1] : The meridional slices 
    self.iters[0:niter-1]                         : The time step numbers stored in this output file
    self.time[0:niter-1]                          : The simulation time corresponding to each time step
    self.version                                  : The version code for this particular output (internal use)
    self.lut                                      : Lookup table for the different diagnostics output

    --Note that the indices (phi_inds) use Python's 0-based array indexing.
    --This means that if phi_inds are 1,2,5, say, then in Rayleigh they correspond to points 2,3,6 
    --on the global grid that runs from 1 through n_phi.
    """

    def __init__(self,filename='none',path='Meridional_Slices/'):
        """filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename)
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
           path      : The directory where the file is located (if full path not in filename)
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
    self.rad_inds[0:nr-1]                         : radial indices of the shell slices output (from the full simulation radial grid) 
                                                  : corresponding to each point in self.radius
    self.inds                                     : same as self.rad_inds (for backwards compatibility)
    self.costheta[0:ntheta-1]                     : cos(theta grid)
    self.sintheta[0:ntheta-1]                     : sin(theta grid)
    self.vals[0:nphi-1,0:ntheta-1,0:nr-1,0:nq-1,0:niter-1] 
                                                  : The shell slices 
    self.iters[0:niter-1]                         : The time step numbers stored in this output file
    self.time[0:niter-1]                          : The simulation time corresponding to each time step
    self.version                                  : The version code for this particular output (internal use)
    self.lut                                      : Lookup table for the different diagnostics output

    --Note that the indices (rad_inds = inds) use Python's 0-based array indexing.
    --This means that if rad_inds are 1,2,5, say, then in Rayleigh they correspond to points 2,3,6 
    --on the global grid that runs from 1 through N_R.
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
        print( 'rad_inds : ', self.rad_inds)
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
           path       : The directory where the file is located (if full path not in filename)
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
        rad_inds = np.reshape(swapread(fd,dtype='int32',count=nr,swap=bs),(nr), order = 'F')
        self.costheta = np.reshape(swapread(fd,dtype='float64',count=ntheta,swap=bs),(ntheta), order = 'F')
        self.sintheta = (1.0-self.costheta**2)**0.5

        # convert from Fortran 1-based to Python 0-based indexing
        rad_inds = rad_inds - 1

        if (len(slice_spec) == 3):

            self.iters = np.zeros(1,dtype='int32')
            self.qv    = np.zeros(1,dtype='int32')
            self.rad_inds  = np.zeros(1,dtype='int32')
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
            self.rad_inds[0]   = rad_inds[rspec]
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
            self.rad_inds   = rad_inds
            self.inds = rad_inds
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
    """Rayleigh SPH Mode Structure
    ----------------------------------
    self.niter                                    : number of time steps
    self.nq                                       : number of diagnostic quantities output
    self.nr                                       : number of shell slices output
    self.nell                                     : number of ell values
    self.qv[0:nq-1]                               : quantity codes for the diagnostics output
    self.radius[0:nr-1]                           : radii of the shell slices output
    self.rad_inds[0:nr-1]                         : radial indices of the shell slices output (from the full simulation radial grid) 
                                                  : corresponding to each point in self.radius
    self.inds                                     : same as self.rad_inds (for backwards compatibility)
    self.lvals[0:nell-1]                          : ell-values output
    self.vals[0:lmax,0:nell-1,0:nr-1,0:nq-1,0:niter-1] 
                                                  : The complex spectra of the SPH modes output
                                                  :  (here lmax denotes the maximum l-value output; not the simulation lmax)
    self.iters[0:niter-1]                         : The time step numbers stored in this output file
    self.time[0:niter-1]                          : The simulation time corresponding to each time step
    self.version                                  : The version code for this particular output (internal use)
    self.lut                                      : Lookup table for the different diagnostics output

    --Note that the indices (rad_inds = inds) use Python's 0-based array indexing.
    --This means that if rad_inds are 1,2,5, say, then in Rayleigh they correspond to points 2,3,6 
    --on the global grid that runs from 1 through N_R.
    """

    def __init__(self,filename='none',path='SPH_Modes/'):
        """
           filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename)
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
        self.rad_inds = np.reshape(swapread(fd,dtype='int32',count=nr,swap=bs),(nr), order = 'F')
        self.lvals = np.reshape(swapread(fd,dtype='int32',count=nell,swap=bs),(nell), order = 'F')
        lmax = np.max(self.lvals)
        nm = lmax+1

        # convert from Fortran 1-based to Python 0-based indexing
        self.rad_inds = self.rad_inds - 1
        self.inds = self.rad_inds

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
                            else:
                                self.vals[0:nm,lv,rr,qv,i].imag = tmp

            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)

        if (self.version < 4):
            # The m>0 --power-- is too high by a factor of 2
            # We divide the --complex amplitude-- by sqrt(2)
            self.vals[1:,:,:,:,:] /= np.sqrt(2.0)

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
    self.rad_inds[0:nr-1]                         : radial indices of the shell slices output (from the full simulation radial grid) 
                                                  : corresponding to each point in self.radius
    self.inds                                     : same as self.rad_inds (for backwards compatibility)
    self.vals[0:lmax,0:mmax,0:nr-1,0:nq-1,0:niter-1] 
                                                  : The complex spectra of the shells output 
    self.lpower[0:lmax,0:nr-1,0:nq-1,0:niter-1,3]    : The power as a function of ell, integrated over m
                                                     :  index indicates (0:total,1:m=0, 2:total-m=0 power)
    self.iters[0:niter-1]                         : The time step numbers stored in this output file
    self.time[0:niter-1]                          : The simulation time corresponding to each time step
    self.version                                  : The version code for this particular output (internal use)
    self.lut                                      : Lookup table for the different diagnostics output

    --Note that the indices (rad_inds = inds) use Python's 0-based array indexing.
    --This means that if rad_inds are 1,2,5, say, then in Rayleigh they correspond to points 2,3,6 
    --on the global grid that runs from 1 through N_R.
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
        print( 'rad_inds : ', self.rad_inds)
        print( '.......................')
        print( 'iters    : ', self.iters)
        print( '.......................')
        print( 'time     : ', self.time)
        print( '.......................')
        print( 'qv       : ', self.qv)


    def __init__(self,filename='none',path='Shell_Spectra/'):
        """
           filename  : The reference state file to read.
           path      : The directory where the file is located (if full path not in filename)
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
        self.rad_inds = np.reshape(swapread(fd,dtype='int32',count=nr,swap=bs),(nr), order = 'F')
        self.vals  = np.zeros((nell,nm,nr,nq,nrec),dtype='complex128')
        
        self.iters = np.zeros(nrec,dtype='int32')
        self.time  = np.zeros(nrec,dtype='float64')
        self.version = version

        # convert from Fortran 1-based to Python 0-based indexing
        self.rad_inds = self.rad_inds - 1
        self.inds = self.rad_inds

        for i in range(nrec):

            tmp = np.reshape(swapread(fd,dtype='float64',count=nq*nr*nell*nm,swap=bs),(nm,nell,nr,nq), order = 'F')
            self.vals[:,:,:,:,i].real = tmp

            tmp2 = np.reshape(swapread(fd,dtype='float64',count=nq*nr*nell*nm,swap=bs),(nm,nell,nr,nq), order = 'F')
            self.vals[:,:,:,:,i].imag = tmp2

            self.time[i] = swapread(fd,dtype='float64',count=1,swap=bs)
            self.iters[i] = swapread(fd,dtype='int32',count=1,swap=bs)

        if (self.version != 4):
            # The m>0 --power-- is too high by a factor of 2
            # We divide the --complex amplitude-- by sqrt(2)
            self.vals[:,1:,:,:,:] /= np.sqrt(2.0)

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
    self.rad_inds[0:nr-1]                         : radial indices of the shell slices output (from the full simulation radial grid) 
                                                  : corresponding to each point in self.radius
    self.inds                                     : same as self.rad_inds (for backwards compatibility)
    self.power[0:lmax,0:nr-1,0:niter-1,0:2]       : the velocity power spectrum.  The third
                                                  : index indicates (0:total,1:m=0, 2:total-m=0 power)
    self.mpower[0:lmax,0:nr-1,0:niter-1,0:2]      : the magnetic power spectrum
    self.iters[0:niter-1]                         : The time step numbers stored in this output file
    self.time[0:niter-1]                          : The simulation time corresponding to each time step
    self.magnetic                                 : True if mpower exists

    --Note that the indices (rad_inds = inds) use Python's 0-based array indexing.
    --This means that if rad_inds are 1,2,5, say, then in Rayleigh they correspond to points 2,3,6 
    --on the global grid that runs from 1 through N_R.
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
        self.rad_inds = self.inds
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
        self.rad_inds = self.inds
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
        self.rad_inds = self.inds
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
            br_index = a.lut[801]
            bt_index = a.lut[802]
            bp_index = a.lut[803]
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
    a = G_Avgs(file_list[0], path = '')
    nfiles = len(file_list)

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
        a = G_Avgs(the_file,path='')
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
    a = AZ_Avgs(file_list[0], path = '')
    nfiles = len(file_list)

    nr = a.nr
    ntheta = a.ntheta
    nq = a.nq
    tmp = np.zeros((ntheta,nr,nq),dtype='float64')
    icount = np.zeros(1,dtype='int32')
    ifinal = np.zeros(1,dtype='int32')
    tfinal = np.zeros(1,dtype='float64')
    icount[0] = 0
    i0 = a.iters[0]
    t0 = a.time[0]
    for i in range(0,nfiles):
        the_file = file_list[i]
        b = AZ_Avgs(the_file,path='')
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
    a = Shell_Avgs(file_list[0], path = '')
    nfiles = len(file_list)

    nr = a.nr
    nq = a.nq
    if (a.version == 1):
        tmp = np.zeros((nr,nq),dtype='float64')
    else:
        tmp = np.zeros((nr,4,nq),dtype='float64')        
    icount = np.zeros(1,dtype='int32')
    ifinal = np.zeros(1,dtype='int32')
    tfinal = np.zeros(1,dtype='float64')
    icount[0] = 0
    i0 = a.iters[0]
    t0 = a.time[0]
    for i in range(0,nfiles):
        the_file = file_list[i]
        b = Shell_Avgs(the_file,path='')
        nrec = b.niter
        for j in range(nrec):
            if (a.version == 1):
                tmp[0:nr,0:nq] += b.vals[0:nr,0:nq,j].astype('float64')
            else:
                tmp[0:nr,0:4,0:nq] += b.vals[0:nr,0:4,0:nq,j].astype('float64')

            tfinal[0] = b.time[j]
            ifinal[0] = b.iters[j]
            icount[0] = icount[0]+1
    div = np.float(icount[0])
    tmp = tmp/div

    # We open the file that we want to store the compiled time traces into and write a header
    ndim=6
    if (a.version < 6):
        ndim = 5
    fd = open(ofile,'wb') #w = write, b = binary
    dims = np.zeros(ndim,dtype='int32')
    dims[0] = 314
    dims[1] = a.version
    dims[2] = 1
    dims[3] = a.nr
    dims[4] = a.nq
    if (a.version >= 6):
        dims[5] = 1
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

class rayleigh_vapor:
    """Generates a vapor dataset from interpolated Rayleigh data"""
    def __init__(self,name=None,varnames=None,varfiles=None, rayleigh_root=None, 
                vapor_bin=None, nxyz=None, grid_file=None, force=False, timeout=300,
                remove_spherical_means=[], rmins=[], rmaxes=[], vapor_version=3,
                vector_names=[],vector_files=[], tempdir='.'):

        self.numts=len(varfiles)
        self.varnames=varnames
        self.data_dir = name+'_data'
        self.nvars=len(varnames)
        self.varfiles=varfiles
        self.timeout=timeout
        self.tempdir=tempdir
        
        self.nvec=len(vector_names)
        if (self.nvec > 0):
            self.vector_names = vector_names
            self.vector_files = vector_files
        
        if (vapor_version == 3):
            self.ccmd='vdccreate '
            self.pcmd='raw2vdc '
            self.vaporfile=name+'.vdc'
        else:
            self.ccmd='vdfcreate '
            self.pcmd='raw2vdf -quiet '
            self.vaporfile=name+'.vdf'
        
        if (len(remove_spherical_means) != self.nvars):
            self.remove_spherical_mean=self.nvars*False
        else:
            self.remove_spherical_mean=remove_spherical_means
        
        if (len(rmins) != self.nvars):
            self.zero_rmin=False
            self.rmins=[None]*self.nvars
        else:
            self.zero_rmin=True
            self.rmins=rmins
            
        if (len(rmaxes) != self.nvars):
            self.zero_rmax=False
            self.rmaxes=[None]*self.nvars
        else:
            self.zero_rmax=True
            self.rmaxes=rmaxes       
        
        varstring=' '
        for i in range(self.nvars):
            varstring=varstring+varnames[i]+':'
        if (self.nvec > 0):
            for vn in self.vector_names:
                for n in vn:
                    varstring=varstring+n+':'
        varstring=varstring[0:len(varstring)-1]  # remove trailing ':'
        print(varstring)
        self.varstring=varstring
        self.vapor_bin=vapor_bin
        self.rayleigh_root=rayleigh_root
        self.nxyz=nxyz
        self.grid_file=grid_file
        if force:
            print('Parameter "force" is set to true.')
            print('Removing: '+self.vaporfile+' > /dev/null')
            print('Removing: '+self.data_dir+' > /dev/null')
            self.destroy_vdc()            
    def create_dataset(self, force=False):
        import subprocess as sp
        res=str(self.nxyz)
        cube_string = res+'x'+res+'x'+res
        cmd1 = 'export PATH=$PATH:'+self.vapor_bin
        cmd2 = ' && ' 
        cmd3 = self.ccmd+' -dimension '+cube_string+' -numts '+str(self.numts)
        cmd3 = cmd3+' -vars3d '+self.varstring+' '+self.vaporfile
        creation_cmd=cmd1+cmd2+cmd3
        s=sp.Popen(creation_cmd,shell=True)
        s.wait(timeout=self.timeout)

    def populate_dataset(self):
        import subprocess as sp
        for i in range(self.numts):
            print('Converting data for timestep '+str(i)+' of '+str(self.numts-1))
            for j in range(self.nvars):
                infile=self.varfiles[i][j]
                ofile=infile+'.cube'
                ofile=self.tempdir+'/temp.cube'
                self.rayleigh_to_cube(infile,ofile,remove_spherical_mean=self.remove_spherical_mean[j], 
                                      rmin=self.rmins[j], rmax=self.rmaxes[j])
                self.cube_to_vdc(ofile,i,j)
            if (self.nvec > 0):
                xfile = self.tempdir+'/x.cube'
                yfile = self.tempdir+'/y.cube'
                zfile = self.tempdir+'/z.cube'
                mfile = self.tempdir+'/m.cube'
                for j in range(self.nvec):
                    vnames = self.vector_names[j]
                    mag=(len(vnames)==4)
                    self.rayleigh_vector_to_cube(self.vector_files[j][i],mag=mag)
                    self.cube_to_vdc(xfile,i, vnames[0])
                    self.cube_to_vdc(yfile,i, vnames[1])
                    self.cube_to_vdc(zfile,i, vnames[2])
                    if (mag):
                        self.cube_to_vdc(mfile,i,vnames[3])
        #Cleanup
        print('Cleaning up temporary files')
        if self.nvars > 0:
            cmd = 'rm -rf '+ofile+' > /dev/null'
            s=sp.Popen(cmd,shell=True)
            s.wait(timeout=self.timeout)
        if (self.nvec > 0):
            for f in [xfile,yfile,zfile,mfile]:
                cmd = 'rm -rf '+f+' > /dev/null'
                s=sp.Popen(cmd,shell=True)
                s.wait(timeout=self.timeout)
        print('Complete.')
                
    def rayleigh_to_cube(self,infile,ofile,remove_spherical_mean=False, rmin=None, rmax=None):
        import subprocess as sp
        cmd1 = 'export PATH=$PATH:'+self.rayleigh_root
        cmd2 = ' &&  interp3d -i '+infile+' -o '+ofile+' -g '+self.grid_file+' -N '+str(self.nxyz)

        if(remove_spherical_mean):
            cmd2=cmd2+" -rsm"

        if(rmin != None):
            cmd2=cmd2+" -rmin "+str(rmin)

        if(rmax != None):
            cmd2=cmd2+" -rmax "+str(rmax)
            
        cmd = cmd1+cmd2
        s=sp.Popen(cmd,shell=True)
        s.wait(timeout=self.timeout)

    def rayleigh_vector_to_cube(self,vfiles, mag=False):
        import subprocess as sp
        rf=vfiles[0]  # r-file
        tf=vfiles[1]  # theta-file
        pf=vfiles[2]  # phi-file
        xfile = self.tempdir+'/x.cube'
        yfile = self.tempdir+'/y.cube'
        zfile = self.tempdir+'/z.cube'
        mfile = self.tempdir+'/m.cube'
        cmd1 = 'export PATH=$PATH:'+self.rayleigh_root
        cmd2 = ' &&  interp3d -ir '+rf+' -it '+tf+' -ip '+pf
        cmd2 = cmd2+' -ox '+xfile+' -oy '+yfile+' -oz '+zfile+' -g '+self.grid_file+' -N '+str(self.nxyz)

        if(mag):
            cmd2=cmd2+" -om "+mfile
            
        cmd = cmd1+cmd2
        #print(cmd)
        s=sp.Popen(cmd,shell=True)
        s.wait(timeout=self.timeout)

    def cube_to_vdc(self,ofile,timeind,varind):
        import subprocess as sp
        if type(varind) is type(1):
            varname=self.varnames[varind]
        else:
            varname=varind  # string was passed
        cmd1 = 'export PATH=$PATH:'+self.vapor_bin
        cmd2 = ' && '
        cmd3 = self.pcmd+' -ts '+str(timeind)+' -varname '+varname
        cmd3 = cmd3+' '+self.vaporfile+' '+ofile
        cmd = cmd1+cmd2+cmd3
        s=sp.Popen(cmd, shell=True)
        s.wait(timeout=self.timeout)

    def destroy_vdc(self):
        import subprocess as sp
        cmd1 = 'rm -rf '+self.vaporfile +' > /dev/null'
        cmd2 = 'rm -rf '+self.data_dir+' > /dev/null'

        try:
            s=sp.Popen(cmd1,shell=True)
            s.wait(timeout=self.timeout)
            print(cmd1)
        except:
            print('cmd1 error', cmd1)
            s.communicate()
            pass

        print(cmd2)
        try:
            s=sp.Popen(cmd2,shell=True)
            s.wait(timeout=self.timeout)
            print(cmd2)
        except:
            print('cmd2 error', cmd2)
            s.communicate()


def gen_3d_filelist( qcodes, diter, istart,iend, directory='Spherical_3D', ndig=8):
    files = []
    for i in range(istart,iend+diter,diter):
        fstring="{:0>"+str(ndig)+"d}"
        istring = directory+"/"+fstring.format(i)
        f = []
        for q in qcodes:
            qfnt="{:0>"+str(4)+"d}"
            qstr= qfnt.format(q)
            f.append(istring+'_'+qstr)
        files.append(f)
    return files

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
    from matplotlib import ticker
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
        img = ax.pcolormesh(yr,xr,field,cmap=mycmap,shading='auto')
    else:
        img = ax.pcolormesh(yr,xr,underlay,cmap=mycmap,shading='auto')
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

    return img

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

###### Function for reading in the checkpoint files (i.e. P,PAB,T,TAB,W,WAB,Z,ZAB)
def checkpoint_read(chk_file_string, nr, ntheta):
    nell = (2*ntheta)//3
    shape = (nell,nell)
    i,j = np.indices(shape)
    m = i <= j
    target_all = np.zeros((nell,nell,nr),dtype="complex")
    length_half = len(np.fromfile(chk_file_string,"f8"))//2
    chunk_length = length_half//nr
    chk_file = np.fromfile(chk_file_string,"f8")
    for i in range(nr):
        target = np.zeros_like((m),dtype="complex")
        target.real[m] = chk_file[:length_half][i*chunk_length:(i+1)*chunk_length] 
        target.imag[m] = chk_file[length_half:][i*chunk_length:(i+1)*chunk_length]
        target_all[:,:,i] = target
    return target_all
