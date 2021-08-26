#Reference State Generation Tools  (reference_tools.py)
import numpy

class equation_coefficients:
    """ equation coeff class  """
    nconst = 10
    nfunc = 14
    version=1
    c_dict    = {'two_omega':1, 'buoy_fact':2, 'p_fact':3, 'lorentz_fact':4, 'visc_fact':5,
                 'diff_fact':6, 'resist_fact':7, 'nu_heat_fact':8, 'ohm_heat_fact':9,
                 'luminosity':10}
    f_dict    = {'density':1, 'buoy':2, 'nu':3, 'temperature':4, 'kappa':5, 'heating':6,
                 'eta':7, 'd_ln_rho':8, 'd2_ln_rho':9, 'd_ln_T':10, 'd_ln_nu':11, 'd_ln_kappa':12,
                 'd_ln_eta':13, 'ds_dr':14}

    def __init__(self,radius=[], file=None):
        if (len(radius) != 0):
            nr = len(radius)
            self.nr = nr
            self.radius = numpy.zeros(nr,dtype='float64')
            self.radius[:] = radius[:]
            self.functions  = numpy.zeros((self.nfunc,nr) , dtype='float64' )
            
            self.constants = numpy.zeros(self.nconst     , dtype='float64' )
            self.cset      = numpy.zeros(self.nconst     , dtype='int32'   )
            self.fset      = numpy.zeros(self.nfunc      , dtype='int32'   )
        elif (file != None):
            self.read(filename=file)

    def __getattr__(self, name):
        if name in self.f_dict:
            return self.functions[self.f_dict[name] - 1]
        elif name in self.c_dict:
            return self.constants[self.c_dict[name] - 1]
        else:
            raise AttributeError("'{}' has no attribute '{}'".format(self.__class__, name))

    def __setattr__(self, name, value):
        if name in self.f_dict:
            self.set_function(value, name)
        elif name in self.c_dict:
            self.set_constant(value, name)
        else:
            super().__setattr__(name, value)

    def set_function(self,y,f_name):

        if (isinstance(f_name,str)):
            fi = self.f_dict[f_name]
        else:
            fi = f_name

        self.functions[fi-1,:] = y
        self.fset[fi-1] = 1
        
    def set_constant(self,c,c_name):
        if (isinstance(c_name,str)):
            ci = self.c_dict[c_name]
        else:
            ci = c_name
        self.constants[ci-1] = c
        self.cset[ci-1] = 1

    def write(self, filename='ecoefs.dat'):
        pi = numpy.array([314],dtype='int32')
        nr = numpy.array([self.nr],dtype='int32')
        version = numpy.array([self.version],dtype='int32')
        fd = open(filename,'wb')
        pi.tofile(fd)
        version.tofile(fd)
        self.cset.tofile(fd)
        self.fset.tofile(fd)
        self.constants.tofile(fd)
        nr.tofile(fd)
        self.radius.tofile(fd)
        self.functions.tofile(fd)
        fd.close() 
        
    def read(self, filename='equation_coefficients'):

        fd = open(filename,'rb')
        picheck = numpy.fromfile(fd,dtype='int32',count=1)[0]
       
        self.version = numpy.fromfile(fd,dtype='int32', count=1)[0]
        self.cset = numpy.fromfile(fd,dtype='int32',count=self.nconst)
        self.fset = numpy.fromfile(fd,dtype='int32',count=self.nfunc)
        self.constants = numpy.fromfile(fd,dtype='float64',count=self.nconst)
        self.nr = numpy.fromfile(fd,dtype='int32',count=1)[0]
        self.radius = numpy.fromfile(fd,dtype='float64',count=self.nr)
        functions=numpy.fromfile(fd,dtype='float64',count=self.nr*self.nfunc)
        self.functions = numpy.reshape(functions, (self.nfunc,self.nr))
        fd.close()

class background_state:
    nr = None
    radius = None
    pressure = None
    temperature = None
    entropy = None
    pressure_gradient = None
    entropy_gradient = None
    density = None
    def __init__(self, radius, pressure = None, temperature = None,
                 density = None, entropy = None, entropy_gradient=None, 
                 pressure_gradient=None
                 ):  # :)
        arg_dict = locals() #dictionary of parameters above ( radius, pressure, etc.)
        kvp = arg_dict.items() # key-value pairs (e.g., 'pressure' : [p1,p2,...pn])
        
        # Ensure that radius was specified correctly.        
        passed_rcheck = False
        try:
            nr = len(radius)
            passed_rcheck = True
        except:
            print('Error:  Radius must be a numpy array or list.')
            print('              Returning empty data structure.')

        # Initialize the class attributes
        if (passed_rcheck):  
            self.nr = nr
            self.radius = radius

            # Loop over the optional keyword parameters passed to __init__.
            # If a parameter was specified consistently, assign its value
            # to the corresponding class attribute.
            for (k,v) in kvp:
                self.set_variable(k,v)
                            
    def set_variable(self,k,v):
        if (k != 'self'  and k != 'radius'):
            attr_consistent = False
            try:
                nv = len(v)
                if (nv == self.nr):
                    attr_consistent = True
                else:
                    print('Error: Number of elements in', k, 'disagrees with number of radial points.')
                    print('Number of elements in ',k,' is: ', nv)
                    print('Number of radial points is    : ', self.nr)
                    print('The',k,'attribute will not be initialized.')

            except:
                if (v != None):
                    print('Error: ', k, 'must be a numpy array or list.')
                    print('The',k,'attribute will not be initialized.')

            if (attr_consistent):
                setattr(self,k,v)  # sets self.k = v
                    
def compute_heating_profile(hpars, radius, htype=0, pressure = 0):
    """  
         Returns a heating profile Q(radius).  Profile returned is 
         determined by value of htype.  Radius is assumed to be in
         ascending order (unlike that used in Rayleigh).
         
         Q(r) is normalized such that:
             integral_rmin_rmax { Q(r)*r^2 * 4pi } = 1

         Profiles:
                 htype = 0:
                     Q(r) = (pressure-pressure[nr-1]) * fil(r),
                            with fil(r) = 1/2 * (tanh( [r-r0]/width + 1), 
                            where r0 = hpars[0]  and width = hpars[1]. 
                            
                            If pressure is left undefined, it is assumed
                            to have a constant value of unity.
                 htype > 0:
                        Currently there are no other types defined..."""
    if (htype == 0):
        nr= len(radius)
        r0 = hpars[0]
        width = hpars[1]
        # Check to see if pressure was provided.
        # If not, assume it is constant
        try:
            plen = len(pressure)
        except:
            pressure = 0*radius+1.0e0
        
        x = (radius-r0)/width

        fil= (numpy.tanh(x)+1)/2  # Filter function used to truncate heating at base of CZ
        profile = fil*(pressure-pressure[nr-1])  # Heating profile
        profile = profile/numpy.max(profile)
        
        ###############################################################
        # Next, we need to integrate the heating profile
        # We normalize it such that it's integral over the volume is 1
        # This way, we can set the luminosity via a constant in the input file

        qint = 0
        lq = numpy.zeros(nr)
        integrand= numpy.pi*4*radius*radius*profile
        # First pass, compute integral to normalize
        for i in range(1,nr):
            lq[i] = numpy.trapz(integrand[0:i+1],x=radius[0:i+1])

        profile = profile/lq[nr-1]
        return profile
    
def gen_poly(radius,n,nrho,mass,rhoi,gconst,cp,rb):
    """
    Returns a thermodynamic background state defined by a polytrope
    with polytropic index {n} and {nrho} density scaleheights between
    {radius[nr-1]} and {rb}, where nr-1 is the number of gridpoints used
    to define the radial grid.  Gravity is assumed to follow a 1/r^2 profile
    when computing the polytropic profile.
    
    Inputs:
            radius -- the radial grid on which to compute the polytrope
            n      -- the polytropic index
            nrho   -- ln(rho(radius[nr-1])/rhoi)
            mass   -- the mass interior to the polytrope
            rhoi   -- the value of density at r = rb 
            gconst -- G (the graviational constant)
                      When computing the poltropic profile, the gravitational 
                      acceleration g is given by g = gconst*mass/radius^2.
            cp     -- The specific heat at constant pressure
            rb     -- Reference radius used for specifying rhoi and enforcing nrho.
            
    Returns:
            A background_state object with the following attribues defined is returned:
                temperature, density, pressure, entropy,
                radial entropy gradient, radial pressure gradient
    """
    nr = len(radius)
    rt   = radius[nr-1]
    beta = rb/rt
    gamma = 5.0/3.0 # microscopic property of the ideal gas (unrelated to n...)
    gas_constant = cp*(1.0-1.0/gamma)  # R
    C = numpy.exp(nrho/n)
    f = (C*beta-1)/(1-C)
    
    btwiddle = mass*gconst/(gas_constant*(n+1))  # b, modulo T_0
    temperature = btwiddle*(1.0/radius+f/rb)
    ti = btwiddle*(1.0/rb+f/rb)  # temperature at inner point

    z = temperature/ti
    b = btwiddle/ti
    dzdr = -b/radius/radius

    zi = btwiddle*(1.0/rb+f/rb)/ti

    rho_0 = rhoi/(zi**n)
    density = rho_0*(z**n)

    Pi = gas_constant*ti*rhoi
    
    P_0 = Pi/(zi**(n+1))
    
    pressure = gas_constant*temperature*density
    
    dpdr = P_0*(n+1)*(z**n)*dzdr
    
    dsdr = (cp/z)*((n+1)/gamma - n)*dzdr
    entropy = cp*numpy.log((pressure**(1.0/gamma))/density)

    new_poly = background_state(radius, pressure = pressure, temperature = temperature,
                                entropy = entropy, entropy_gradient=dsdr, 
                                pressure_gradient= dpdr, density=density)
    return new_poly

