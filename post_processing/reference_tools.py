#Reference State Generation Tools  (reference_tools.py)
import numpy

class equation_coefficients:
    """ equation coeff class  """
    nconst = 10
    nfunc = 16
    nr = 0
    functions = numpy.zeros((nfunc,1)  , dtype='float64' )
    radius    = numpy.zeros(1          , dtype='float64' )
    constants = numpy.zeros(nconst     , dtype='float64' )
    cset      = numpy.zeros(nconst     , dtype='int32'   )
    fset      = numpy.zeros(nfunc      , dtype='int32'   )
    def __init__(self,radius):
        nr = len(radius)
        self.nr = nr
        self.radius = numpy.zeros(nr,dtype='float64')
        self.radius[:] = radius[:]
        self.functions  = numpy.zeros((self.nfunc,nr) , dtype='float64' )
    def set_function(self,y,fi):
        self.functions[fi-1,:] = y[:]
        self.fset[fi-1] = 1
    def set_constant(self,c,ci):
        self.constants[ci-1] = c
        self.cset[ci-1] = 1

    def write(self, filename='ecoefs.dat'):
        pi = numpy.array([314],dtype='int32')
        nr = numpy.array([self.nr],dtype='int32')
        fd = open(filename,'wb')
        pi.tofile(fd)
        self.cset.tofile(fd)
        self.fset.tofile(fd)
        self.constants.tofile(fd)
        nr.tofile(fd)
        self.radius.tofile(fd)
        self.functions.tofile(fd)
        fd.close() 





                            

                    

    
     
