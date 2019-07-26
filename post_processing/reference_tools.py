#Reference State Generation Tools  (reference_tools.py)
import numpy

class equation_coefficients:
    """ equation coeff class  """
    nconst = 10
    nfunc = 14
    version=1


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
            
    def set_function(self,y,fi):
        self.functions[fi-1,:] = y[:]
        self.fset[fi-1] = 1
    def set_constant(self,c,ci):
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





                            

                    

    
     
