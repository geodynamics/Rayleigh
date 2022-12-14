import reference_tools as rt
import numpy
# Grid Parameters
nr    = 128     # Number of radial points

#aspect ratio
#beta=ri/ro 
beta = 0.35e0


# shell depth 
d=1.e0 

#outer radial boundary
ro=d/(1-beta)

#inner radial boundary
ri=beta*ro

#radial grid
radius=numpy.linspace(ri,ro,nr)

# Initialize an equation_coefficients structure consistent with the case 0 benchmark test
my_ref = rt.equation_coefficients(radius)

ones = numpy.ones(nr,dtype='float64')


gravity_power=1e0
buoy = (radius/radius[nr-1])**gravity_power # buoyancy term calculation
my_ref.set_function(buoy,2)   # buoyancy term
my_ref.set_function(ones,6)  # heating function 

c10 = numpy.zeros(1,dtype='float64')
c2 = numpy.zeros(1,dtype='float64')


Qamp  = numpy.zeros(1,dtype='float64')
Ra = numpy.zeros(1,dtype='float64')
Pr = numpy.zeros(1,dtype='float64')

Qamp[0]=1e1
Ra[0]=1.2e5
Pr[0]=1

c10 = Qamp
c2 = Ra/Pr

my_ref.set_constant(c10[0],10)  # multiplies the Heating function
my_ref.set_constant(c2[0],2)    # multiplies the buoyancy
my_ref.write('with_custom.dat') # Write the data file


