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
zeros = numpy.zeros(nr,dtype='float64')




gravity_power=1e0
buoy = (radius/radius[nr-1])**gravity_power # buoyancy term calculation
my_ref.set_function(ones,1)   # denisty rho
my_ref.set_function(buoy,2)   # buoyancy term
my_ref.set_function(ones,3)   # nu(r)
my_ref.set_function(ones,4)   # temperature T
my_ref.set_function(ones,5)   # kappa(r)
my_ref.set_function(zeros,6)  # heating function 

c1 = numpy.zeros(1,dtype='float64')
c2 = numpy.zeros(1,dtype='float64')
c3 = numpy.zeros(1,dtype='float64')
c6 = numpy.zeros(1,dtype='float64')

E  = numpy.zeros(1,dtype='float64')
Ra = numpy.zeros(1,dtype='float64')
Pr = numpy.zeros(1,dtype='float64')

E[0]=1e-3
Ra[0]=1e5
Pr[0]=1

c1 = 2/E
c2 = Ra/Pr
c3 = 1/E
c6 = 1/Pr

my_ref.set_constant(c1[0],1)  # multiplies the Coriolis term
my_ref.set_constant(c2[0],2)  # multiplies the buoyancy
my_ref.set_constant(c3[0],3)  # multiplies the pressure gradient
my_ref.set_constant(1 ,5)     # multiplies the viscosity
my_ref.set_constant(c6[0],6)  # multiplies kappa, here it is: E/Pr
my_ref.set_constant(0,10)     # multiplies the heating, here it is 0, since we have assumed that there is no heating function!
my_ref.write('case0.dat')     # Here we write our data file to be used to run our simulation


