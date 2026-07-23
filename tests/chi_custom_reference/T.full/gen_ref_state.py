import reference_tools as rt
import numpy
# Grid Parameters
nr    = 192     # Number of radial points - this is intentionally not the same as main_input as I want to test the interpolation of the reference state

#aspect ratio
ar = 0.35e0
# shell depth 
d=1.e0 

#outer radial boundary
ro=d/(1-ar)

#inner radial boundary
ri=ar*ro

#radial grid
radius=numpy.linspace(ri,ro,nr)

ones = numpy.ones(nr,dtype='float64')
Q0 = 10.0

# Initialize an equation_coefficients structure consistent with the case 0 benchmark test
my_ref = rt.equation_coefficients(radius)

my_ref.set_function(ones,'heating')  # heating function
my_ref.set_constant(Q0,'luminosity')  # multiplies the Heating function
my_ref.write('ref_state.dat') # Write the data file

chi_ref = rt.scalar_equation_coefficients(1, radius)

chi_ref.set_function(ones,'source_chi')
chi_ref.set_constant(Q0, 'source_chi_scale')
chi_ref.write('chi_ref_state.dat')


