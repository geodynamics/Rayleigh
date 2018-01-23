# VII.3  Shell Slices
# --------------------------
# 
# **Summary:**    2-D, spherical profiles of selected output variables sampled in at discrete radii. 
# 
# **Subdirectory:**  Shell_Slices
# 
# **main_input prefix:** shellslice
# 
# **Python Class:** Shell_Slices
# 
# **Additional Namelist Variables:**  
# 
# * shellslice_levels (indicial) : indices along radial grid at which to output spherical surfaces.
# 
# * shellslice_levels_nrm (normalized) : normalized radial grid coordinates at which to output spherical surfaces.
# 
# 
# The shell-slice output type allows us to examine how the fluid properties vary on spherical surfaces.
# 
# Examining the *main_input* file, we see that the following output values have been denoted for the Shell Slices (see *rayleigh_output_variables.pdf* for mathematical formulae):
# 
# 
# | Menu Code  | Description |
# |------------|-------------|
# | 1          | Radial Velocity |
# | 2          | Theta Velocity |
# | 3          | Phi Velocity  |
# 
# 
# 
# 
# In the examples that follow, we demonstrate how to create a 2-D plot of radial velocity:
# * on a Cartesian, lat-lon grid
# * projected onto a spherical surface using Basemap   (MUST set use Basemap=True below)
# 
# 
# 
# Plotting on a lat-lon grid is straightforward and illustrated below.   The shell-slice data structure is also displayed via the help() function in the example below and contains information needed to define the spherical grid for plotting purposes.
# 
# It is worth noting the *slice_spec* keyword (described in the docstring) that can be passed to the init method.  When reading large shell slices, a user can save time and memory during the read process by specifying the slice they want to read.

# In[21]:

#####################################
#  Shell Slice
from rayleigh_diagnostics import Shell_Slices
import numpy
import matplotlib.pyplot as plt
from matplotlib import ticker, font_manager
# Read the data

istring = '00040000'
ss = Shell_Slices(istring)

ntheta = ss.ntheta
nphi = ss.nphi
costheta = ss.costheta
theta = numpy.arccos(costheta)


tindex =1 # All example quantities were output with same cadence.  Grab second time-index from all.
rindex = 0 # only output one radius
sizetuple=(8,8)

vr = ss.vals[:,:,rindex,ss.lut[1],tindex]
fig, ax = plt.subplots(figsize=sizetuple)


img = plt.imshow(numpy.transpose(vr), extent=[0,360,-90,90])
ax.set_xlabel( 'Longitude')
ax.set_ylabel( 'Latitude')
ax.set_title(  'Radial Velocity')

plt.tight_layout()
savefile = 'Shell_Slices_LatLon.pdf'
print('Saving figure to: ', savefile)
plt.savefig(savefile)


