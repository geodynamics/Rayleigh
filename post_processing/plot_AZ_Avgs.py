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

# VI.  Azimuthal Averages
# =============
# 
# 
# **Summary:**    Azimuthal averages of requested output variables.  Each output variable is stored as a 2-D function of radius and latitude.
# 
# **Subdirectory:**  AZ_Avgs
# 
# **main_input prefix:** azavg
# 
# **Python Class:** AZ_Avgs
# 
# **Additional Namelist Variables:**  
# None
# 
# Azimuthally-Averaged outputs are particularly useful for examining a system's mean flows (i.e., differential rotation and meridional circulation). 
# 
# Examining the *main_input* file, we see that the following output values have been denoted for the Azimuthal Averages (see *rayleigh_output_variables.pdf* for mathematical formulae):
# 
# 
# | Menu Code  | Description |
# |------------|-------------|
# | 1          | Radial Velocity |
# | 2          | Theta Velocity |
# | 3          | Phi Velocity  |
# | 201        | Radial Mass Flux |
# | 202        | Theta Mass Flux |
# | 501        | Temperature Perturbation |
# 
# 
# 
# In the example that follows, we demonstrate how to plot azimuthal averages, including how to generate streamlines of mass flux.   Note that since the benchmark is Boussinesq, our velocity and mass flux fields are identical.  This is not the case when running an anelastic simulation.
# 
# We begin with the usual preamble and also import two helper routines used for displaying azimuthal averages.
# 
# Examining the data structure, we see that the vals array is dimensioned to account for latitudinal variation, and that we have new attributes costheta and sintheta used for referencing locations in the theta direction.

# In[13]:

from rayleigh_diagnostics import AZ_Avgs, build_file_list, plot_azav, streamfunction
import matplotlib.pyplot as plt
import pylab
import numpy
#from azavg_util import *
files = build_file_list(30000,40000,path='AZ_Avgs')
az = AZ_Avgs(files[0],path='')



# ***
# Before creating our plots, let's time-average over the last two files that were output (thus sampling the equilibrated phase).

# In[14]:




nfiles = len(files)
tcount=0
for i in range(nfiles):
    az=AZ_Avgs(files[i],path='')

    if (i == 0):
        nr = az.nr
        ntheta = az.ntheta
        nq = az.nq
        azavg=numpy.zeros((ntheta,nr,nq),dtype='float64')

    for j in range(az.niter):
        azavg[:,:,:] += az.vals[:,:,:,j]
        tcount+=1
azavg = azavg*(1.0/tcount)  # Time steps were uniform for this run, so a simple average will suffice

lut = az.lut
vr = azavg[:,:,lut[1]]
vtheta = azavg[:,:,lut[2]]
vphi = azavg[:,:,lut[3]]
rhovr = azavg[:,:,lut[201]]
rhovtheta = azavg[:,:,lut[202]]
temperature = azavg[:,:,lut[501]]
radius = az.radius
costheta = az.costheta
sintheta = az.sintheta


# Before we render, we need to do some quick post-processing:
# 1. Remove the spherical mean temperature from the azimuthal average.
# 2. Convert v_phi into omega
# 3. Compute the magnitude of the mass flux vector
# 4. Compute stream function associated with the mass flux field

# In[15]:

#Subtrace the ell=0 component from temperature at each radius
for i in range(nr):
    temperature[:,i]=temperature[:,i] - numpy.mean(temperature[:,i])

#Convert v_phi to an Angular velocity
omega=numpy.zeros((ntheta,nr))
for i in range(nr):
    omega[:,i]=vphi[:,i]/(radius[i]*sintheta[:])

#Generate a streamfunction from rhov_r and rhov_theta
psi = streamfunction(rhovr,rhovtheta,radius,costheta,order=0)
#contours of mass flux are overplotted on the streamfunction PSI
rhovm = numpy.sqrt(rhovr**2+rhovtheta**2)*numpy.sign(psi)    


# Finally, we render the azimuthal averages.  
# **NOTE:**  If you want to save any of these figures, you can mimic the saveplot logic at the bottom of this example.

# In[16]:

#   We do a single row of 3 images 
#   Spacing is default spacing set up by subplot
figdpi=300
sizetuple=(5.5*3,3*3)


tsize = 20     # title font size
cbfsize = 10   # colorbar font size
fig, ax = plt.subplots(ncols=3,figsize=sizetuple,dpi=figdpi)
plt.rcParams.update({'font.size': 14})

#temperature
#ax1 = f1.add_subplot(1,3,1)
units = '(nondimensional)'
plot_azav(fig,ax[0],temperature,radius,costheta,sintheta,mycmap='RdYlBu_r',boundsfactor = 2, 
          boundstype='rms', units=units, fontsize = cbfsize)
ax[0].set_title('Temperature',fontsize=tsize)

#Differential Rotation
#ax1 = f1.add_subplot(1,3,2)
units = '(nondimensional)'
plot_azav(fig,ax[1],omega,radius,costheta,sintheta,mycmap='RdYlBu_r',boundsfactor = 1.5, 
          boundstype='rms', units=units, fontsize = cbfsize)
ax[1].set_title(r'$\omega$',fontsize=tsize)

#Mass Flux
#ax1 = f1.add_subplot(1,3,3)
units = '(nondimensional)'
plot_azav(fig,ax[2],psi,radius,costheta,sintheta,mycmap='RdYlBu_r',boundsfactor = 1.5, 
          boundstype='rms', units=units, fontsize = cbfsize, underlay = rhovm)
ax[2].set_title('Mass Flux',fontsize = tsize)

savefile = 'AZ_Avgs.pdf'
print('Saving figure to: ', savefile)
plt.savefig(savefile)

