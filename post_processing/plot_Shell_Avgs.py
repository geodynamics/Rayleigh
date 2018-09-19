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

# V.  Shell Averages
# ==========
# 
# **Summary:**   Spherical averages of requested output variables.  Each output variable is stored as a 1-D function of radius.
# 
# **Subdirectory:**  Shell_Avgs
# 
# **main_input prefix:** shellavg
# 
# **Python Class:** Shell_Avgs
# 
# **Additional Namelist Variables:**  
# None
# 
# The Shell-Averaged outputs are useful for examining how quantities vary as a function of radius.  They are particularly useful for examining the distribution of energy as a function of radius, or the heat flux balance established by the system.
# 
# Examining the *main_input* file, we see that the following output values have been denoted for the Shell Averages (see *rayleigh_output_variables.pdf* for mathematical formulae):
# 
# 
# | Menu Code  | Description |
# |------------|-------------|
# | 1          | Radial Velocity |
# | 2          | Theta Velocity |
# | 3          | Phi Velocity  |
# | 501        | Temperature Perturbation |
# | 1438       | Radial Convective Heat Flux|
# | 1468       |Radial Conductive Heat Flux |
# 
# 
# In the example that follows, we will plot the spherically-averaged velocity field as a function of radius, the mean temperature profile, and the radial heat flux.  We begin with a preamble similar to that used for the Global Averages.  Using the help function, we see that the Shell_Avgs data structure is similar to that of the G_Avgs.  There are three important differences:
# *  There is a radius attribute (necessary if we want to plot anything vs. radius)
# *  The dimensionality of the values array has changed;  radial index forms the first dimension.
# *  The second dimension of the values array has a length of 4.  In addition to the spherical mean, the 1st, 2nd and 3rd moments are stored in indices 0,1,2, and 3 respectively.

# In[8]:

from rayleigh_diagnostics import Shell_Avgs, build_file_list
import matplotlib.pyplot as plt
import numpy

# Build a list of all files ranging from iteration 0 million to 1 million
files = build_file_list(0,1000000,path='Shell_Avgs')
a = Shell_Avgs(filename=files[0], path='')



# 
# ***
# 
# While it can be useful to look at instaneous snapshots of Shell Averages, it's often useful to examine these outputs in a time-averaged sense.    Let's average of all 200 snapshots in the last file that was output.  We could average over data from multiple files, but since the benchmark run achieves a nearly steady state, a single file will do in this case.

# In[9]:

nfiles = len(files)

nr = a.nr
nq = a.nq
nmom = 4
niter = a.niter
radius = a.radius
savg=numpy.zeros((nr,nmom,nq),dtype='float64')
for i in range(niter):
    savg[:,:,:] += a.vals[:,:,:,i]
savg = savg*(1.0/niter)

lut = a.lut
vr = lut[1]        # Radial Velocity
vtheta = lut[2]    # Theta Velocity
vphi = lut[3]      # Phi Velocity
thermal = lut[501] # Temperature


eflux = lut[1438]  # Convective Heat Flux (radial)
cflux = lut[1468]  # Conductive Heat Flux (radial)


# Velocity vs. Radius
# ---------------------
# Next, we plot the mean velocity field, and its first moment, as a function of radius.   Notice that the radial and theta velocity components have a zero spherical mean.  Since we are running an incompressible model, this is a good sign!

# In[10]:

sizetuple = (7,7)
fig, ax = plt.subplots(nrows=2, ncols =1, figsize=sizetuple)

ax[0].plot(radius,savg[:,0,vr],label=r'$v_r$')
ax[0].plot(radius,savg[:,0,vtheta], label=r'$v_\theta$')
ax[0].plot(radius,savg[:,0,vphi], label=r'$v_\phi$')
ax[0].legend(shadow=True,loc='lower right')
ax[0].set_xlabel('Radius')
ax[0].set_ylabel('Velocity')
ax[0].set_title('Spherically-Averaged Velocity Components')

ax[1].plot(radius,savg[:,1,vr],label=r'$v_r$')
ax[1].plot(radius,savg[:,1,vtheta], label=r'$v_\theta$')
ax[1].plot(radius,savg[:,1,vphi], label=r'$v_\phi$')
ax[1].legend(shadow=True,loc='upper left')
ax[1].set_xlabel('Radius')
ax[1].set_ylabel('Velocity')
ax[1].set_title('Velocity Components: First Spherical Moment')


plt.tight_layout()
savefile = 'velocity_variation.pdf'
print('Saving figure to: ', savefile)
plt.savefig(savefile)



# Radial Temperature Profile
# ------------------------------
# We might also look at temperature ...

# In[11]:

fig, ax = plt.subplots()

ax.plot(radius,savg[:,0,thermal],label='Temperature (mean)')
ax.plot(radius,savg[:,1,thermal]*10, label='Temperature (standard dev.)')
ax.legend(shadow=True,loc='upper right')
ax.set_xlabel('Radius')
ax.set_ylabel('Velocity')
ax.set_title('Temperature')


savefile = 'temperature_variation.pdf'
print('Saving figure to: ', savefile)
plt.savefig(savefile)


# Heat Flux Contributions
# --------------------------
# We can also examine the balance between convective and conductive heat flux.  In this case, before plotting these quantities as a function of radius, we normalize them by the surface area of the sphere to form a luminosity.

# In[12]:

fpr=4.0*numpy.pi*radius*radius
elum = savg[:,0,eflux]*fpr
clum = savg[:,0,cflux]*fpr
tlum = elum+clum
fig, ax = plt.subplots()
ax.plot(radius,elum,label='Convection')
ax.plot(radius,clum, label='Conduction')
ax.plot(radius,tlum, label='Total')
ax.set_title('Flux Balance')
ax.set_ylabel(r'Energy Flux ($\times 4\pi r^2$)')
ax.set_xlabel('Radius')
ax.legend(shadow=True)
plt.tight_layout()
savefile = 'heat_flux.pdf'
print('Saving figure to: ', savefile)
plt.savefig(savefile)

