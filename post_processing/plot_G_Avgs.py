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

# IV.  Global Averages
# ------------
# 
# **Summary:**   Full-volume averages of requested output variables over the full, spherical shell
# 
# **Subdirectory:**  G_Avgs
# 
# **main_input prefix:** globalavg
# 
# **Python Class:** G_Avgs
# 
# **Additional Namelist Variables:**  
# None
# 
# *Before proceeding, ensure that you have copied Rayleigh/etc/analysis/rayleigh_diagnostics.py to your simulation directory.  This Python module is required for reading Rayleigh output into Python.*
# 
# Examining the *main_input* file, we see that the following output values have been denoted for the Global Averages (see *rayleigh_output_variables.pdf* for the mathematical formulae):
# 
# 
# | Menu Code | Description |
# |-----------|-------------|
# | 401       | Full Kinetic Energy Density (KE) |
# | 402       | KE (radial motion) |
# | 403       | KE (theta motion)  |
# | 404       | KE (phi motion) |
# | 405       | Mean Kinetic Energy Density (MKE) |
# | 406       | MKE (radial motion) |
# | 407       | MKE (theta motion) |
# | 408       | MKE (phi motion) |
# | 409       | Fluctuating Kinetic Energy Density (FKE) |
# | 410       | FKE (radial motion) |
# | 411       | FKE (theta motion) |
# | 412       | FKE (phi motion) |
# 
# In the example that follows, we will plot the time-evolution of these different contributions to the kinetic energy budget.  We begin with the following preamble:

# In[1]:

from rayleigh_diagnostics import G_Avgs, build_file_list
import matplotlib.pyplot as plt
import numpy


# The preamble for each plotting example will look similar to that above.  We import the numpy and matplotlib.pyplot modules, aliasing the latter to *plt*.   We also import two items from *rayleigh_diagnostics*: a helper function *build_file_list* and the *GlobalAverage* class. 
# 
# The *G_Avgs* class is the Python class that corresponds to the full-volume averages stored in the *G_Avgs* subdirectory of each Rayleigh run.
# 
# We will use the build_file_list function in many of the examples that follow.  It's useful when processing a time series of data, as opposed to a single snapshot.  This function accepts three parameters: a beginning time step, an ending time step, and a subdirectory (path).  It returns a list of all files found in that directory that lie within the inclusive range [beginning time step, ending time step].  The file names are prepended with the subdirectory name, as shown below.

# In[2]:

# Build a list of all files ranging from iteration 0 million to 1 million
files = build_file_list(0,1000000,path='G_Avgs')
print(files)


# We can create an instance of the G_Avgs class by initializing it with a filename.  The optional keyword parameter *path* is used to specify the directory.  If *path* is not specified, its value will default to the subdirectory name associated with the datastructure (*G_Avgs* in this instance).  
# 
# Each class was programmed with a **docstring** describing the class attributes.   Once you created an instance of a rayleigh_diagnostics class, you can view its attributes using the help function as shown below.

# In[3]:

a = G_Avgs(filename=files[0],path='')  # Here, files[0]='G_Avgs/00010000'
#a= G_Avgs(filename='00010000') would yield an equivalent result

#help(a)


# Examining the docstring, we see a few important attributes that are common to the other outputs discussed in this document:
# 1.  niter -- the number of time steps in the file
# 2.  nq   -- the number of output variables stored in the file
# 3.  qv   -- the menu codes for those variables
# 4.  vals -- the actual data
# 5.  time -- the simulation time corresponding to each output dump
# 
# The first step in plotting a time series is to collate the data.

# In[4]:

# Loop over all files and concatenate their data into a single array
nfiles = len(files)
for i,f in enumerate(files):
    a = G_Avgs(filename=f,path='')
    if (i == 0):
        nq = a.nq
        niter = a.niter
        gavgs = numpy.zeros((niter*nfiles,nq),dtype='float64')
        iters = numpy.zeros(niter*nfiles,dtype='int32')
        time = numpy.zeros(niter*nfiles,dtype='float64')
    i0 = i*niter
    i1 = (i+1)*niter
    gavgs[i0:i1,:] = a.vals
    time[i0:i1] = a.time
    iters[i0:i1] = a.iters



# The Lookup Table (LUT)
# ------------------
# 
# The next step in the process is to identify where within the *gavgs* array our deisired output variables reside.  Every Rayleigh file object possesses a lookup table (lut).  The lookup table is a python list used to identify the index within the vals array where a particular menu code resides.  For instance, the menu code for the theta component of the velocity is 2.  The location of v_theta in the vals array is then stored in lut[2].  
# 
# Note that you should never assume that output variables are stored in any particular order.  Moreover, the lookup table is unique to each file and is likely to change during a run if you modify the output variables in between restarts.  When running the benchmark, we kept a consistent set of outputs throughout the entirety of the run.  This means that the lookup table did not change between outputs and that we can safely use the final file's lookup table (or any other file's table) to reference our data.
# 
# Plotting Kinetic Energy
# ---------------------------
# Let's examine the different contributions to the kinetic energy density in our models.  Before we can plot, we should use the lookup table to identify the location of each quantity we are interested in plotting.

# In[5]:

#The indices associated with our various outputs are stored in a lookup table
#as part of the GlobalAverage data structure.  We define several variables to
#hold those indices here:

lut = a.lut
ke  = lut[401]  # Kinetic Energy (KE)
rke = lut[402]  # KE associated with radial motion
tke = lut[403]  # KE associated with theta motion
pke = lut[404]  # KE associated with azimuthal motion

#We also grab some energies associated with the mean (m=0) motions
mke  = lut[405]
mrke = lut[406]  # KE associated with mean radial motion
mtke = lut[407]  # KE associated with mean theta motion
mpke = lut[408]  # KE associated with mean azimuthal motion

#We also output energies associated with the fluctuating/nonaxisymmetric
#motions (e.g., v- v_{m=0})
fke  = lut[409]
frke = lut[410]  # KE associated with mean radial motion
ftke = lut[411]  #KE associated with mean theta motion
fpke = lut[412]  # KE associated with mean azimuthal motion


# To begin with, let's plot the total, mean, and fluctuating kinetic energy density during the initial transient phase, and then during the equilibrated phase.

# In[6]:

sizetuple=(10,3)
fig, ax = plt.subplots(ncols=2, figsize=sizetuple)
ax[0].plot(time, gavgs[:,ke], label='KE')
ax[0].plot(time, gavgs[:,mke],label='MKE')
ax[0].plot(time, gavgs[:,fke], label='FKE')
ax[0].legend(loc='center right', shadow=True)
ax[0].set_xlim([0,0.2])
ax[0].set_title('Equilibration Phase')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Energy')

ax[1].plot(time, gavgs[:,ke], label='KE')
ax[1].plot(time, gavgs[:,mke], label = 'MKE')
ax[1].plot(time,gavgs[:,fke],label='FKE')
ax[1].legend(loc='center right', shadow=True)
ax[1].set_title('Entire Time-Trace')
ax[1].set_xlabel('Time')
ax[1].set_ylabel('Energy')

saveplot = True # Plots appear in the notebook and are not written to disk (set to True to save to disk)
savefile = 'energy_trace.pdf'  #Change .pdf to .png if pdf conversion gives issues
plt.tight_layout()

print('Saving figure to: ', savefile)
plt.savefig(savefile)


# We can also look at the energy associated with each velocity component.
# Note that we log scale in the last plot.  There is very little mean radial or theta kinetic energy; it is mostly phi energy.

# In[7]:

sizetuple=(5,10)
xlims=[0,0.2]
fig, ax = plt.subplots(ncols=1, nrows=3, figsize=sizetuple)
ax[0].plot(time, gavgs[:,ke], label='KE')
ax[0].plot(time, gavgs[:,rke],label='RKE')
ax[0].plot(time, gavgs[:,tke], label='TKE')
ax[0].plot(time, gavgs[:,pke], label='TKE')
ax[0].legend(loc='center right', shadow=True)
ax[0].set_xlim(xlims)
ax[0].set_title('Total KE Breakdown')
ax[0].set_xlabel('Time')
ax[0].set_ylabel('Energy')

ax[1].plot(time, gavgs[:,fke],  label='FKE')
ax[1].plot(time, gavgs[:,frke], label='FRKE')
ax[1].plot(time, gavgs[:,ftke], label='FTKE')
ax[1].plot(time, gavgs[:,fpke], label='FPKE')
ax[1].legend(loc='center right', shadow=True)
ax[1].set_xlim(xlims)
ax[1].set_title('Fluctuating KE Breakdown')
ax[1].set_xlabel('Time')
ax[1].set_ylabel('Energy')

ax[2].plot(time, gavgs[:,mke],  label='MKE')
ax[2].plot(time, gavgs[:,mrke], label='MRKE')
ax[2].plot(time, gavgs[:,mtke], label='MTKE')
ax[2].plot(time, gavgs[:,mpke], label='MPKE')
ax[2].legend(loc='lower right', shadow=True)
ax[2].set_xlim(xlims)
ax[2].set_title('Mean KE Breakdown')
ax[2].set_xlabel('Time')
ax[2].set_ylabel('Energy')
ax[2].set_yscale('log')

plt.tight_layout()
savefile = 'energy_breakdown.pdf'
print('Saving figure to: ', savefile)
plt.savefig(savefile)


