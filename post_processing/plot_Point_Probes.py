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

# IX.  Point Probes
# ============
# 
# **Summary:**    Point-wise sampling of desired output variables.
# 
# **Subdirectory:**  Point_Probes
# 
# **main_input prefix:** point_probe
# 
# **Python Class:**  Point_Probes 
# 
# 
# **Additional Namelist Variables:**  
# 
# * point_probe_r  : radial indices for point-probe output
# 
# * point_probe_theta  : theta indices for point-probe output
# 
# * point_probe_phi : phi indices for point-probe output
# 
# * point_probe_r_nrm  : normalized radial coordinates for point-probe output
# 
# * point_probe_theta_nrm : normalized theta coordinates for point-probe output
# 
# * point_probe_phi_nrm : normalized phi coordinates for point-probe output
# 
# * point_probe_cache_size : number of time-samples to save before accessing the disk 
# 
# 
# Point-probes allow us to sample a simulation at an arbitrary set of points.  This output type serves two purposes:
# 1.  It provides an analog to laboratory measurements where slicing and averaging are difficult, but taking high-time-cadence using (for example) thermistors is common-practice.
# 2.  It provides an alternative method of slicing a model ( for when equatorial, meridional, or shell slices do yield the desired result).
# 
# IX.1  Specifying Point-Probe Locations
# ---------
# 
# Point-probe locations are indicated by specifying a grid.   The user does not supply a set of ordered coordinates (r,theta,phi).  Instead, the user specifies nodes on the grid using the namelist variables described above.   Examples follow.
# 
# **Example 1:  4-point Coarse Grid **
# 
# point_probe_r_nrm = 0.25, 0.5  
# point_probe_theta_nrm = 0.5  
# point_probe_phi_nrm = 0.2, 0.8
# 
# This example would produce point probes at the four coordinates  { (0.25, 0.5, 0.2), (0.25, 0.5, 0.8), (0.5, 0.5, 0.2), (0.5,0.5,0.8) }   (r,theta,phi; normalized coordinates).
# 
# **Example 2:  "Ring" in Phi **
# 
# point_probe_r_nrm = 0.5  
# point_probe_theta_nrm = 0.5  
# point_probe_phi_nrm = 0.0, -1.0
# 
# This example describes a ring in longitude, sampled at mid-shell, in the equatorial plane.  We have made use of the positional range feature here by indicating normalized phi coordinates of 0.0, -1.0.  Rayleigh intreprets this as an instruction to sample all phi coordinates.
# 
# ** Example 3: 2-D Surface in (r,phi)  **
# 
# point_probe_r_nrm = 0, -1.0  
# point_probe_theta_nrm = 0.25  
# point_probe_phi_nrm = 0, -1.0
# 
# This example uses the positional range feature along with normalized coordinates to generate a 2-D slice in r-phi at theta = 45 degrees (theta_nrm = 0.25).  Using the syntax 0,-1.0 instructs *Rayleigh* to grab all r and phi coordinates.
# 
# ** Example 4:  3-D Meridional "Wedges" **
# 
# point_probe_r_nrm = 0.0, -1.0  
# point_probe_theta_nrm =  0.0, -1.0  
# point_probe_phi_nrm =  0.20, -0.30, 0.7, -0.8
# 
# This example generates two 3-D wedges described by all r,theta points and all longitudes in the ranges [72 deg, 108 deg]  and [252 deg, 288 deg].
# 
# IX.2  Point-Probe Caching
# -----------------------------
# 
# When performing sparse spatial sampling using point-probes, it may be desireable to output with a high-time cadence.  As this may cause disk-access patterns characterized by frequent, small writes, the point-probes are programmed with a caching feature.  This feature is activated by specifing the **point_probe_cache_size** variable in the output namelist.
# 
# This variable determines how many time-samples are saved in memory before a write is performed.  Its default value is 1, which means that the disk is accessed with a frequency of **point_probe_frequency**.  If the cache size is set to 10 (say), then samples are still peformed at **point_probe_frequency** but they are only written to disk after 10 have been collected in memory.   
# 
# **NOTE:**  Be sure that **point_probe_cache_size** divides evenly into **point_probe_nrec**.
# 
# 
# IX.3  Example:  Force-Balance with Point Probes
# -------------------------------------------------------
# 
# Our example input file specifies a coarse, six-point grid.  Examining the *main_input* file, we see that all variables necessary to examine the force balance in each direction have been specified. (see *rayleigh_output_variables.pdf* for mathematical formulae):
# 
# 
# | Menu Code  | Description |
# |------------|-------------|
# | 1          | Radial Velocity |
# | 2          | Theta Velocity |
# | 3          | Phi Velocity  |
# | 1201       | Radial Advection (v dot grad v) |
# | 1202       | Theta Advection |
# | 1203       | Phi Advection  |
# | 1216          | Buoyancy Force (ell=0 component subtracted) |
# | 1219  | Radial Coriolis Force |
# | 1220  | Theta Coriolis Force  |
# | 1221  | Phi Coriolis Force |
# | 1228  | Radial Viscous Force |
# | 1229  | Theta Viscous Force  |
# | 1230  | Phi Viscous Force|
# 
# 
# **Note that the pressure force appears to be missing.** This is not an oversight.  The diagnostic nature of the Pressure equation in incompressible/anelastic models, coupled with the second-order Crank-Nicolson time-stepping scheme, means that the pressure field can exhibit an even/odd sawtoothing in time.  The *effective* pressure force (as implemented through the Crank-Nicolson scheme) is always a weighted average over two time steps **and is always well-resolved in time**.  
# 
# When sampling at regular intervals as we have here, if we directly sample the pressure force, we will sample either the high or low end of the sawtooth envelope, and the force balance will be off by a large factor.  The easiest fix is to output the velocity field and compute its time derivative.  This, in tandem with the sum of all other forces, can be used to calculate the effective pressure as a post-processing step.  The (undesireable)  alternative is to output once every time step and compute the effective pressure using the Crank-Nicolson weighting.    
# 
# We demonstrate how to compute the effective pressure force via post-processing in the example below.

# In[25]:

from rayleigh_diagnostics import Point_Probes,  build_file_list
import numpy
from matplotlib import pyplot as plt

#Decide which direction you want to look at (set direction = {radial,theta, or phi})
#This is used to determine the correct quantity codes below
radial = 0
theta  = 1
phi    = 2
direction=radial
# Build a list of all files ranging from iteration 0 million to 1 million
files = build_file_list(0,1000000,path='Point_Probes')
nfiles = len(files)-1


for i in range(nfiles):
    pp = Point_Probes(files[i],path='')
    if (i == 0):
        nphi = pp.nphi
        ntheta = pp.ntheta
        nr = pp.nr
        nq = pp.nq
        niter = pp.niter
        vals=numpy.zeros( (nphi,ntheta,nr,nq,niter*nfiles),dtype='float64')
        time=numpy.zeros(niter*nfiles,dtype='float64')
    vals[:,:,:,:, i*niter:(i+1)*niter] = pp.vals
    time[i*niter:(i+1)*niter]=pp.time
istring='00040000' # iteration to examine

##################################################
# We choose the coordinate indices **within**
# the Point-Probe array that we want to examine
# These indices start at zero and run to n_i-1
# where n_i is the number of points sampled in 
# the ith direction

# Use help(pp) after loading the Point-Probe file
# to see the Point-Probe class structure

pind = 0           # phi-index to examine
rind = 0           # r-index to examine
tind = 0           # theta-index to examine


pp   = Point_Probes(istring)
lut  = pp.lut

nt   = pp.niter


#######################################################################
#  Grab velocity from the point probe data
u  = vals[pind,0,rind,pp.lut[1+direction],:]
dt=time[1]-time[0]


###########################################################################
# Use numpy to compute time-derivative of u
# (necessary to compute a smooth effective pressure without outputing every timestep)

#Depending on Numpy version, gradient function takes either time (array) or dt (scalar)
try:
    dudt = numpy.gradient(u,time)
except:
    dt = time[1]-time[0]  # Assumed to be constant...
    dudt = numpy.gradient(u,dt) 



################################################################
#      Forces (modulo pressure)
# Note the minus sign for advection.  Advective terms are output as u dot grad u, not -u dot grad u
advec = -vals[ pind, tind, rind, lut[1201 + direction], :] 
cor   = vals[ pind, tind, rind, lut[1219 + direction], :]
visc  = vals[ pind, tind, rind, lut[1228 + direction], :]
forces = visc+cor+advec
if (direction == radial):
    buoy  = vals[ pind, tind, rind, lut[1216], :]
    forces = forces+buoy


############################################3
# Construct effective pressure force
pres = dudt-forces
forces = forces+pres
############################################################
# Set up the plot
yfsize='xx-large'  # size of y-axis label

ustrings = [r'u_r', r'u_\theta', r'u_\phi']
ustring=ustrings[direction]
dstring = r'$\frac{\partial '+ustring+'}{\partial t}$'
fstrings = [r'$\Sigma\,F_r$' , r'$\Sigma\,F_\theta$' , r'$\Sigma\,F_\phi$' ]
fstring = fstrings[direction]
diff_string = dstring+' - '+fstring

pstring = 'pressure'
cstring = 'coriolis'
vstring = 'viscous'
bstring = 'buoyancy'
fig, axes = plt.subplots(nrows=2, figsize=(7*2.54, 9.6))
ax0 = axes[0]
ax1 = axes[1]


########################################
# Upper: dur/dt and F_total
#mpl.rc('xtick', labelsize=20) --- still trying to understand xtick label size etc.
#mpl.rc('ytick', labelsize=20)

ax0.plot(time,forces, label = fstring)
ax0.plot(time,pres,label=pstring)
ax0.plot(time,cor,label=cstring)
ax0.plot(time,visc,label=vstring)
if (direction == radial):
    ax0.plot(time,buoy,label=bstring)
ax0.set_xlabel('Time', size=yfsize)

ax0.set_ylabel('Acceleration', size=yfsize)
ax0.set_title('Equilibration Phase',size=yfsize)
ax0.set_xlim([0,0.1])
leg0 = ax0.legend(loc='upper right', shadow=True, ncol = 1, fontsize=yfsize) 

##########################################
# Lower:  Numpy Gradient Approach
ax1.plot(time,forces,label=fstring)
ax1.plot(time,pres,label=pstring)
ax1.plot(time,cor,label=cstring)
ax1.plot(time,visc,label=vstring)
if (direction == radial):
    ax1.plot(time,buoy,label=bstring)
ax1.set_title('Late Evolution',size=yfsize)
ax1.set_xlabel('Time',size=yfsize)
ax1.set_ylabel('Acceleration', size =yfsize)
ax1.set_xlim([0.2,4])
leg1 = ax1.legend(loc='upper right', shadow=True, ncol = 1, fontsize=yfsize)


plt.tight_layout()
savefile = 'Point_Probes.pdf'
print('Saving figure to: ', savefile)
plt.savefig(savefile)

