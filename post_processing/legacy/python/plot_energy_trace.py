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

###############################################
#
#  Global-Averages (G_Avgs) plotting example
#  Reads in time steps 2.7 million through 4 million
#  Plots average KE vs time
#
#  This example routine makes use of the GlobalAverage
#  data structure associated with the G_Avg output.
#  Upon initializing a GlobalAverage object, the 
#  object will contain the following attributes:
#
#    ----------------------------------
#    self.niter                  : number of time steps
#    self.nq                     : number of diagnostic quantities output
#    self.qv[0:nq-1]             : quantity codes for the diagnostics output
#    self.vals[0:niter-1,0:nq-1] : Globally averaged diagnostics as function of time and quantity index
#    self.iters[0:niter-1]       : The time step numbers stored in this output file
#    self.time[0:niter-1]        : The simulation time corresponding to each time step
#    self.lut                    : Lookup table for the different diagnostics output


from diagnostic_reading import GlobalAverage, build_file_list
import matplotlib.pyplot as plt
import numpy


# Set saveplot to True to save to a .png file. 
# Set to False to view plots interactively on-screen.
saveplot = False
savefile = 'energy_trace.pdf'  #Change .pdf to .png if pdf conversion gives issues

# Build a list of all files ranging from iteration 2.7 million to 4 million
files = build_file_list(2700000,4000000,path='G_Avgs')
#The indices associated with our various outputs are stored in a lookup table
#as part of the GlobalAverage data structure.  We define several variables to
#hold those indices here:
a = GlobalAverage(filename=files[0],path='')  # read a file to get the lookup table
ke_index  = a.lut[125]  # Kinetic Energy (KE)
rke_index = a.lut[126]  # KE associated with radial motion
tke_index = a.lut[127]  # KE associated with theta motion
pke_index = a.lut[128]  # KE associated with azimuthal motion

#We also grab some energies associated with the mean (m=0) motions
mrke_index = a.lut[130]  # KE associated with mean radial motion
mtke_index = a.lut[131]  # KE associated with mean theta motion
mpke_index = a.lut[132]  # KE associated with mean azimuthal motion

#Next, we intialize some empy lists to hold those variables
ke = []
rke = []
tke = []
pke = []

mrke = []
mtke = []
mpke = []

# We also define an empty list to hold the time.  This time
# will be absolute time since the beginning of the run, 
# not since the beginning of the first output file.  
# In this case, the run was dimensional, so the time is recorded 
# in seconds.  We will change that to days as we go.
alldays = []

# Next, we loop over all files and grab the desired data
# from the file, storing it in appropriate lists as we go.
for f in files:
    a = GlobalAverage(filename=f,path='')


    days = a.time/(3600.0*24.0)
    for i,d in enumerate(days):
        alldays.append(d)
        ke.append(a.vals[i,ke_index])
        rke.append(a.vals[i,rke_index])
        tke.append(a.vals[i,tke_index])
        pke.append(a.vals[i,pke_index])

        mrke.append(a.vals[i,mrke_index])
        mtke.append(a.vals[i,mtke_index])
        mpke.append(a.vals[i,mpke_index])


#Create two plots.
# Plot 1:  Total KE and its breakdown.
# Plot 2:  Breakdown of the mean KE

#Changing the plotting parameter somewhat depending on saveplot
if (saveplot):
    plt.figure(1,figsize=(7.5, 4.0), dpi=300)
    plt.rcParams.update({'font.size': 12})
else:
    plt.figure(1,figsize=(15,5),dpi=100)
    plt.rcParams.update({'font.size': 14})

plt.subplot(121)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.plot(alldays,ke,'black', label = 'KE'+r'$_{total}$')
plt.plot(alldays,rke,'r', label = 'KE'+r'$_{rad}$')
plt.plot(alldays,tke,'g', label = 'KE'+r'$_\theta$')
plt.plot(alldays,pke,'b', label = 'KE'+r'$_\phi$')
plt.yscale('log')
plt.xlabel('Time (days)')
plt.ylabel('Energy Density '+r'( erg cm$^{-3}$)')
if (saveplot):
    legend = plt.legend(loc='upper right', shadow=True, ncol = 2, fontsize = 'x-small') 
else:
    legend = plt.legend(loc='upper right', shadow=True, ncol = 2) 

t1 = r'$\frac{1}{2}\rho$'
plt.subplot(122)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.plot(alldays,mrke,'r', label = t1+r'$<v^2_r>$')
plt.plot(alldays,mtke,'g', label = t1+r'$<v^2_\theta>$')
plt.plot(alldays,mpke,'b', label = t1+r'$<v^2_\phi>$')
plt.yscale('log')
plt.ylim([1,1e7])
plt.xlabel('Time (days)')
plt.ylabel('Energy Density '+r'( erg cm$^{-3}$)')
if (saveplot):
    legend = plt.legend(loc='lower right', shadow=True, ncol = 2, fontsize = 'x-small') 
else:
    legend = plt.legend(loc='lower right', shadow=True, ncol = 2) 
plt.tight_layout()

if (saveplot):
    plt.savefig(savefile)
else:
    plt.show()
