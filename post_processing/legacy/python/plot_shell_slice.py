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

####################################################################################################
#
#  Shell-Slice (Shell_Slices) plotting example
#  - Reads in a single Shell_Slice file.
#  - Plots vr, vphi, and entropy
#
#  This example routine makes use of the ShellSlice
#  data structure associated with the Shell_Slices output.
#  Upon initializing a ShellSlice object, the 
#  object will contain the following attributes:
#    ----------------------------------
#    self.niter                                    : number of time steps
#    self.nq                                       : number of diagnostic quantities output
#    self.nr                                       : number of shell slices output
#    self.ntheta                                   : number of theta points
#    self.nphi                                     : number of phi points
#    self.qv[0:nq-1]                               : quantity codes for the diagnostics output
#    self.radius[0:nr-1]                           : radii of the shell slices output
#    self.inds[0:nr-1]                             : radial indices of the shell slices output
#    self.costheta[0:ntheta-1]                     : cos(theta grid)
#    self.sintheta[0:ntheta-1]                     : sin(theta grid)
#    self.vals[0:nphi-1,0:ntheta-1,0:nr-1,0:nq-1,0:niter-1] 
#                                                  : The shell slices 
#    self.iters[0:niter-1]                         : The time step numbers stored in this output file
#    self.time[0:niter-1]                          : The simulation time corresponding to each time step
#    self.version                                  : The version code for this particular output (internal use)
#    self.lut                                      : Lookup table for the different diagnostics output
#    -------------------------------------
import pylab as p 
import numpy as np

from diagnostic_reading import ShellSlice
import matplotlib.pyplot as plt
from matplotlib import ticker

# Set saveplot to True to save to a .png file. 
# Set to False to view plots interactively on-screen.
saveplot = False
savefile = 'shell_slice.png'


#Read in our shell slice
a = ShellSlice(filename='03280000',path='Shell_Slices/')

#Identify the variables indices for vr,vphi, and entropy
var_inds = [a.lut[1],a.lut[3],a.lut[64]]

rad_inds = [0,3,6]  # pick 3 depth to plot at
tind = 0 # grab time index 0 (the first record of the file)

#Tex can be enclosed in dollar signs within a string.  The r in front of the string is necessary for strings enclosing Tex
units = [r'm s$^{-1}$', r'm s$^{-1}$', r'erg g$^{-1}$ K$^{-1}$']  
vnames = [r'v$_r$', r'v$_\phi$', "S'"]
ncol = len(var_inds)
nrow = len(rad_inds)
nplots = ncol*nrow
        

#Create the plots.  This first portion only handles the projection and colorbars.
#Labeling comes further down.
ind = 1
f1 = p.figure(figsize=(5.5*3, 5*3), dpi=80)
for  j in range(nrow):
    for i in range(ncol):
        sslice = a.vals[:,:,rad_inds[j],var_inds[i],tind].reshape(a.nphi,a.ntheta)
        if (i == 2):
            sslice = sslice-np.mean(sslice) #subtract ell=0 from entropy
        else:
            sslice = sslice/100.0 # convert the velocity field to m/s
        sslice = np.transpose(sslice)
        
        # Set the projection
        ax1 = f1.add_subplot(nrow,ncol,ind, projection="mollweide")
        
        ind = ind+1


        twosigma = 2*np.std(sslice)  #limits are set to +/- twosigma

        contour_levels = twosigma*np.linspace(-1,1,256)
        image1 = ax1.imshow(sslice,vmin=-twosigma, vmax=twosigma,
                 extent=(-np.pi,np.pi,-np.pi/2,np.pi/2), clip_on=False,
                 aspect=0.5, interpolation='bicubic')
        image1.axes.get_xaxis().set_visible(False)  # Remove the longitude and latitude axis labels
        image1.axes.get_yaxis().set_visible(False)
        image1.set_cmap('RdYlBu_r')  # Red/Blue map


        cbar = f1.colorbar(image1,orientation='horizontal', shrink=0.5)
        cbar.set_label(units[i])

        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()

# Next, set various parameters that control the plot layout
# These parameters can also be set interactively when using plt.show()
pbottom = 0.05
pright = 0.95
pleft = 0.15
ptop = 0.95
phspace = 0.1
pwspace = 0.03

p.subplots_adjust(left = pleft, bottom = pbottom, right = pright, top = ptop, wspace=pwspace,hspace=phspace) 

# Add information about radial depth to the left margin
rspace = (ptop-pbottom)/nrow
rnorm = 6.96e10
lilbit = 0.05
xpos = 0.1*pleft
for i in range(nrow):
    ypos = ptop-rspace*0.4-i*(rspace+phspace/(nrow-1))
    ratio = float(a.radius[rad_inds[i]]/rnorm)
    r_str = r'r/R$_\odot$ = %1.3f' % ratio
    f1.text(xpos,ypos,r_str,fontsize=16)
    ypos = ypos-lilbit
    ind_str = 'radial index: '+str(a.inds[rad_inds[i]])
    f1.text(xpos,ypos,ind_str,fontsize=16)

#Label the plots (entropy, v_r etc.)
cspace = (pright-pleft)/ncol
for i in range(ncol):
    ypos = ptop
    xpos = pleft+cspace*0.47+i*cspace 
    f1.text(xpos,ypos,vnames[i],fontsize=20)

#Save or display the figure
if (saveplot):
    p.savefig(savefile)  
else:
    plt.show()

