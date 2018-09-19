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

# VII.2  Meridional Slices
# --------------------------
# 
# **Summary:**    2-D profiles of selected output variables sampled in meridional planes. 
# 
# **Subdirectory:**  Meridional_Slices
# 
# **main_input prefix:** meridional
# 
# **Python Class:** Meridional_Slices
# 
# **Additional Namelist Variables:**  
# 
# * meridional_indices (indicial) : indices along longitudinal grid at which to output meridional planes.
# 
# * meridional_indices_nrm (normalized) : normalized longitudinal grid coordinates at which to output
# 
# 
# The meridional-slice output type allows us to examine how the fluid properties vary in latitude and radius.
# 
# Examining the *main_input* file, we see that the following output values have been denoted for the Meridional Slices (see *rayleigh_output_variables.pdf* for mathematical formulae):
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
# In the example that follows, we demonstrate how to create a 2-D plot of radial velocity in a meridional plane.  The procedure is similar to that used to plot an azimuthal average.
# 
# 
# We begin with the usual preamble and import the *plot_azav* helper function.  Examining the data structure, we see that it is similar to the AZ_Avgs data structure.  The *vals* array possesses an extra dimension relative to its AZ_Avgs counterpart to account for the multiple longitudes that may be output, we see attributes *phi* and *phi_indices* have been added to reference the longitudinal grid.  

# In[19]:

#####################################
#  Meridional Slice
from rayleigh_diagnostics import Meridional_Slices, plot_azav
import numpy
import matplotlib.pyplot as plt
from matplotlib import ticker, font_manager
# Read the data

istring = '00040000'
ms = Meridional_Slices(istring)
tindex =1 # All example quantities were output with same cadence.  Grab second time-index from all.



# In[20]:


radius = ms.radius
costheta = ms.costheta
sintheta = ms.sintheta
phi_index = 0  # We only output one Meridional Slice
vr_ms = ms.vals[phi_index,:,:,ms.lut[1],tindex]
units = 'nondimensional'

# Plot
sizetuple=(8,5)
fig, ax = plt.subplots(figsize=(8,8))
tsize = 20     # title font size
cbfsize = 10   # colorbar font size
ax.axis('equal')  # Ensure that x & y axis ranges have a 1:1 aspect ratio
ax.axis('off')    # Do not plot x & y axes
plot_azav(fig,ax,vr_ms,radius,costheta,sintheta,mycmap='RdYlBu_r',boundsfactor = 4.5, 
          boundstype='rms', units=units, fontsize = cbfsize)
ax.set_title('Radial Velocity',fontsize=tsize)
plt.tight_layout()
savefile = 'Meridional_Slices.pdf'
print('Saving figure to: ', savefile)
plt.savefig(savefile)
