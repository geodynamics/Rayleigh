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

###########################################################################
#  This file contains utilities useful for plotting the AZ-Average files

import numpy as np
from pylab import *
import pylab as p 
import matplotlib.pyplot as plt
from matplotlib import ticker, font_manager
def get_lims(arr,boundstype='minmax',boundsfactor=1,themin=True):
    if (themin):
        if (boundstype == 'minmax'):
            val=arr.min()*boundsfactor
        elif (boundstype == 'rms'):
            val = -np.std(arr)*boundsfactor
    else:
        if(boundstype == 'minmax'):
            val=max(abs(arr.min()),arr.max())*boundsfactor
        elif (boundstype == 'rms'):
            val = np.std(arr)*boundsfactor
    return val
def plot_azav(fig,axis,field,radius,costheta,sintheta,r_bcz=0.71,mini=-1,maxi=-1,mycmap='jet',cbar=True, 
    boundsfactor = 1, boundstype = 'minmax', units = '',fontsize = 12, underlay = [0], nlevs = 6):
    #Modified version of Antoine Strukarek's routine
    #r = radius/6.9599e10
    r = radius/np.max(radius)
    n_r=len(r)
    n_t=len(costheta)
    rtmp = r.reshape(1,n_r)
    cthtmp = costheta.reshape(n_t,1)
    sthtmp = sintheta.reshape(n_t,1)
    xr = p.dot(cthtmp,rtmp)
    yr = p.dot(sthtmp,rtmp)

    if (mini == -1):
        mini = get_lims(field,boundsfactor=boundsfactor,boundstype=boundstype,themin = True)
    if (maxi == -1):
        maxi = get_lims(field,boundsfactor=boundsfactor,boundstype=boundstype,themin = False)
    if (len(underlay)!=1):
        umini = get_lims(underlay,boundsfactor=boundsfactor,boundstype=boundstype,themin = True)
        umaxi = get_lims(underlay,boundsfactor=boundsfactor,boundstype=boundstype,themin = False)


    plt.hold(True)
    if (len(underlay) == 1):
        img = plt.pcolormesh(yr,xr,field,cmap=mycmap)
    else:
        img = plt.pcolormesh(yr,xr,underlay,cmap=mycmap)
    plt.plot(r_bcz*sintheta,r_bcz*costheta,'k--',[0,1],[0,0],'k--')
    plt.axis('equal')
    plt.axis('off')
    if (cbar):
        cbar = fig.colorbar(img,orientation='horizontal', shrink=0.5, aspect = 15)
        cbar.set_label(units)
        
        tick_locator = ticker.MaxNLocator(nbins=5)
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.ax.tick_params(labelsize=fontsize)   #font size for the ticks

        t = cbar.ax.xaxis.label
        t.set_fontsize(fontsize)  # font size for the axis title
    if (len(underlay) == 1):
        plt.clim((mini,maxi))
    else:
        plt.clim((umini,umaxi))
    plt.xlim((0,1))
    plt.ylim((-1,1))
    plt.plot(r[0]*sintheta,r[0]*costheta,'k')
    plt.plot(r[n_r-1]*sintheta,r[n_r-1]*costheta,'k')
    plt.plot([0,0],[-r[n_r-1],r[n_r-1]],'k--')
    plt.plot([0,0],[-r[0],-r[n_r-1]],'k',[0,0],[r[n_r-1],r[0]],'k')
    plt.hold(False)

    plt.hold(True)
    levs=mini+np.linspace(1,nlevs,nlevs)/float(nlevs)*(maxi-mini)
    plt.contour(yr,xr,field,colors='w',levels=levs)
    plt.hold(False)
def streamfunction(vr,vt,r,cost,order=0):
    """------------------------------------------------------------
    This routine takes as input a divergenceless axisymmetric 
    vector field in spherical coordinates and computes from 
    it a streamfunction (a.k.a. a flux flunction).  The grid
    is decribed by r and costheta and can be non-uniform.
   ------------------------------------------------------------
    INPUTS:
   
    Vr, Vtheta = the 2-d vector velocity (or magnetic) field.
                 Dimensions are (N_Theta,N_R)
    r,cost     = the radius and cos(colatitude) of the grid.
                 r is assumed to vary from rmax to rmin and 
                 costheta from  1 to -1 (i.e. 90 degrees
                 to -90 degrees in latitude).
                 Dimensions are r(N_R), costheta(N_Theta)
    order      = If greater than zero, integration begins at the
                 outer shell and the north pole and proceeds
                 inward and southward.  If less than zero,
                 integration begins at the inner shell and 
                 south pole and proceeds upward and northward.
                 If equal to zero, both are done and an average
                 is taken.
   ------------------------------------------------------------
    OUTPUTS:
   
    psi = the streamfunction
   ------------------------------------------------------------
    """

    [n_t,n_r]=shape(vr)
    dtheta = np.zeros(n_t)
    dr     = np.zeros(n_r)

    psi = np.zeros((n_t,n_r))

    dpsi_dr = np.zeros((n_t,n_r))
    dpsi_dt = np.zeros((n_t,n_r))

    theta = np.arccos(cost)
    sint  = sqrt(1.0-cost**2)

    for i in r_[0:n_t]:
        dpsi_dr[i,:] = -r*sint[i]*vt[i,:]
        dpsi_dt[i,:] = r*r*sint[i]*vr[i,:]

    if (order >= 0):
        # double precision accumulation
        dtheta[1:n_t] = theta[1:n_t]-theta[0:n_t-1]
        dr[1:n_r] = r[1:n_r]-r[0:n_r-1]

        dtheta[0]=0 
        dr[0]=0

        for i in r_[1:n_r]:
            psi[1:n_t,i] = psi[1:n_t,i-1] + dpsi_dr[1:n_t,i]*dr[i]
        for i in r_[1:n_t]:
            psi[i,1:n_r] = psi[i-1,1:n_r] + dpsi_dt[i,1:n_r]*dtheta[i]

    if (order <= 0):
        psi2=np.zeros((n_t,n_r))
        
        dtheta[0:n_t-1] = theta[0:n_t-1]-theta[1:n_t]
        dr[0:n_r-1] = r[0:n_r-1]-r[1:n_r]
        
        dtheta[n_t-1]=0 
        dr[n_r-1]=0
        
        for i in r_[0:n_r-1][::-1]:
            psi[0:n_t-1,i] = psi[0:n_t-1,i+1] + dpsi_dr[0:n_t-1,i]*dr[i]
        for i in r_[0:n_t-1][::-1]:
            psi[i,0:n_r-1] = psi[i+1,0:n_r-1] + dpsi_dt[i,0:n_r-1]*dtheta[i]
        
        if (order < 0):
            return psi2
        else:
            psi=0.5*(psi+psi2)
            
    return psi
