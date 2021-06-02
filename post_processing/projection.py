import numpy
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import griddata

def sph_to_cart(theta,phi):
    """ Assumes r = 1 for now.
        Converts theta,phi to x,y,z
        theta and phi should have same dimension.
        
        returns x,y,z """
    x = numpy.sin(theta)*numpy.cos(phi)
    y = numpy.sin(theta)*numpy.sin(phi)
    z = numpy.cos(theta)
    return x,y,z

def rotate_cart_z(x1,z1,alpha):
    """Rotates Cartesian coordinate system about y-axis
     by angle alpha in positive theta direction.
     As y is unchanged by this rotation,
     this routine returns x2,z2 """
    z2 = z1*numpy.cos(alpha)+x1*numpy.sin(alpha)
    x2 = -z1*numpy.sin(alpha)+x1*numpy.cos(alpha)
    return x2,z2

def cart_to_sph(x,y,z,eps=1e-4):
    """ Converts Cartesian x,y,z to spherical coordinates.
        Returns theta,phi (r assumed to be 1).
        
        x,y,z should have same dimension.
        
        Optional keyword eps is to avoid divide by zero issues when
        computing 1/sin(theta).
        Phi runs from -pi to pi.  Phi is taken to be angle
        relative to x-axis [i.e., x = 0 -> phi = 0."""
    theta = numpy.arccos(z)

    sintheta = numpy.sqrt(1-z*z)
    phi_temp = numpy.arccos(x/(sintheta+eps))

    #Note that acos only returns values from 0 through pi
    # Thus acos(pi+delta) = acos(pi-delta)
    # Since we want phi to run from -pi to pi,
    # this also means acos(pi+delta) = acos(-pi-delta)
    # Use the value of y to distinguish which angle we want
    condition = (y > 0)
    phi = numpy.where(condition,phi_temp,-phi_temp) #-numpy.pi)    
    return theta, phi

def find_viewable(phi,y,z,data=[]):
    """Returns data at coordinates viewable from the selected vantage point.
       Return values are y_v,z_v,data_v, where _v is for 'viewable'"""
    condition = (numpy.abs(phi) <= numpy.pi/2)

    y_v = numpy.extract(condition,y)
    z_v = numpy.extract(condition,z)
    if (len(data)==0):
        return y_v, z_v
    else:
        data_v = numpy.extract(condition, data)
        return y_v,   z_v, data_v

def transform_coordinates(theta,phi,alpha,data=[]):
    """ Takes theta/phi coordinates and rotates them by
        angle alpha in theta.  Returns y,z coordinates of those
        not hidden behind the sphere.  If data is provided,
        it's visible elements are also returned."""

    x1,y1,z1 = sph_to_cart(theta,phi) # Convert to Cartesian

    x2,z2 = rotate_cart_z(x1,z1,alpha) # Rotate about y-axis, by angle alpha in theta

    # Use the value of phi in the new coordinate system to identify
    # those points which are viewable, and those which are on the
    # opposite side of the sphere.  Viewable points lie in the range
    # -90 < phi < 90
    theta2, phi2 = cart_to_sph(x2,y1,z2) #Find phi in the new coordinate system
    #Return those coordinates and associated data that are viewable
    return find_viewable(phi2,y1,z2,data=data)

def rescale_data(data,x,y,ifactor=2):
    """Data is striped data[y,x]"""
    ny = data.shape[0] 
    nx = data.shape[1]
    
    xmin = numpy.min(x)
    xmax = numpy.max(x)
    ymin = numpy.min(y)
    ymax = numpy.max(y)
    
    nxnew = int(ifactor*nx)
    nynew = int(ifactor*ny)
    
    xnew = numpy.linspace(xmin,xmax,nxnew,dtype='float64')
    ynew = numpy.linspace(ymin,ymax,nynew,dtype='float64')

    # Create a 2-D grid of input coordinates
    y2di = numpy.zeros((ny,nx),dtype='float64')
    x2di = numpy.zeros((ny,nx),dtype='float64')
    for i in range(ny):
        x2di[i,:] = x[:]
    for j in range(nx):
        y2di[:,j] = y[:]      

    # Convert input grid to 1-D, reshaped array
    x1d = x2di.reshape(nx*ny)
    y1d = y2di.reshape(nx*ny)
    d1d = data.reshape(nx*ny)
        
    # Create a 2-D grid of output coordinates
    y2d = numpy.zeros((nynew,nxnew),dtype='float64')
    x2d = numpy.zeros((nynew,nxnew),dtype='float64')
    for i in range(nynew):
        x2d[i,:] = xnew[:]
    for j in range(nxnew):
        y2d[:,j] = ynew[:]  

    data_new = griddata((x1d, y1d), d1d, (x2d, y2d), method='cubic',fill_value=0)

    return data_new, x2d, y2d

def axis_range(ax): # gets subplot coordinates on a figure in "normalized"
        # coordinates
    pos = plt.get(ax, 'position')
    bottom_left = pos.p0
    top_right = pos.p1
    xmin, xmax = bottom_left[0], top_right[0]
    ymin, ymax = bottom_left[1], top_right[1]
    
    return xmin, xmax, ymin, ymax
    


def plot_ortho(data_in,costheta,fig, ax, cax, latcen=45, loncen=0, nyz=512, colormap="RdYlBu",
               nlats = 9, nlons = 8, nplot=200, gridcolor='black', 
               filename='', hwidth=2, gridstyle='-',plot_grid=True,
               gridwidth=1, scale_type=[None, None],pole_style='truncated'):
    """
    z
    |
    |   Coordinate system used.
    |   x-direction is out of page/screen.
    |   Rotation in +theta is about y-axis, toward viewer
    -------------y
    """
    import copy
    nphi = data_in.shape[0]
    ntheta = data_in.shape[1]
    ntcheck = costheta.shape[0]
    
    # Because we may want to re-interpolate onto a finer grid,
    # make a copy
    data = numpy.zeros((nphi,ntheta),dtype='float64')

    if (ntcheck != ntheta):
        print("Error: costheta and 2nd dimension of data must be equal")
        return
    
    # Initialize the input coordinate system
    theta = numpy.arccos(costheta)
    twopi = 2*numpy.pi
    dphi = twopi/nphi
    phi = numpy.linspace(0,twopi-dphi,nphi)

    # Shift phi-vantage-point to loncen
    phi = phi-loncen*numpy.pi/180.0

    # Based on desired latcen, calculate needed rotation angle in theta (alpha).
    alpha = -latcen*numpy.pi/180.0

    pre_interpolate=False
    ifactor = 2
    if (pre_interpolate):
        data, theta2d, phi2d = rescale_data(data_in,theta,phi,ifactor=ifactor)
        print('interp out: ', data.shape,theta2d.shape,phi2d.shape, numpy.max(data), numpy.min(data))
    else:
        data[:,:] = data_in[:,:]
        # Create a 2-D grid of input coordinates
        phi2d = numpy.zeros((nphi,ntheta),dtype='float64')
        theta2d = numpy.zeros((nphi,ntheta),dtype='float64')
        for i in range(nphi):
            theta2d[i,:] = theta[:]
        for j in range(ntheta):
            phi2d[:,j] = phi[:]  


        
    # Calculate y,z positions in plane-page after rotation by alpha
    # Transform_coordinates performs the following actions
    # 1.  Convert input spherical to input Cartesian coordinates
    # 2.  Rotate input coordinate system by alpha
    # 3.  Return the y,z coordinates and associated data values for those
    #     coordinates that are visible from this vantage point (i.e.,
    #     throw out anything that is on the far side of the sphere.)
    # y,z, and data_extract are 1-D arrays
    y,z,data_extract = transform_coordinates(theta2d,phi2d,alpha,data=data)
    
    
    ############################################################
    # Next, create a 2-D grid of z,y values for interpolation

    z1d = numpy.linspace(-1,1,nyz,dtype='float64')
    y1d = numpy.linspace(-1,1,nyz,dtype='float64')

    y2d = numpy.zeros((nyz,nyz),dtype='float64')
    z2d = numpy.zeros((nyz,nyz),dtype='float64')

    for j in range(nyz):
        y2d[:,j] = y1d[:]
    for i in range(nyz):
        z2d[i,:] = z1d[:]

    ################################################################
    # Interpolate 1-D data onto theta 2-D y,z grid
    val = griddata((y, z), data_extract, (y2d, z2d), method='cubic',fill_value=0)

    if (scale_type[0]=='rms'):
        dmax = scale_type[1]*numpy.sqrt(numpy.mean(val*val))
        dmin = -dmax
        print('scale: ', dmin, dmax)
    elif (scale_type[0]=='abs'):
        dmax = scale_type[1]*numpy.max(numpy.abs(val))
        dmin = -dmax
    elif (scale_type[0]=='force'):
        dmax = scale_type[1][1]
        dmin = scale_type[1][0]
    else:
        dmin = numpy.min(val)
        dmax = numpy.max(val)
    # Set anything outside the unit circle to NaN
    # This allows us to mask those points out using set_bad,
    # making it easy to show them in white regardless of color table.
    rsq = y2d*y2d+z2d*z2d
    val = numpy.where(rsq>1,numpy.nan,val)


    # Copy the color map, so we don't modify a global variable
    my_cmap = copy.copy(matplotlib.cm.get_cmap(colormap))
    my_cmap.set_bad('white')


    #fig, ax = plt.subplots(figsize=(10,10))
    dshow = ax.imshow(numpy.flip(numpy.transpose(val),axis=0), interpolation='bicubic', cmap=my_cmap,
                      vmin = dmin, vmax = dmax)


    ##########################################
    # Plot latitude lines
    if (nlats > 0 and plot_grid):
        dlat= numpy.pi/(nlats+1)
        lat_plot = numpy.zeros(nplot,dtype='float64')
        phi_plot = numpy.linspace(-numpy.pi,numpy.pi,nplot,dtype='float64')

        latvals = []
        for i in range(nlats):
            latval = dlat*(i+1)
            latvals.append(latval)
            lat_plot[:] = latval

            # Transform our lat/lon coords just like we did with the data
            y_plot,z_plot= transform_coordinates(lat_plot,phi_plot,alpha)

            ax.plot(y_plot*nyz/2+nyz/2, nyz/2-z_plot*nyz/2, color=gridcolor, 
                    linewidth=gridwidth, linestyle=gridstyle)


    ################################
    # Plot longitude lines
    if (nlons > 0 and plot_grid):
        if (pole_style == 'truncated'):
            lat_plot = numpy.linspace(numpy.min(latvals),numpy.max(latvals),nplot,dtype='float64')
        elif (pole_style == 'polar'):
            lat_plot = numpy.linspace(0,numpy.pi,nplot,dtype='float64')
        else:
            print('Unknown polar option selected -- using "truncated"')
            lat_plot = numpy.linspace(numpy.min(latvals),numpy.max(latvals),nplot,dtype='float64')
        dlon = numpy.pi/(nlons)
        for i in range(nlons):
            lonval = dlon*(i+1)
            phi_plot[:] = lonval

            y_plot,z_plot= transform_coordinates(lat_plot,phi_plot,alpha)
            ax.plot(y_plot*nyz/2+nyz/2, nyz/2-z_plot*nyz/2, color=gridcolor,
                    linewidth=gridwidth,linestyle=gridstyle)

            phi_plot[:] = lonval-numpy.pi
            y_plot,z_plot= transform_coordinates(lat_plot,phi_plot,alpha)
            ax.plot(y_plot*nyz/2+nyz/2, nyz/2-z_plot*nyz/2, color=gridcolor,
                    linewidth=gridwidth, linestyle=gridstyle)

    ##################################3
    # Finally, plot the bounding circle
    sf=0.99
    hp = 1/nyz
    yplot = numpy.linspace(-1,1,nplot,dtype='float64')
    yplot = yplot-hp*0
    zplot1 = numpy.sqrt(1-yplot*yplot) -hp
    zplot2 = -zplot1

    px=1  # data coordinates need extra pixel shift
    #print(len(yplot), len(zplot1), len(zplot2))
    ax.plot(yplot*nyz/2+nyz/2-0*px/2, nyz/2-zplot1*nyz/2, color=gridcolor, linewidth = hwidth)
    ax.plot(yplot*nyz/2+nyz/2-0*px/2, nyz/2-zplot2*nyz/2, color=gridcolor, linewidth = hwidth)

#    plt.show()

    ax.set_xlim([-10,nyz+10])
    ax.set_ylim([nyz+10,-10])
    #ax.set_xlim([-1.1,1.1])
    #ax.set_ylim([1.1,-1.1])
    ax.axis('off')
   
    # Add a color bar
    #fig.colorbar(dshow, ax=ax, orientation='horizontal', shrink=0.5)
    ax_xmin, ax_xmax, ax_ymin, ax_ymax = axis_range(ax)
    ax_delta_x = ax_xmax - ax_xmin
    ax_delta_y = ax_ymax - ax_ymin


    
    ax_center_x = ax_xmin + 0.5*ax_delta_x

    
    cbar_aspect = 1./10.
    fig_aspect = ax_delta_x/ax_delta_y # assumes subplot aspect ratio is 1
        # otherwise, must multiply by proper subplot aspect ratio (=0.5 for
        # Mollweide))
    cbar_width = 0.5*ax_delta_x # make cbar half as long as plot is wide
    cbar_height = cbar_width*cbar_aspect/fig_aspect
    cbar_bottom = ax_ymin - 1.5*cbar_height
    cbar_left = ax_xmin + 0.5*ax_delta_x - 0.5*cbar_width
    
    caxi = fig.add_axes((cbar_left, cbar_bottom, cbar_width, cbar_height))
    # Note that if you want to pin the color bar to a specific axes, you can modify
    # the caxi that is used here -- though you will likely lose control over the bar's
    # sizing.
    cbar = plt.colorbar(dshow, cax=caxi, orientation='horizontal') # , shrink=0.15, aspect=0.1)

