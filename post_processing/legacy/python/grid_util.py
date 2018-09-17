import numpy as np
def gen_cheby_grid(nx,xmin,xmax):
    """Returns a Chebyshev grid (using zeros, not extrema)
       spanning xmin to xmax, with nx points.
    """
    xgrid = np.zeros(nx,dtype='float64')
    theta = np.zeros(nx,dtype='float64')
    ctheta = np.zeros(nx,dtype='float64')
    arg = np.zeros(1,dtype='float64')
    tmp = np.zeros(1,dtype='float64')
    tmp[0] = nx
    dctheta = np.pi/tmp[0]
    tmp[0] = 0.5
    ctheta0 = dctheta*tmp[0]
    arg = ctheta0
    for i in range(nx):
        ctheta[i] = arg
        xgrid[i] = np.cos(arg)
        arg = arg+dctheta 
    dx = xmax-xmin
    xgrid = (xgrid+1.0)*0.5*dx+xmin

############
#  OLD FORTRAN
#		Integer :: i
#		Real*8 :: arg
#        ! Calculate the colocation points X { -1 , 1
#        ! Also calculate the theta grid ctheta { 0 , pi 
#	    Allocate(x(1:N_max))
#        Allocate(ctheta(1:N_max))
#        If (use_extrema) Then  !Use the extrema of T_{N_max}  
#            dctheta = pi/(N_max-1)
#            ctheta0 = 0.0d0
#            arg = ctheta0
#            Do i = 1, N_max
#                ctheta(i) = arg
#                x(i) = cos(arg)
#                arg = arg+dctheta
#            Enddo


#        Else !Use the zeroes of T_{N_max}
#
#            dctheta = pi/N_max
#            ctheta0 = dctheta*0.5d0
#		    arg = ctheta0
#		    Do i = 1, N_Max
#                ctheta(i) = arg
#			    x(i) = cos(arg)
#			    arg = arg+dctheta
#		    Enddo
#
#       Endif

##########3    
    return xgrid
