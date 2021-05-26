"""
Module to interface with Rayleigh MathLayer, including:

    + Fourier Transforms/Derivatives
    + Legendre Transforms/Derivatives
    + Chebyshev Transforms/Derivatives

R. Orvedahl May 24, 2021
"""
from __future__ import print_function
import numpy as np
import scipy

def swap_axis(array, axis0, axis1):
    """
    Swap axes of array, other axes remain untouched

    Args
    ----
    array : ndarray
        The input array
    axis0 : int
        The original axis that will be swapped
    axis1 : int
        The final destination axis

    Returns
    -------
    array : ndarray
        The input array with the swapped axes

    Examples
    --------
    >>> B.shape
    (2,3,4,5,6)
    >>> A = swap_axis(B, -1, 1)   # swap the 2nd axis with the last axis
    >>> A.shape
    (2,6,4,5,3)
    >>> A = swap_axis(A, 1, -1)   # undo the previous swap
    >>> A.shape
    (2,3,4,5,6)
    """
    array = np.asarray(array)
    dim = len(np.shape(array))
    if (dim != 0):
        axis0 = pos_axis(axis0, dim)
        axis1 = pos_axis(axis1, dim)
        if (axis0 == axis1):
            return array
        return np.swapaxes(array, axis0, axis1)
    else:
        return array

def pos_axis(axis, dim):
    """
    Convert axis to a positive integer into array of size dim

    Args
    ----
    axis : int
        The axis to verify
    dim : int
        Number of dimensions where axis is valid

    Returns
    -------
    axis : int
        The positive integer reference to the given axis

    Example
    -------
    >>> x.shape
    (128,64,8)
    >>> dim = len(x.shape)
    >>> axis = -2
    >>>
    >>> Axis = pos_axis(axis, dim)
    >>> Axis
    1
    """
    if (axis >= 0):
        if (axis >= dim):
            raise ValueError("axis must be < dim, axis={}, dim={}".format(axis, dim))
        pos_axis = axis
    else:
        if (abs(axis) > dim):
            raise ValueError("|axis| must be <= dim, axis={}, dim={}".format(axis, dim))
        pos_axis = np.arange(dim)[axis]
    return pos_axis

def more_dimensions(C, Ndim, prepend=False, axis=None):
    """
    Increase dimensionality of array, newly added axes will be of size 1.

    Args
    ----
    C : ndarray
        Array to modify, must have dimension <= Ndim
    Ndim : int
        New dimensionality of array
    axis : int, optional
        Convert "old" axis so it points to the same data set in the new array.
    prepend : bool, optional
        Add any necessary dimensions at beginning of shape

    Returns
    -------
    C : ndarray
        Modified array with Ndim dimensions. Newly formed axes have size 1.
    added_axes : tuple of ints
        Collection of axes into new array that were added
    axis : int
        Modified axis that points to same data as the given axis, only returned
        if axis was provided.

    Example
    -------
    >>> x.shape
    (128,64,8)
    >>> X,new = more_dimensions(x, 5)
    >>> X.shape
    (128,64,8,1,1)
    >>> new
    (3,4)
    >>>
    >>> x.shape
    (128,64,8)
    >>> axis = 1  # first axis of x has 64 elements
    >>>
    >>> X, new, Axis = more_dimensions(x, 5, axis=axis, prepend=True)
    >>> X.shape
    (1,1,128,64,8)
    >>> Axis      # Axis references same data as before with 64 elements
    3
    >>> new
    (0,1)
    """
    C = np.asarray(C)
    shp = np.shape(C); dim = len(shp)

    if (axis is not None):
        paxis = pos_axis(axis, dim)

    if (Ndim == dim): # don't need to do anything
        added = ()
        if (axis is not None):
            return C, added, paxis
        else:
            return C, added

    if (Ndim < dim):
        raise ValueError("Cannot decrease dimensionality, dim={}, Ndim={}".format(dim,Ndim))

    Ndiff = Ndim - dim

    new = (1,)*Ndiff
    if (prepend):
        shape = new + shp
        added = tuple(np.arange(Ndiff))
        new_axis = paxis + Ndiff
    else:
        shape = shp + new
        added = tuple(np.arange(len(shape))[-Ndiff:])
        new_axis = paxis

    C = np.reshape(C, shape)

    if (axis is not None):
        return C, added, new_axis
    else:
        return C, added

def grid_size(N, spectral, dealias=1.0):
    """
    Find size of physical and spectral grids

    Args
    ----
    N : int
        Number of grid points or coefficients
    spectral : bool
        Does the incoming `N` refer to physical space (Ngrid) or spectral (Npoly_max)
    dealias : float, optional
        Amount to dealias: N_grid = dealias*(N_poly_max + 1)

    Returns
    -------
    Ngrid : int
        Number of physical space grid points
    Npoly_max : int
        Maximum degree in polynomial expansion
    """
    # N_grid = dealias*(N_poly_max + 1)
    if (spectral):
        Npoly_max = N
        Ngrid = int(dealias*(Npoly_max+1))
    else:
        Ngrid = N
        Npoly_max = int(Ngrid/dealias - 1)
    return Ngrid, Npoly_max

class Fourier:

    def __init__(self, N):
        """
        initialize the Fourier grid and FFT (assumes real data)

        Args
        ----
        N : int
            Number of physical space grid points
        """
        self.N = N
        self.x = self.grid()
        self.dx = np.mean(self.x[1:]-self.x[:-1])

        # alias
        self.phi = self.x

        # compute frequencies
        self.freq = self.frequencies()
        self.low_freq = self.freq.min()
        self.high_freq = self.freq.max()

        self.angular_freq = 2*np.pi*self.freq

    def frequencies(self):
        """
        calculate frequencies on the Fourier grid

        Returns
        -------
        freq : 1D array
            Array of positive frequencies
        """
        freq = np.fft.fftfreq(self.N, d=self.dx) # all frequencies: positive & negative
        freq = np.abs(freq[0:int(self.N/2)+1]) # only keep positive frequencies
        return freq

    def grid(self):
        """
        calculate Fourier grid points in [0,2pi)

        Returns
        -------
        x : 1D array
            Physical space grid points
        """
        self.period = 2.*np.pi
        dphi = self.period/self.N
        xgrid = dphi*np.arange(self.N)
        return xgrid

    def ToSpectral(self, data_in, axis=0, window=None):
        """
        FFT from physical space to spectral space

        Args
        ----
        data_in : ndarray
            Input data array of real values
        axis : int, optional
            The axis along which the FFT will be taken
        window : 1D array, optional
            Apply a window in physical space before the FFT

        Returns
        -------
        data_out : ndarray
            Complex Fourier coefficients
        """
        data_in = np.asarray(data_in)
        shp = np.shape(data_in); dim = len(shp)
        axis = pos_axis(axis, dim)

        if (shp[axis] != self.N):
            e = "FFT expected length={} along axis={}, found N={}"
            raise ValueError(e.format(self.N, axis, shp[axis]))

        if (window is None): # build window of 1.0 if not supplied
            window = np.ones((self.N))
        elif (np.shape(window)[0] != self.N):
            e = "FFT window length={}, expected={}"
            raise ValueError(e.format(np.shape(window)[0], self.N))

        if (len(np.shape(window)) > 1): raise ValueError("FFT window must be 1D")

        # reshape window for compatibility with incoming data
        nshp = [1]*dim; nshp[axis] = -1; nshp = tuple(nshp)
        window = np.reshape(window, nshp)

        norm = 2./self.N # two because we neglect negative freqs
        data_out = norm*np.fft.rfft(data_in*window, axis=axis)

        return data_out

    def ToPhysical(self, data_in, axis=0):
        """
        FFT from spectral space to physical space

        Args
        ----
        data_in : ndarray
            Spectral space coefficients
        axis : int, optional
            The axis along which the FFT will be taken

        Returns
        -------
        data_out : ndarray
            Array of physical space values
        """
        data_in = np.asarray(data_in)
        shp = np.shape(data_in); dim = len(shp)
        axis = pos_axis(axis, dim)

        # determine size of physical space for normalization
        Nm = shp[axis]
        if (Nm % 2 == 0):
            npts = 2*Nm - 1
        else:
            npts = 2*Nm - 2

        if (shp[axis] != self.freq.shape[0]):
            e = "FFT expected length={} along axis={}, found N={}"
            raise ValueError(e.format(self.freq.shape[0], axis, shp[axis]))

        norm = 2./npts
        data_out = np.fft.irfft(data_in/norm, axis=axis)

        # throw out the imaginary part, since data is assumed real
        data_out = data_out.real

        return data_out

    def d_dphi(self, data_in, axis=0, physical=True):
        """
        Take a phi derivative

        Args
        ----
        data_in : ndarray
            Input data array in spectral space
        axis : int,optional
            Axis along which the derivative will be computed
        physical : bool,optional
            Specify the incoming data as being in physical space or spectral

        Returns
        -------
        data_out : ndarray
            Output data containing d/dphi
        """
        data_in = np.asarray(data_in)
        shp = np.shape(data_in); dim = len(shp)
        axis = pos_axis(axis, dim)

        if (physical):
            if (self.N != shp[axis]):
                e = "d_dphi shape error: expected N={}, found={}"
                raise ValueError(e.format(np.shape(self.freq)[0], shp[axis]))
        else:
            if (np.shape(self.freq)[0] != shp[axis]):
                e = "d_dphi shape error: expected N={}, found={}"
                raise ValueError(e.format(np.shape(self.freq)[0], shp[axis]))

        if (physical): # transform to spectral space first, if necessary
            data_in = self.ToSpectral(data_in, axis=axis)

        nshp = [1]*dim; nshp[axis] = -1; nshp = tuple(nshp)
        freq = np.reshape(self.freq, nshp)

        # multiply by 2*pi*i to compute derivative
        data_out = 2*np.pi*1j*freq*data_in

        if (physical): # transform back to incoming space
            data_out = self.ToPhysical(data_out, axis=axis)

        return data_out

def _evaluate_Pl(x, n):
    """
    Evaluate n-th Legendre polynomial at the given grid point:

        P_{n+1} = (2*n+1)*x/(n+1)*P_n - n/(n+1)*P_n-1

        P_{m} = ( (2*m-1)*x*P_{m-1} - (m-1)*P_{m-2} ) /m

    Args
    ----
    x : float or (N,) ndarray
        The evaluation points
    n : int
        The order of the Legendre function

    Returns
    -------
    pn : float or (N,) ndarray
        Legendre function evaluated at x
    deriv_pn : float or (N,) ndarray
        Derivative of the Legendre function evaluated at x
    """
    x = np.asarray(x)
    if (len(np.shape(x)) > 0): # x is array
        length = len(x)
    else:
        length = 1             # x is scalar

    pn_minus1 = 0.0
    pn = 1.0

    # use recursion relation
    for j in range(n):
        pn_minus2 = pn_minus1
        pn_minus1 = pn
        pn = ((2.*j + 1.)*x*pn_minus1 - j*pn_minus2)/float(j+1.)

    # get derivative
    deriv_pn = n*(x*pn-pn_minus1)/(x*x - 1.)

    return pn, deriv_pn

def compute_Pl(x, lmax):
    """
    Compute array of modified associated Legendre functions for m=0. The computed
    values include the spherical harmonic normalization:

        Y_l^m = A_l^m * P_l^m(x) * exp(i*m*phi)

    This routine computes A_l^m*P_l^m(x) for all available values of x and l, but only m=0.

    Args
    ----
    x : (Nth,) ndarray
        The Legendre grid points
    lmax : int
        The maximum order of the Legendre polynomials that will be included

    Returns
    -------
    Pl : (Nth, lmax+1) ndarray
        The l-th Legendre polynomial evaluated at x[i]
    """
    x = np.asarray(x)
    n = np.shape(x)[0]

    Pl = np.zeros((n,lmax+1))

    mv = 0 # azimuthal wavenumber

    # start with the l=m & l=m+1 pieces

    # compute factorial ratio = sqrt[ (2m)! / 4**m / m! / m!] for m=0
    ratio = 1.0

    amp = np.sqrt((mv+0.5)/(2.*np.pi))
    amp *= ratio
    for i in range(n):
        xi = x[i]
        tmp = 1. - xi*xi
        if (mv%2 == 1):
            Pl[i,mv] = -amp*tmp**(mv/2+0.5) # odd m
        else:
            Pl[i,mv] = amp*tmp**(mv/2) # even m

        # l=m+1 part
        if (mv < lmax):
            Pl[i,mv+1] = Pl[i,mv]*xi*np.sqrt(2.*mv+3)

    # l>m+1 part
    for l in range(mv+2,lmax+1):
        amp = np.sqrt( ((l-1)**2 - mv*mv) / (4.*(l-1)**2 - 1.) )
        amp2 = np.sqrt( (4.*l*l-1.) / (l*l-mv*mv) )
        for i in range(n):
            xi = x[i]
            tmp = Pl[i,l-1]*xi - amp*Pl[i,l-2]
            Pl[i,l] = tmp*amp2

    return Pl

class Legendre:

    def __init__(self, N, spectral=False, dealias=1.5, parity=True):
        """
        initialize the Legendre grid and transform

        Args
        ----
        N : int
            Resolution of the theta grid
        spectral : bool, optional
            Does N refer to physical space or spectral space. If spectral=True,
            N would be the maximum polynomial degree (N_poly_max). The default
            is that N refers to the physical space resolution (N_grid)
        dealias : float, optional
            Amount to dealias: N_grid = dealias*(N_poly_max + 1)
        parity : bool, optional
            Exploit parity properties when computing the Legendre transform
        """
        self.nth, self.lmax = grid_size(N, spectral, dealias)
        self.nell    = self.lmax + 1
        self.dealias = dealias
        self.parity  = parity

        # generate grid and weights
        self.x, self.w = self.grid(self.nth)

        # alias
        self.sinth = np.sqrt(1.- self.x*self.x)
        self.costh = self.x
        self.theta = np.arccos(self.costh)

        # build array of P_l(x)
        self.Pl = compute_Pl(self.x, self.lmax)

    def grid(self, Npts):
        """
        Calculate the Legendre grid points ordered as x[i] < x[i+1] with x in (-1,1)

        Args
        ----
        Npts : int
            Number of grid points

        Returns
        -------
        x : (Npts,) ndarray
            The Legendre grid points
        w : (Npts,) ndarray
            The Legendre integration weights
        """
        x = np.zeros((Npts)); w = np.zeros((Npts))
        midpoint = 0.0
        scaling =  1.0
        n_roots = int((Npts + 1)/2)

        eps = 3.e-15

        # this is a Newton-Rhapson find for the roots
        for i in range(n_roots):
            ith_root = np.cos(np.pi*(i+0.75)/(Npts+0.5))
            converged = False; iters = 0
            while (not converged):
                pn, deriv_pn = _evaluate_Pl(ith_root, Npts)

                new_guess = ith_root - pn/deriv_pn
                delta = np.abs(ith_root - new_guess)
                ith_root = new_guess
                if (delta <= eps):
                    converged = True

                x[i] = midpoint - scaling*ith_root
                x[Npts-1-i] = midpoint + scaling*ith_root

                w[i] = 2.*scaling/((1.-ith_root*ith_root)*deriv_pn*deriv_pn)
                w[Npts-1-i] = w[i]
                iters += 1

        return x, w

    def ToSpectral(self, data_in, axis=0):
        """
        Transform from physical space to spectral space

        Args
        ----
        data_in : ndarray
            Input data to be transformed, must be dimension 4 or less
        axis : int,optional
            The axis over which the transform will take place

        Returns
        -------
        data_out : ndarray
            Transformed data, same shape as input, except along the transformed axis
        """
        data_in = np.asarray(data_in)
        shp = np.shape(data_in); dim = len(shp)
        axis = pos_axis(axis, dim)

        if (shp[axis] != self.nth):
            e = "LegendreTransform expected length={} along axis={}, found N={}"
            raise ValueError(e.format(self.nth, axis, shp[axis]))

        # increase dimensionality of input data, as necessary...magic 4 refers to number of for loops below
        data_in, added_axes, axis = more_dimensions(data_in, 4, prepend=False, axis=axis)

        # make the transform axis first
        transform_axis = 0
        data_in = swap_axis(data_in, axis, transform_axis)

        # build proper weights for transform
        iPl = 2*np.pi*np.reshape(self.w, (self.nth,1))*self.Pl # (nth,lmax+1)

        # perform the transform with calls to scipy BLAS
        alpha = 1.0; beta = 0.0

        shape = list(data_in.shape); shape[transform_axis] = self.lmax+1; shape = tuple(shape)
        data_out = np.zeros((shape))

        n, ny, nz, nt = np.shape(data_in)
        for j in range(nt):
            for i in range(nz):
                data_out[:,:,i,j] = scipy.linalg.blas.dgemm(alpha=alpha, beta=beta,
                                                            trans_a=1, trans_b=0,
                                                            a=iPl, b=data_in[:,:,i,j])
        # restore axis order
        data_out = swap_axis(data_out, transform_axis, axis)

        # remove any dimensions that were added
        data_out = np.squeeze(data_out, axis=added_axes)

        return data_out

    def ToPhysical(self, data_in, axis=0):
        """
        Transform from spectral space to physical space

        Args
        ----
        data_in : ndarray
            Input data to be transformed, must be dimension 4 or less
        axis : int,optional
            The axis over which the transform will take place

        Returns
        -------
        data_out : ndarray
            Transformed data, same shape as input, except along the transformed axis
        """
        data_in = np.asarray(data_in)
        shp = np.shape(data_in); dim = len(shp)
        axis = pos_axis(axis, dim)

        if (shp[axis] != self.lmax+1):
            e = "LegendreTransform expected length={} along axis={}, found N={}"
            raise ValueError(e.format(self.nth, axis, shp[axis]))

        # increase dimensionality of input data, as necessary...magic 4 refers to number of for loops below
        data_in, added_axes, axis = more_dimensions(data_in, 4, prepend=False, axis=axis)

        # make the transform axis first
        transform_axis = 0
        data_in = swap_axis(data_in, axis, transform_axis)

        # perform the transform with calls to scipy BLAS
        alpha = 1.0; beta = 0.0

        shape = list(data_in.shape); shape[transform_axis] = self.nth; shape = tuple(shape)
        data_out = np.zeros((shape))

        n, ny, nz, nt = np.shape(data_in)
        for j in range(nt):
            for i in range(nz):
                data_out[:,:,i,j] = scipy.linalg.blas.dgemm(alpha=alpha, beta=beta,
                                                            trans_a=0, trans_b=0,
                                                            a=self.Pl, b=data_in[:,:,i,j])
        # restore axis order
        data_out = swap_axis(data_out, transform_axis, axis)

        # remove any dimensions that were added
        data_out = np.squeeze(data_out, axis=added_axes)

        return data_out

    def d_dtheta(self, data_in, axis=0, physical=True):
        """
        Compute derivative with respect to theta

        Args
        ----
        data_in : ndarray
            Input data array, must be dimension 4 or less
        axis : int,optional
            The axis over which the derivative will take place
        physical : bool,optional
            Specify the incoming data as being in physical space or spectral

        Returns
        -------
        data_out : ndarray
            Output data containing d/dtheta
        """
        data_in = np.asarray(data_in)
        shp = np.shape(data_in); dim = len(shp)
        axis = pos_axis(axis, dim)

        if (physical):
            if (shp[axis] != (self.nth)):
                e = "d_dtheta expected length={} along axis={}, found N={}"
                raise ValueError(e.format(self.nth, axis, shp[axis]))
        else:
            if (shp[axis] != (self.lmax+1)):
                e = "d_dtheta expected length={} along axis={}, found N={}"
                raise ValueError(e.format(self.lmax+1, axis, shp[axis]))

        # increase dimensionality of input data, as necessary
        data_in, added_axes, axis = more_dimensions(data_in, 4, prepend=False, axis=axis)

        # make the d/dx axis first
        transform_axis = 0
        data_in = swap_axis(data_in, axis, transform_axis)

        if (physical): # transform to spectral space first, if necessary
            data_in = self.ToSpectral(data_in, axis=transform_axis)

        # compute the derivative: sin(th)*dF/dth with m=0
        #     sin(th)*dP_l^m/dth = l*D_{l+1}^m*P_{l+1}^m - (l+1)*D_l^m*P_{l-1}^m
        #     where P_l^m includes the A_l^m normalization of the spherical harmonics
        #     and D_l^m = sqrt[ ((l-m)*(l+m)) / ((2l-1)*(2l+1)) ]
        Dl0 = np.zeros((self.lmax+1))
        for l in range(self.lmax+1):
            num = l*l
            den = 4.*l*l-1.
            Dl0[l] = np.sqrt(num/den)

        shape = list(data_in.shape); shape[transform_axis] = self.lmax+1; shape = tuple(shape)
        data_out = np.zeros((shape))

        data_out[0,...] = 0.0 # the l=0 component

        for l in range(self.lmax): # interior l values
            data_out[l,...] = l*Dl0[l+1]*data_in[l+1,...] - (l+1)*Dl0[l]*data_in[l-1,...]

        # l=lmax component
        data_out[self.lmax,...] = -(self.lmax+1)*Dl0[self.lmax]*data_in[self.lmax-1,...]

        # convert back to physical space
        if (physical):
            data_out = self.ToPhysical(data_out, axis=transform_axis)

        # undo the sin(th) part
        nshp = [1]*len(data_out.shape); nshp[axis] = -1; nshp = tuple(nshp)
        sinth = np.reshape(self.sinth, nshp)
        data_out /= sinth

        # restore axis order
        data_out = swap_axis(data_out, transform_axis, axis)

        # remove any dimensions that were added
        data_out = np.squeeze(data_out, axis=added_axes)

        return data_out

