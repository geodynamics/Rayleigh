"""
Interface to the Rayleigh grid and spectral transforms, including:

    + Fourier grid/transforms/derivatives
    + Legendre grid/transforms/derivatives
    + Chebyshev grid/transforms/derivatives

R. Orvedahl May 24, 2021
"""
from __future__ import print_function
import numpy as np
from scipy.linalg.blas import dgemm as DGEMM
from scipy.linalg.blas import zgemm as ZGEMM

def check_dims(array, max_dims):
    """
    Verify array has less dimensions than allowed maximum

    Args
    ----
    array : ndarray
        The input array
    max_dims : int
        The maximum allowed number of dimensions
    """
    array = np.asarray(array)
    shp = array.shape
    ndims = len(shp)
    if (ndims > max_dims):
        d = "Array must have <= {} dimensions, given array has dims = {}"
        raise ValueError(e.format(max_dims, ndims))

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
        paxis = axis
    else:
        if (abs(axis) > dim):
            raise ValueError("|axis| must be <= dim, axis={}, dim={}".format(axis, dim))
        paxis = np.arange(dim)[axis]
    return paxis

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
        e = "Cannot decrease dimensionality, Ndim must be >= {}, Ndim = {}"
        raise ValueError(e.format(dim, Ndim))

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
    """
    Handle real data on a Fourier grid, i.e., longitude/phi grid

    Attributes
    ----------
    nphi : integer
        Physical space grid resolution. This can also be accessed as n_phi.
    phi : ndarray (nphi,)
        The longitude grid points
    dphi : float
        The grid spacing
    mvals : ndarray (nphi/2 + 1,)
        The angular frequencies, also accessed as angular_freq.

    Methods
    -------
    to_spectral(data, axis=0, window=None)
        Transform to spectral space along the given axis. A window function can be
        applied in physical space before doing the transform by setting the window
        to a ndarray of shape (nphi,).
    to_physical(data, axis=0)
        Transform to physical space along the given axis
    d_dphi(data, axis=0, physical=True)
        Compute a derivative with respect to phi along the given axis. Data is assumed
        to be in physical space (physical=True)
    """

    def __init__(self, N):
        """
        Initialize the Fourier grid and FFT (assumes real data)

        Args
        ----
        N : int
            Number of physical space grid points
        """
        self.N = N
        self.x = self.grid()
        self.dx = np.mean(self.x[1:]-self.x[:-1])

        # alias
        self.dphi = self.dx
        self.phi = self.x
        self.nphi = self.n_phi = self.N

        # compute frequencies
        self.freq = self.frequencies()
        self.low_freq = self.freq.min()
        self.high_freq = self.freq.max()

        self.angular_freq = 2*np.pi*self.freq
        self.mvals = self.angular_freq

    def frequencies(self):
        """
        Calculate frequencies on the Fourier grid

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
        Calculate Fourier grid points in [0,2pi)

        Returns
        -------
        x : 1D array
            Physical space grid points
        """
        self.period = 2.*np.pi
        dphi = self.period/self.N
        xgrid = dphi*np.arange(self.N)
        return xgrid

    def to_spectral(self, data_in, axis=0, window=None):
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

    def to_physical(self, data_in, axis=0):
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
            Input data array
        axis : int, optional
            Axis along which the derivative will be computed
        physical : bool, optional
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
            data_in = self.to_spectral(data_in, axis=axis)

        nshp = [1]*dim; nshp[axis] = -1; nshp = tuple(nshp)
        freq = np.reshape(self.freq, nshp)

        # multiply by 2*pi*i to compute derivative
        data_out = 2*np.pi*1j*freq*data_in

        if (physical): # transform back to incoming space
            data_out = self.to_physical(data_out, axis=axis)

        return data_out

def legendre_grid(Npts, quad=False):
    """
    Calculate the Legendre grid points ordered as x[i] < x[i+1] with x in (-1,1)

    Args
    ----
    Npts : int
        Number of grid points
    quad : bool, optional
        Return arrays using quad precision, default is double

    Returns
    -------
    x : (Npts,) ndarray
        The Legendre grid points
    w : (Npts,) ndarray
        The Legendre integration weights
    """
    x = np.zeros((Npts), dtype=np.float64); w = np.zeros((Npts), dtype=np.float64)
    _x = np.zeros((Npts), dtype=np.float128); _w = np.zeros((Npts), dtype=np.float128)
    midpoint = np.asarray(0.0, dtype=np.float128)
    scaling = np.asarray(1.0, dtype=np.float128)
    n_roots = int((Npts + 1)/2)

    eps = np.asarray(3.e-15, dtype=np.float128)

    # this is a Newton-Rhapson find for the roots
    for i in range(n_roots):
        ith_root = np.asarray(np.cos(np.pi*(i+0.75)/(Npts+0.5)), dtype=np.float128)
        converged = False; iters = 0
        while (not converged):
            pn, deriv_pn = _evaluate_Pl(ith_root, Npts)

            new_guess = ith_root - pn/deriv_pn
            delta = np.abs(ith_root - new_guess)
            ith_root = new_guess
            if (delta <= eps):
                converged = True

            _x[i] = midpoint - scaling*ith_root
            _x[Npts-1-i] = midpoint + scaling*ith_root

            _w[i] = 2.*scaling/((1.-ith_root*ith_root)*deriv_pn*deriv_pn)
            _w[Npts-1-i] = _w[i]
            iters += 1

    # store in double precision
    x[:] = _x[:]
    w[:] = _w[:]

    if (quad):
        return _x, _w
    else:
        return x, w

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
    x = np.asarray(x, dtype=np.float128)
    if (len(np.shape(x)) > 0): # x is array
        length = len(x)
    else:
        length = 1             # x is scalar

    pn_minus1 = np.asarray(0.0, dtype=np.float128)
    pn = np.asarray(1.0, dtype=np.float128)

    # use recursion relation
    for j in range(n):
        pn_minus2 = pn_minus1
        pn_minus1 = pn
        pn = ((2.*j + 1.)*x*pn_minus1 - j*pn_minus2)/(j+1.)

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
    x = np.asarray(x, dtype=np.float128)
    n = np.shape(x)[0]

    Pl = np.zeros((n,lmax+1), dtype=np.float64)

    # compute in "quad" precision...it will actually be closer to 80 bit
    _Pl = np.zeros((n,lmax+1), dtype=np.float128)

    mv = 0 # azimuthal wavenumber

    # start with the l=m & l=m+1 pieces

    # compute factorial ratio = sqrt[ (2m)! / 4**m / m! / m!] for m=0
    ratio = np.asarray(1.0, dtype=np.float128)

    amp = np.sqrt((mv+0.5)/(2.*np.pi))
    amp *= ratio
    for i in range(n):
        xi = x[i]
        tmp = 1. - xi*xi
        if (mv%2 == 1):
            _Pl[i,mv] = -amp*tmp**(mv/2+0.5) # odd m
        else:
            _Pl[i,mv] = amp*tmp**(mv/2) # even m

        # l=m+1 part
        if (mv < lmax):
            _Pl[i,mv+1] = _Pl[i,mv]*xi*np.sqrt(2.*mv+3)

    # l>m+1 part
    for l in range(mv+2,lmax+1):
        amp = np.sqrt( ((l-1)**2 - mv*mv) / (4.*(l-1)**2 - 1.) )
        amp2 = np.sqrt( (4.*l*l-1.) / (l*l-mv*mv) )
        for i in range(n):
            xi = x[i]
            tmp = _Pl[i,l-1]*xi - amp*_Pl[i,l-2]
            _Pl[i,l] = tmp*amp2

    # store in double precision
    Pl[:,:] = _Pl[:,:]

    return Pl

class Legendre:
    """
    Handle data on a Legendre grid, i.e., latitude/theta grid

    Attributes
    ----------
    nth : integer
        Physical space grid resolution. This can also be accessed as ntheta or n_theta.
    lmax : integer
        Maximum polynomial degree in spectral space. This can also be accessed as l_max.
    nl : integer
        Number of polynomials in spectral space. This can also be accessed as nell.
    theta : ndarray (nth,)
        The co-latitude grid points
    costh : ndarray (nth,)
        The cosine of the co-latitude points. This can also be accessed as costheta.
    sinth : ndarray (nth,)
        The sine of the co-latitude points. This can also be accessed as sintheta.

    Methods
    -------
    to_spectral(data, axis=0)
        Transform to spectral space along the given axis
    to_physical(data, axis=0)
        Transform to physical space along the given axis
    d_dtheta(data, axis=0, physical=True)
        Compute a derivative with respect to theta along the given axis. Data is assumed
        to be in physical space (physical=True)
    """

    def __init__(self, N, spectral=False, dealias=1.5):
        """
        Initialize the Legendre grid and transform

        Args
        ----
        N : int
            Resolution of the theta grid
        spectral : bool, optional
            Does N refer to physical space or spectral space. If spectral=True,
            N would be the maximum polynomial degree (l_max). The default
            is that N refers to the physical space resolution (N_theta)
        dealias : float, optional
            Amount to dealias: N_theta = dealias*(l_max + 1)
        """
        self.nth, self.lmax = grid_size(N, spectral, dealias)
        self.nell    = self.lmax + 1
        self.dealias = dealias
        self.parity  = False

        # generate grid, weights, and Pl array
        _x, _w = legendre_grid(self.nth, quad=True)
        self.x = np.zeros((self.nth), dtype=np.float64)
        self.w = np.zeros((self.nth), dtype=np.float64)
        self.x[:] = _x[:]
        self.w[:] = _w[:]

        # build array of P_l(x)
        self.Pl = compute_Pl(_x, self.lmax)

        # alias
        self.ntheta = self.nth
        self.n_theta = self.ntheta
        self.l_max = self.lmax
        self.nl = self.nell
        self.sinth = np.sqrt(1.- self.x*self.x)
        self.costh = self.x
        self.theta = np.arccos(self.costh)
        self.costheta = self.costh
        self.sintheta = self.sinth

    def to_spectral(self, data_in, axis=0):
        """
        Legendre transform from physical space to spectral space

        Args
        ----
        data_in : ndarray
            Input data to be transformed, must be dimension 4 or less
        axis : int, optional
            The axis over which the transform will take place

        Returns
        -------
        data_out : ndarray
            Transformed data, same shape as input, except along the transformed axis
        """
        data_in = np.asarray(data_in)
        shp = np.shape(data_in); dim = len(shp)
        axis = pos_axis(axis, dim)

        check_dims(data_in, 4) # only coded for 4D or less

        if (shp[axis] != self.nth):
            e = "Legendre transform expected length={} along axis={}, found N={}"
            raise ValueError(e.format(self.nth, axis, shp[axis]))

        # increase dimensionality of input data, as necessary...assumes max of 4 dimensions
        data_in, added_axes, axis = more_dimensions(data_in, 4, prepend=False, axis=axis)

        # make the transform axis first
        transform_axis = 0
        data_in = swap_axis(data_in, axis, transform_axis)

        # build proper weights for transform
        iPl = 2*np.pi*np.reshape(self.w, (self.nth,1))*self.Pl # (nth,lmax+1)

        # perform the transform with calls to BLAS
        alpha = 1.0; beta = 0.0

        shape = list(data_in.shape); shape[transform_axis] = self.lmax+1; shape = tuple(shape)
        data_out = np.zeros((shape))

        # data assumed to be 4D, so two for loops plus a matrix multiply
        n, ny, nz, nt = np.shape(data_in)
        for j in range(nt):
            for i in range(nz):
                data_out[:,:,i,j] = DGEMM(alpha=alpha, beta=beta,
                                          trans_a=1, trans_b=0,
                                          a=iPl, b=data_in[:,:,i,j])
        # restore axis order
        data_out = swap_axis(data_out, transform_axis, axis)

        # remove any dimensions that were added
        data_out = np.squeeze(data_out, axis=added_axes)

        return data_out

    def to_physical(self, data_in, axis=0):
        """
        Legendre transform from spectral space to physical space

        Args
        ----
        data_in : ndarray
            Input data to be transformed, must be dimension 4 or less
        axis : int, optional
            The axis over which the transform will take place

        Returns
        -------
        data_out : ndarray
            Transformed data, same shape as input, except along the transformed axis
        """
        data_in = np.asarray(data_in)
        shp = np.shape(data_in); dim = len(shp)
        axis = pos_axis(axis, dim)

        check_dims(data_in, 4) # only coded for 4D or less

        if (shp[axis] != self.lmax+1):
            e = "Legendre transform expected length={} along axis={}, found N={}"
            raise ValueError(e.format(self.nth, axis, shp[axis]))

        # increase dimensionality of input data, as necessary...assumes max of 4 dimensions
        data_in, added_axes, axis = more_dimensions(data_in, 4, prepend=False, axis=axis)

        # make the transform axis first
        transform_axis = 0
        data_in = swap_axis(data_in, axis, transform_axis)

        # perform the transform with calls to BLAS
        alpha = 1.0; beta = 0.0

        shape = list(data_in.shape); shape[transform_axis] = self.nth; shape = tuple(shape)
        data_out = np.zeros((shape))

        # data assumed to be 4D, so two for loops plus a matrix multiply
        n, ny, nz, nt = np.shape(data_in)
        for j in range(nt):
            for i in range(nz):
                data_out[:,:,i,j] = DGEMM(alpha=alpha, beta=beta,
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
        axis : int, optional
            The axis over which the derivative will take place
        physical : bool, optional
            Specify the incoming data as being in physical space or spectral space

        Returns
        -------
        data_out : ndarray
            Output data containing d/dtheta
        """
        data_in = np.asarray(data_in)
        shp = np.shape(data_in); dim = len(shp)
        axis = pos_axis(axis, dim)

        check_dims(data_in, 4) # only coded for 4D or less

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

        if (physical): # go to spectral space if necessary
            data_in = self.to_spectral(data_in, axis=transform_axis)

        # compute the derivative: sin(th)*dF/dth with m=0
        #     sin(th)*dP_l^m/dth = l*D_{l+1}^m*P_{l+1}^m - (l+1)*D_l^m*P_{l-1}^m
        #     where P_l^m includes the A_l^m normalization of the spherical harmonics
        #     and D_l^m = sqrt[ ((l-m)*(l+m)) / ((2l-1)*(2l+1)) ]
        Dl0 = np.zeros((self.lmax+1))
        for l in range(1,self.lmax+1):
            num = l*l
            den = 4.*l*l-1.
            Dl0[l] = np.sqrt(num/den)

        shape = list(data_in.shape); shape[transform_axis] = self.lmax+1; shape = tuple(shape)
        data_out = np.zeros((shape))

        # F = sum C_l*P_l
        # sin(th)dF/dth = sum C_l*sin(th)dP_l/dth
        #               = sum_{l=0}^lmax C_l * [ a_l*P_{l+1} + b_l*P_{l-1} ]
        #               = sum_{l=0}^lmax C_l*a_l*P_{l+1} + sum_{l=0}^lmax C_l*b_l*P_{l-1}
        #        redo sum indices: k=l+1 for first one and k=l-1 for second one
        #               = sum_{k=1}^lmax C_{k-1}*a_{k-1}*P_k
        #                  + sum_{k=0}^{lmax-1} C_{k+1}*b_{k+1}*P_k
        #               = sum_{l=1}^{lmax-1} [C_{l-1}*a_{l-1} + C_{l+1}*b_{l+1}]*P_l
        #                  + C_{lmax-1}*a_{lmax-1}*P_lmax
        # the l=0 component is zero
        for l in range(1,self.lmax): # interior l values
            data_out[l,...] = -(l+2)*Dl0[l+1]*data_in[l+1,...] + (l-1)*Dl0[l]*data_in[l-1,...]

        # l=lmax component
        data_out[self.lmax,...] = (self.lmax-1)*Dl0[self.lmax]*data_in[self.lmax-1,...]

        # convert back to physical space
        data_out = self.to_physical(data_out, axis=transform_axis)

        # undo the sin(th) part in physical space
        nshp = [1]*len(data_out.shape); nshp[transform_axis] = -1; nshp = tuple(nshp)
        sinth = np.reshape(self.sinth, nshp)
        data_out /= sinth

        if (not physical): # go back to spectral space if necessary
            data_out = self.to_spectral(data_out, axis=transform_axis)

        # restore axis order
        data_out = swap_axis(data_out, transform_axis, axis)

        # remove any dimensions that were added
        data_out = np.squeeze(data_out, axis=added_axes)

        return data_out

def chebyshev_zeros(Npts, reverse=False):
    """
    generate the Chebysheve zeros grid

    Args
    ----
    Npts : int
        The number of grid points
    reverse : bool, optional
        Reverse the grid order to be x[i] > x[i+1], i.e., x = (-1,1)

    Returns
    -------
    grid : 1D array
        The array of grid points
    """
    grid = np.zeros((Npts))
    dcth = np.pi/Npts
    arg = dcth*0.5
    for i in range(Npts):
        grid[i] = np.cos(arg)
        arg += dcth
    if (not reverse):
        grid = grid[::-1] # this will produce x[i] < x[i+1]
    return grid

class Chebyshev:
    """
    Handle data on a Chebyshev grid, i.e., radius grid

    Attributes
    ----------
    nr : integer
        The total number of radial grid points. This can also be accessed as n_r.
    radius : ndarray (nr,)
        The radial grid points
    n_domains : integer
        Number of separate radial domains
    nr_domains : ndarray (n_domains,)
        The number of radial grid points in each domain
    npoly : ndarray (n_domains,)
        The number of polynomials used in each domain
    npoly_dealias : ndarray (n_domains,)
        The number of dealiased polynomials used in each domain
    boundaries : ndarray (n_domains+1,)
        The domain boundaries
    rmin : float
        The minimum radius, i.e., the lower boundary
    rmax : float
        The maximum radius, i.e., the upper boundary
    aspect_ratio : float
        The ratio of rmin/rmax
    shell_depth : float
        The depth of the shell, rmax-rmin
    dealias : ndarray (n_domains,)
        The amount of dealiasing for each domain.

    Methods
    -------
    to_spectral(data, axis=0)
        Transform to spectral space along the given axis
    to_physical(data, axis=0)
        Transform to physical space along the given axis
    d_dr(data, axis=0, physical=True)
        Compute a derivative with respect to radius along the given axis. Data is assumed
        to be in physical space (physical=True)
    """

    def __init__(self, nr_domains,
                 rmin=None, rmax=None, aspect_ratio=None, shell_depth=None,
                 boundaries=None,
                 n_uniform_domains=1, uniform_bounds=False,
                 dealias=None,
                 dmax=3):
        """
        Initialize the Chebyshev grid and transform

        Support for three types of grids:
            a) single Chebyshev domain [default]
            b) uniform set of N Chebyshev domains
            c) N Chebyshev domains with different resolutions

        Args
        ----
        nr_domains : int or 1D array_like of ints
            Resolution of the radial grid(s). If using uniform domains, these
            entries refer to the resolution per domain.
        rmin : float, optional
            The lower boundary of the global domain
        rmax : float, optional
            The upper boundary of the global domain
        aspect_ratio : float, optional
            Set the aspect ratio of the domain, defined as rmin/rmax
        shell_depth : float, optional
            The total domain shell thickness, defined as rmax-rmin
        boundaries : 1D array_like, optional
            Set the boundaries for multiple subdomains. This overrides
            the rmin/rmax and aspect_ratio/shell_depth arguments, since the first
            element of boundaries is assumed to be rmin and the last element is
            assumed to be rmax. The number of elements in boundaries must be one
            more than the number of domain resolutions found in nr_domains.
        n_uniform_domains : int, optional
            Choose the number of Chebyshev domains, each having the same
            resolution and uniformly distributed boundaries. Only the first
            element of nr_domains will be used and is treated as the resolution
            of each subdomain. Similarly, only the first element of the dealias
            option is used and applied to each subdomain.
        uniform_bounds : bool, optional
            Uniformly distribute the boundaries across the domain, allowing different
            resolutions and dealias values for each subdomain.
        dealias : int or 1D array_like of ints
            Specify the number of polynomials to use in each subdomain by choosing
            the amount of dealiasing:
                n_polynomials = n_r - dealias
            where n_r is the grid resolution of the subdomain and n_polynomials is
            the number of Chebyshev polynomials that will be used in that domain.
            The default behavior is the standard 2/3 rule:
                n_polynomials = 2*n_r/3
        dmax : int, optional
            Allocate array storage space for up to and including the dmax-th derivative
        """
        if (not hasattr(nr_domains, "__len__")): # user gave a single value
            nr_domains = [nr_domains]

        # error check that global domain bounds were chosen
        if ((rmin is None) and (rmax is None) \
            and (aspect_ratio is None) and (shell_depth is None)):
            if (boundaries is None):
                msg = "Must specifiy global boundaries using one of the following:"
                msg += "\n\t1) set rmin & rmax"
                msg += "\n\t2) set aspect_ratio & shell_depth"
                msg += "\n\t3) set the boundaries array"
                raise ValueError(msg)
            else:
                # boundaries was specified, make sure it is consistent with nr_domains
                if (not hasattr(boundaries, "__len__")):
                    e = "The boundaries value must be 1D array-like: array/list/tuple"
                    raise ValueError(e)
                else:
                    if (len(nr_domains) != (len(boundaries)-1)):
                        e = "Number of domains must be one less than number of boundaries"
                        e += "\n\t# domains    = len(nr_domains) = {}".format(len(nr_domains))
                        e += "\n\t# boundaries = len(boundaries) = {}".format(len(boundaries))
                        raise ValueError(e)

        if (boundaries is not None):
            boundaries = np.sort(boundaries) # just to be sure
            rmin = boundaries[0]
            rmax = boundaries[-1]
            aspect_ratio = rmin/rmax
            shell_depth = rmax - rmin
        elif ((aspect_ratio is not None) and (shell_depth is not None)):
            rmax = shell_depth/(1.-aspect_ratio)
            rmin = rmax*aspect_ratio
        elif ((rmin is not None) and (rmax is not None)):
            aspect_ratio = rmin/rmax
            shell_depth = rmax - rmin
        else:
            msg = "Must specifiy global boundaries using one of the following:"
            msg += "\n\t1) set both rmin & rmax"
            msg += "\n\t2) set both aspect_ratio & shell_depth"
            msg += "\n\t3) set the boundaries array (rmin=first element, rmax=last element)"
            raise ValueError(msg)

        if (dealias is None): # default
            dealias = [-1]
        elif (not hasattr(dealias, "__len__")): # user gave a single value
            dealias = [dealias]

        n_domains = len(nr_domains)

        if (len(dealias) < n_domains): # extend unspecified entries using default
            Ndealias = len(dealias)
            diff = n_domains - Ndealias
            new = [-1]*diff
            dealias = list(dealias) + new

        if ((n_domains == 1) and (n_uniform_domains < 2)): # default case (a)
            boundaries = [rmin, rmax]
            n_r = nr_domains[0]

        if (n_uniform_domains > 1): # case (b)
            n_r = n_uniform_domains*nr_domains[0]
            n_domains = n_uniform_domains
            nr_domains = [nr_domains[0]]*n_domains
            dealias = [dealias[0] for i in range(n_domains)] # set all to same
            boundaries = [0]*(n_domains+1)
            boundaries[0] = rmin
            dr = 1.*(shell_depth)/n_domains
            for i in range(1,n_domains+1):
                boundaries[i] = boundaries[i-1] + dr

        else: # case (c)
            n_r = np.sum(nr_domains)

            if (uniform_bounds): # same as (b) above, but allows different resolutions
                n_uniform_domains = n_domains
                boundaries = [0]*(n_domains+1)
                boundaries[0] = rmin
                dr = (shell_depth)/n_domains
                for i in range(1,n_domains+1):
                    boundaries[i] = boundaries[i-1] + dr

        # store some variables
        self.rmin = boundaries[0]
        self.rmax = boundaries[-1]
        self.aspect_ratio = self.rmin/self.rmax
        self.shell_depth = self.rmax - self.rmin
        self.n_domains = n_domains
        self.nr_domains = nr_domains
        self.n_r = sum(self.nr_domains) # total grid resolution
        self.dealias = dealias
        self.boundaries = boundaries

        self.build_grid(dmax=dmax) # build the grid

        # alias
        self.radius = self.grid
        self.nr = self.n_r

        self.npoly_dealias = self.rda # number of dealiased polynomials per domain

    def build_grid(self, dmax=3):
        """
        Build the grid(s)

        Args
        ----
        dmax : int, optional
            Allocate array storage space for up to and including the dmax-th derivative
        """
        self.npoly = np.zeros((self.n_domains), dtype=np.int32)
        self.rda = np.zeros((self.n_domains), dtype=np.int32)

        for i in range(self.n_domains):
            n = self.nr_domains[i]
            self.npoly[i] = n
            db = int(2*n/3)
            self.rda[i] = db
            if (self.dealias[i] > 0):
                db = self.dealias[i]
                if ((db >= 1) and (db < n)):
                    self.rda[i] = n - db

        self.npoly = self.npoly[::-1]
        self.rda = self.rda[::-1]
        self.boundaries = self.boundaries[::-1]
        self.dealias = self.dealias[::-1]
        self.nr_domains = self.nr_domains[::-1]

        gmax = self.boundaries[0]
        gmin = self.boundaries[-1]

        self.max_npoly = self.npoly.max()
        self.ntotal = self.npoly.sum()

        self.x = np.zeros((self.max_npoly, self.n_domains), dtype=np.float64)
        self.theta = np.zeros((self.max_npoly, self.n_domains), dtype=np.float64)
        self.deriv_scaling = np.zeros((self.ntotal), dtype=np.float64)

        # use list of ndarray to mimic the rmcontainers
        self.cheby_even = [0]*self.n_domains    # cheby_even[i] = 2D array
        self.cheby_odd = [0]*self.n_domains     # cheby_odd[i] = 2D array
        self.dcheby = [0]*self.n_domains        # dcheby[i] = 3D array
        self.n_even = np.zeros((self.n_domains), dtype=np.int32)
        self.n_odd = np.zeros((self.n_domains), dtype=np.int32)

        self.find_colocation_pts()
        self.find_Tn()
        self.find_Tn_deriv_array(dmax)

        # global grid and rescale derivative arrays
        self.deriv_scaling = np.zeros((self.n_r), dtype=np.float64)
        self.integration_weights = np.zeros((self.n_r), dtype=np.float64)
        self.grid = np.zeros((self.n_r), dtype=np.float64)
        ind = 0
        for n in range(self.n_domains):
            n_max = self.npoly[n]
            ind2 = ind + n_max - 1

            upperb = self.boundaries[n]
            lowerb = self.boundaries[n+1]
            xmin = self.x[n_max-1,n]
            xmax = self.x[0,n]

            domain_delta = upperb - lowerb

            scaling = domain_delta/(xmax - xmin)
            self.grid[ind:ind2+1] = lowerb + scaling*(self.x[:n_max, n] - xmin)

            int_scale = 3*np.pi*scaling/( (gmax**3 - gmin**3)*n_max )

            scaling = 1./scaling

            for i in range(dmax):
                self.dcheby[n][:,:,i] *= scaling**(i+1)

            for i in range(n_max):
                gind = ind + i
                xx = self.x[i,n]
                self.integration_weights[gind] = \
                                   int_scale*self.grid[gind]**2*np.sqrt(1.-xx*xx)
                self.deriv_scaling[gind] = scaling
            self.integration_weights[ind] *= 0.5 # boundaries x 1/2
            self.integration_weights[ind2] *= 0.5

            ind = ind2 + 1

    def find_colocation_pts(self):
        """
        Compute the colocation grid points
        """
        for n in range(self.n_domains):
            n_max = self.npoly[n]

            # zeros-based grid
            dth = np.pi/n_max
            arg = dth*0.5
            for i in range(n_max):
                self.theta[i,n] = arg
                self.x[i,n] = np.cos(arg)
                arg += dth

    def find_Tn(self):
        """
        Compute the Chebyshev polynomials evaluated at the colocation grid points
        """
        cheby = np.zeros((self.max_npoly, self.max_npoly), dtype=np.float64)

        for n in range(self.n_domains):
            n_max = self.npoly[n]
            for r in range(n_max):
                for i in range(n_max):
                    arg = i*self.theta[r,n]
                    cheby[r,i] = np.cos(arg)

            n_odd = int(n_max/2)
            n_even = n_odd + (n_max%2)
            n_x = n_even
            self.n_odd[n] = n_odd
            self.n_even[n] = n_even

            self.cheby_even[n] = np.zeros((n_x,n_even), dtype=np.float64)
            self.cheby_odd[n] = np.zeros((n_x,n_odd), dtype=np.float64)

            for i in range(n_even):
                self.cheby_even[n][:,i] = cheby[:n_x,2*i]
            for i in range(n_odd):
                self.cheby_odd[n][:,i] = cheby[:n_x,2*i+1]

            if (n_x != n_odd): # have x=0 point, don't double count power when using parity
                self.cheby_even[n][n_x-1,:] *= 0.5
                self.cheby_odd[n][n_x-1,:] *= 0.5

        cheby = None

    def find_Tn_deriv_array(self, dmax):
        """
        Compute derivative of Chebyshev polynomials evaluated at the colocation grid points
        """
        alpha = np.zeros((self.max_npoly, self.max_npoly), dtype=np.float64)

        for m in range(self.n_domains):
            n_max = self.npoly[m]

            alpha[:,:] = 0.0
            alpha[n_max-1,:] = 0.0
            alpha[n_max-2,n_max-1] = 2.*(n_max-2)

            for k in range(n_max-3, -1, -1): # includes n_max-3 & 0
                alpha[k,k+1] = 2.*k
                alpha[k,:] = alpha[k,:] + alpha[k+2,:]

            self.dcheby[m] = np.zeros((n_max,n_max,dmax), dtype=np.float64)
            for r in range(n_max):
                for i in range(n_max):
                    arg = k*self.theta[r,m]
                    self.dcheby[m][r,i,0] = np.cos(arg)

            self.dcheby[m][:,0,0] = self.dcheby[m][:,0,0] - 0.5

            if (dmax >= 1):
                for d in range(1,dmax):
                    for n in range(n_max):
                        for k in range(n_max):
                            for i in range(n_max):
                                self.dcheby[m][k,n,d] += self.dcheby[m][k,i,d-1]*alpha[i,n]

        alpha = None

    def dealias_buffer(self, data_in, axis=0):
        """
        Dealias the data array

        Args
        ----
        data_in : ndarray
            Input data array, must be dimension 4 or less
        axis : int, optional
            The axis over which the dealiasing will take place

        Returns
        -------
        data_out : ndarray
            Dealiased data
        """
        data_in = np.asarray(data_in)
        shp = np.shape(data_in); dim = len(shp)
        axis = pos_axis(axis, dim)

        check_dims(data_in, 4) # only coded for 4D or less

        # increase dimensionality of input data, as necessary
        data_in, added_axes, axis = more_dimensions(data_in, 4, prepend=False, axis=axis)

        transform_axis = 0
        data_in = swap_axis(data_in, axis, transform_axis)

        _, ni, nj, nk = data_in.shape
        for k in range(nk):
            for j in range(nj):
                for i in range(ni):
                    offset = 0
                    for n in range(self.n_domains):
                        ind1 = self.rda[n] + offset
                        ind2 = self.npoly[n] + offset
                        data_in[ind1:ind2,i,j,k] = 0.0
                        offset += self.npoly[n]

        data_out = swap_axis(data_in, transform_axis, axis)

        # remove any dimensions that were added
        data_out = np.squeeze(data_in, axis=added_axes)

        return data_out

    def to_spectral(self, data_in, axis=0):
        """
        Chebyshev transform from physical space to spectral space

        Args
        ----
        data_in : ndarray
            Input data to be transformed, must be dimension 4 or less
        axis : int, optional
            The axis over which the transform will take place

        Returns
        -------
        data_out : ndarray
            Transformed data, same shape as input, except along the transformed axis
        """
        data_in = np.asarray(data_in)
        shp = np.shape(data_in); dim = len(shp)
        axis = pos_axis(axis, dim)

        check_dims(data_in, 4) # only coded for 4D or less

        if (shp[axis] != self.n_r):
            e = "Chebyshev transform expected length={} along axis={}, found N={}"
            raise ValueError(e.format(self.n_r, axis, shp[axis]))

        # increase dimensionality of input data, as necessary...assumes max of 4 dimensions
        data_in, added_axes, axis = more_dimensions(data_in, 4, prepend=False, axis=axis)

        # make the transform axis first
        transform_axis = 0
        data_in = swap_axis(data_in, axis, transform_axis)

        # perform the transform with calls to BLAS
        beta = 0.0

        shape = list(data_in.shape); shape[transform_axis] = self.ntotal; shape = tuple(shape)
        data_out = np.zeros((shape))

        ng, n2, n3, n4 = data_in.shape
        hoff = 0
        for hh in range(self.n_domains):
            n_even = self.n_even[hh]
            n_odd = self.n_odd[hh]
            n_max = self.npoly[hh]
            n_x = self.n_even[hh]
            alpha = 2./n_max

            #c_temp = np.zeros((n_even, n2, n3, n4)) # done below with call to BLAS
            f_even = np.zeros((n_x, n2, n3, n4))
            f_odd  = np.zeros((n_x, n2, n3, n4))

            # build even/odd input
            for kk in range(n4):
                for k in range(n3):
                    for j in range(n2):
                        for i in range(n_x):
                            f_even[i,j,k,kk] = data_in[hoff+i,j,k,kk] + \
                                               data_in[hoff+n_max-1-i,j,k,kk]
                            f_odd[i,j,k,kk]  = data_in[hoff+i,j,k,kk] - \
                                               data_in[hoff+n_max-1-i,j,k,kk]

            # call BLAS
            c_temp = DGEMM(alpha=alpha, beta=beta,
                           trans_a=1, trans_b=0,
                           a=self.cheby_even[hh],
                           b=f_even[:,:,:,:])
            c_temp = np.reshape(c_temp, (n_even, n2, n3, n4))

            # unpack the output
            for kk in range(n4):
                for k in range(n3):
                    for j in range(n2):
                        for i in range(n_even):
                            data_out[hoff+2*i,j,k,kk] = c_temp[i,j,k,kk]

            #if (n_even != n_odd):
            #    c_temp = np.zeros((n_odd, n2, n3, n4))

            # call BLAS
            c_temp = DGEMM(alpha=alpha, beta=beta,
                           trans_a=1, trans_b=0,
                           a=self.cheby_odd[hh],
                           b=f_odd[:,:,:,:])
            c_temp = np.reshape(c_temp, (n_odd, n2, n3, n4))

            # unpack the output
            for kk in range(n4):
                for k in range(n3):
                    for j in range(n2):
                        for i in range(n_odd):
                            data_out[hoff+2*i+1,j,k,kk] = c_temp[i,j,k,kk]

            hoff += self.npoly[hh]

        # restore axis order
        data_out = swap_axis(data_out, transform_axis, axis)

        # remove any dimensions that were added
        data_out = np.squeeze(data_out, axis=added_axes)

        return data_out

    def to_physical(self, data_in, axis=0):
        """
        Chebyshev transform from spectral space to physical space

        Args
        ----
        data_in : ndarray
            Input data to be transformed, must be dimension 4 or less
        axis : int, optional
            The axis over which the transform will take place

        Returns
        -------
        data_out : ndarray
            Transformed data, same shape as input, except along the transformed axis
        """
        data_in = np.asarray(data_in)
        shp = np.shape(data_in); dim = len(shp)
        axis = pos_axis(axis, dim)

        check_dims(data_in, 4) # only coded for 4D or less

        if (shp[axis] != self.ntotal):
            e = "Chebyshev transform expected length={} along axis={}, found N={}"
            raise ValueError(e.format(self.ntotal, axis, shp[axis]))

        # increase dimensionality of input data, as necessary...assumes max of 4 dimensions
        data_in, added_axes, axis = more_dimensions(data_in, 4, prepend=False, axis=axis)

        # make the transform axis first
        transform_axis = 0
        data_in = swap_axis(data_in, axis, transform_axis)

        # perform the transform with calls to BLAS
        alpha = 1.0; beta = 0.0

        shape = list(data_in.shape); shape[transform_axis] = self.n_r; shape = tuple(shape)
        data_out = np.zeros((shape))

        ng, n2, n3, n4 = data_in.shape
        hoff = 0
        for hh in range(self.n_domains):
            n_even = self.n_even[hh]
            n_odd = self.n_odd[hh]
            n_max = self.npoly[hh]
            n_x = self.n_even[hh]

            c_temp = np.zeros((n_even, n2, n3, n4))
            #f_temp = np.zeros((n_x, n2, n3, n4))

            for kk in range(n4):
                for k in range(n3):
                    for j in range(n2):
                        for i in range(n_even):
                            c_temp[i,j,k,kk] = data_in[hoff+2*i,j,k,kk]

            # call BLAS
            f_temp = DGEMM(alpha=alpha, beta=beta,
                           trans_a=0, trans_b=0,
                           a=self.cheby_even[hh],
                           b=c_temp[:,:,:,:])
            f_temp = np.reshape(f_temp, (n_x, n2, n3, n4))

            # unpack the output
            for kk in range(n4):
                for k in range(n3):
                    for j in range(n2):
                        for i in range(n_even):
                            data_out[hoff+i,j,k,kk] = f_temp[i,j,k,kk]
                            data_out[hoff+n_max-i-1,j,k,kk] = f_temp[i,j,k,kk]

            if (n_even != n_odd):
                data_out[hoff+n_x-1,:,:,:] *= 2.0
                c_temp = np.zeros((n_odd, n2, n3, n4))

            for kk in range(n4):
                for k in range(n3):
                    for j in range(n2):
                        for i in range(n_odd):
                            c_temp[i,j,k,kk] = data_in[hoff+2*i+1,j,k,kk]

            # call BLAS
            f_temp = DGEMM(alpha=alpha, beta=beta,
                           trans_a=0, trans_b=0,
                           a=self.cheby_odd[hh],
                           b=c_temp[:,:,:,:])
            f_temp = np.reshape(f_temp, (n_x, n2, n3, n4))

            # unpack the output
            for kk in range(n4):
                for k in range(n3):
                    for j in range(n2):
                        for i in range(n_odd):
                            data_out[hoff+i,j,k,kk] += f_temp[i,j,k,kk]
                            data_out[hoff+n_max-i-1,j,k,kk] -= f_temp[i,j,k,kk]

            if (n_even != n_odd):
                data_out[hoff+n_x-1,:,:,:] += 2.0*f_temp[hoff+n_x-1,:,:,:]

            hoff += self.npoly[hh]

        hoff = 0
        for hh in range(self.n_domains):
            n_max = self.npoly[hh]
            istart = hoff
            iend = istart + n_max
            for kk in range(n4):
                for k in range(n3):
                    for j in range(n2):
                        data_out[istart:iend,j,k,kk] -= 0.5*data_in[hoff,j,k,kk]
 
            hoff += self.npoly[hh]

        # restore axis order
        data_out = swap_axis(data_out, transform_axis, axis)

        # remove any dimensions that were added
        data_out = np.squeeze(data_out, axis=added_axes)

        return data_out

    def d_dr(self, data_in, axis=0, physical=True):
        """
        Compute derivative with respect to radius

        Args
        ----
        data_in : ndarray
            Input data array, must be dimension 4 or less
        axis : int, optional
            The axis over which the derivative will take place
        physical : bool, optional
            Specify the incoming data as being in physical space or spectral space

        Returns
        -------
        data_out : ndarray
            Output data containing d/dr
        """
        data_in = np.asarray(data_in)
        shp = np.shape(data_in); dim = len(shp)
        axis = pos_axis(axis, dim)

        check_dims(data_in, 4) # only coded for 4D or less

        if (physical):
            if (shp[axis] != (self.n_r)):
                e = "d_dr expected length={} along axis={}, found N={}"
                raise ValueError(e.format(self.n_r, axis, shp[axis]))
        else:
            if (shp[axis] != self.ntotal):
                e = "d_dr expected length={} along axis={}, found N={}"
                raise ValueError(e.format(self.ntotal, axis, shp[axis]))

        # increase dimensionality of input data, as necessary
        data_in, added_axes, axis = more_dimensions(data_in, 4, prepend=False, axis=axis)

        # make the d/dx axis first
        transform_axis = 0
        data_in = swap_axis(data_in, axis, transform_axis)

        if (physical): # go to spectral space if necessary
            data_in = self.to_spectral(data_in, axis=transform_axis)

        data_out = np.zeros_like(data_in)

        n, n2, n3, n4 = data_in.shape
        for k in range(n3):
            for j in range(n2):
                hoff = 0
                for hh in range(self.n_domains):
                    n = self.npoly[hh]

                    data_out[hoff+n-1,j,k,:] = 0.0
                    data_out[hoff+n-2,j,k,:] = 2.*(n-1)*data_in[hoff+n-1,j,k,:]

                    for i in range(n-3,-1,-1):
                        data_out[hoff+i,j,k,:] = data_out[hoff+i+2,j,k,:] + \
                                      2.*(i+1)*data_in[hoff+i+1,j,k,:]

                    hoff += self.npoly[hh]

        # rescale
        scaling = np.reshape(self.deriv_scaling, (-1,1,1,1))
        data_out[:,:,:,:] *= scaling

        if (physical): # go back to physical space if necessary
            data_out = self.to_physical(data_out, axis=transform_axis)

        # restore axis order
        data_out = swap_axis(data_out, transform_axis, axis)

        # remove any dimensions that were added
        data_out = np.squeeze(data_out, axis=added_axes)

        return data_out

class SHT:
    """
    Handle full spherical harmonic transforms of (theta,phi) to/from (l,m)

    Attributes
    ----------
    theta : ndarray (nth,)
        The co-latitude grid points
    costh : ndarray (nth,)
        The cosine of the co-latitude points. This can also be accessed as costheta.
    sinth : ndarray (nth,)
        The sine of the co-latitude points. This can also be accessed as sintheta.
    phi : ndarray (nphi,)
        The longitude grid points

    Methods
    -------
    to_spectral(data, laxis=0, maxis=1)
        Transform to spectral space along the given axes
    to_physical(data, laxis=0, maxis=1)
        Transform to physical space along the given axes
    """

    def __init__(self, N, spectral=False, dealias=1.5, parity=True):
        """
        Initialize the SHT grid and transform

        Args
        ----
        N : int
            Resolution of the theta grid
        spectral : bool, optional
            Does N refer to physical space or spectral space. If spectral=True,
            N would be the maximum polynomial degree (l_max). The default
            is that N refers to the physical space resolution (N_theta)
        dealias : float, optional
            Amount to dealias: N_theta = dealias*(l_max + 1)
        parity : bool, optional
            Exploit parity during the Legendre transforms
        """
        self.nth, self.lmax = grid_size(N, spectral, dealias)
        self.nell    = self.lmax + 1
        self.dealias = dealias
        self.parity  = parity

        self.nphi = 2*self.nth

        # build Fourier grid/transforms
        self.fourier = Fourier(self.nphi)

        # generate grid, weights (xq = x in quad precision)
        self.xq, self.wq = legendre_grid(self.nth, quad=True)
        self.x = np.zeros((self.nth), dtype=np.float64)
        self.w = np.zeros((self.nth), dtype=np.float64)
        self.x[:] = self.xq[:]
        self.w[:] = self.wq[:]

        # alias
        self.ntheta = self.nth
        self.n_theta = self.ntheta
        self.l_max = self.lmax
        self.nl = self.nell
        self.nm = self.nl
        self.sinth = np.sqrt(1.- self.x*self.x)
        self.costh = self.x
        self.theta = np.arccos(self.costh)
        self.costheta = self.costh
        self.sintheta = self.sinth

        self.dphi = self.fourier.dphi
        self.phi = self.fourier.phi
        self.cosphi = np.cos(self.phi)
        self.sinphi = np.sin(self.phi)

        # build arrays of P_l^m(x)
        self.compute_Plm()

    def compute_Plm(self):
        """
        Compute various arrays holding the A_l^m*P_l^m data evaluated on the grid
        """
        n_m = self.nm

        n_l = np.zeros((n_m), dtype=np.int32)
        p_lm = [0]*n_m
        p_lmq = [0]*n_m
        ip_lm = [0]*n_m
        ntmax = self.nth
        if (self.parity):
            p_lm_odd = [0]*n_m
            p_lm_even = [0]*n_m
            ip_lm_odd = [0]*n_m
            ip_lm_even = [0]*n_m
            n_l_odd = np.zeros((n_m), dtype=np.int32)
            n_l_even = np.zeros((n_m), dtype=np.int32)
            lvals = [{'even':None, 'odd':None}]*n_m
            lvalsi = [{'even':None, 'odd':None}]*n_m
            ntmax = int(self.nth / 2)

        piq = np.asarray(np.pi, dtype=np.float128)
        one = np.asarray(1.0, dtype=np.float128)
        two = np.asarray(2.0, dtype=np.float128)
        four = np.asarray(4.0, dtype=np.float128)
        half = np.asarray(0.5, dtype=np.float128)

        for m in range(n_m): # m loop index is also actual m value
            n_l[m] = self.lmax - m + 1
            p_lmq[m] = np.zeros((ntmax, self.lmax+1)) # really should be (nt, m:lmax)

            factorial_ratio = np.asarray(1.0, dtype=np.float128)
            for i in range(1,m):
                factorial *= ((i-half)/i)**half

            amp = ((m+half)/(2*piq))**half
            amp *= factorial_ratio

            for i in range(ntmax):
                x = self.xq[i]
                tmp = one - x*x

                # l=m pieces
                if (m%2 == 1):
                    p_lmq[m][i,m] = -amp*tmp**(int(m/2)+half)
                else:
                    p_lmq[m][i,m] = amp*tmp**(int(m/2))

                # l = m+1 pieces
                if (m < self.lmax):
                    p_lmq[m][i,m+1] = p_lmq[m][i,m]*x*(two*m+3)**half

            # general recursion for l > m+1
            for l in range(m+2,lmax+1):
                for i in range(ntmax):
                    x = self.xq[i]
                    amp = (l-1)**2 - m*m
                    amp = amp / (four*(l-1)**2-one)
                    amp = amp**half
                    tmp = p_lmq[m][i,l-1]*x - amp*p_lmq[m][i,l-2]
                    amp = (four*l*l-one)/(l*l-m*m)
                    p_lmq[m][i,l] = tmp*amp**half

            if (self.parity): # parity resort
                if (m == 0):
                    pts_norm = 1./(2*self.nth)
                    stp_norm = 1.0
                else:
                    pts_norm = 1./(self.nth)
                    stp_norm = 0.5

                n_l_even[m] = 0
                n_l_odd[m] = 0
                for l in range(m, self.lmax+1):
                    parity_test = l-m
                    if (parity_test % 2 == 1):
                        n_l_odd[m] += 1
                    else:
                        n_l_even[m] += 1

                if (n_l_even[m] > 0):
                    ip_lm_even[m] = np.zeros((self.nth/2, n_l_even[m]))
                    p_lm_even[m] = np.zeros((n_l_even[m], self.nth/2))
                    lvals[m]['even'] = np.zeros((n_l_even[m]))
                    lvalsi[m]['even'] = np.zeros((n_l_even[m]), dtype=np.int32)
                if (n_l_odd[m] > 0):
                    ip_lm_odd[m] = np.zeros((self.nth/2, n_l_odd[m]))
                    p_lm_odd[m] = np.zeros((n_l_odd[m], self.nth/2))
                    lvals[m]['odd'] = np.zeros((n_l_odd[m]))
                    lvalsi[m]['odd'] = np.zeros((n_l_odd[m]), dtype=np.int32)

                indeven = 1; indodd = 1
                for l in range(m, self.lmax+1):
                    parity_test = l-m
                    if (parity_test % 2 == 1):
                        lvals[m]['odd'][indodd] = l
                        lvalsi[m]['odd'][indodd] = l
                        for i in range(self.nth/2):
                            renorm = two*piq*self.wq[i]
                            tmp = p_lmq[m][i,l]*renorm
                            ip_lm_odd[m][i,indodd] = tmp*pts_norm
                            p_lm_odd[m][indodd,i] = p_lmq[m][i,l]*stp_norm
                        indodd += 1
                    else:
                        lvals[m]['even'][indeven] = l
                        lvalsi[m]['even'][indeven] = l
                        for i in range(self.nth/2):
                            renorm = two*piq*self.wq[i]
                            tmp = p_lmq[m][i,l]*renorm
                            ip_lm_even[m][i,indeven] = tmp*pts_norm
                            p_lm_even[m][indeven,i] = p_lmq[m][i,l]*stp_norm
                        indeven += 1

            else:
                # store in double precision
                p_lm[m]  = np.zeros((self.lmax+1, ntmax), dtype=np.float64)

                p_lm[m] = np.transpose(p_lmq[m])
                renorm = np.reshape(two*piq*self.wq, (self.nth,1))

                ip_lm[m] = np.asarray(p_lmq[m]*renorm, dtype=np.float64)

            # try to release some memory
            p_lmq[m] = None
            if (self.parity):
                p_lm[m] = None
                ip_lm[m] = None

        # try to release some memory
        p_lmq = None
        if (self.parity):
            p_lm = None
            ip_lm = None

        # make available to class
        self.ntmax = ntmax
        if (self.parity):
            self.lvals = lvals
            self.lvalsi = lvalsi
            self.p_lm_even = p_lm_even
            self.p_lm_odd = p_lm_odd
            self.ip_lm_even = ip_lm_even
            self.ip_lm_odd = ip_lm_odd
            self.n_l_even = n_l_even
            self.n_l_odd = n_l_odd
        else:
            self.p_lm = p_lm
            self.ip_lm = ip_lm

    def to_spectral(self, data_in, th_axis=0, phi_axis=1, fft=True):
        """
        SHT transform from physical space (theta,phi) to spectral space (l,m)

        Args
        ----
        data_in : ndarray
            Input data to be transformed
        th_axis : int, optional
            The axis over which the Legendre transform (theta/l) will take place
        phi_axis : int, optional
            The axis over which the Fourier transform (phi/m) will take place
        fft : bool, optional
            If fft=True, a full SHT transform from (th,phi) to (l,m) will occur,
            starting with a FFT in the phi direction. If False, the FFT will not
            occur and the Legendre transform will assume m=0.

        Returns
        -------
        data_out : ndarray
            Transformed data, same shape as input, except along the transformed axes
        """
        data_in = np.asarray(data_in)
        shp = np.shape(data_in); dim = len(shp)
        th_axis = pos_axis(th_axis, dim)
        phi_axis = pos_axis(phi_axis, dim)

        if (dim < 2 and fft):
            e = "Full SHT requires data to be at least 2d. Found {}d"
            raise ValueError(e.format(dim))

        if (shp[th_axis] != self.nth):
            e = "SHT: Legendre transform expected length={} along axis={}, found N={}"
            raise ValueError(e.format(self.nth, th_axis, shp[th_axis]))

        if (fft and (shp[phi_axis] != self.nphi)):
            e = "SHT: Fourier transform expected length={} along axis={}, found N={}"
            raise ValueError(e.format(self.nphi, phi_axis, shp[phi_axis]))

        if (fft): # first do FFT
            data_in = self.fourier.to_spectral(data_in, axis=phi_axis)

            # trim the m-axis to be consistent with triangular truncation
            slc = [slice(None)]*dim
            slc[phi_axis] = slice(0, self.nm)
            data_in = data_in[tuple(slc)]

        # make the transform axis first
        transform_axis = 0
        data_in = swap_axis(data_in, axis, transform_axis)

        # allocate space
        shape = list(data_in.shape); shape[transform_axis] = self.nl
        data_out = np.zeros(tuple(shape))

        shape[transform_axis] = self.nth/2
        f_even = np.zeros_like(tuple(shape))
        f_odd  = np.zeros_like(tuple(shape))

        # build even/odd input data
        for i in range(int(self.nth/2)):
            f_even[i,...] = data_in[i,...] + data_in[self.nth-1-i,...]
            f_odd[i,...]  = data_in[i,...] - data_in[self.nth-1-i,...]

        alpha = 1.0; beta = 0.0
        if (fft):
            _nm = self.nm
        else:
            _nm = [0] # assume only m=0 for "pure" Legendre transform

        slc = [slice(None)]*dim
        for m in range(_nm):
            if (fft):
                slc[phi_axis] = m
            if (self.n_l_even[m] > 0):
                out_shape = list(f_even.shape)
                out_shape[transform_axis] = self.n_l_even[m]
                c_temp = ZGEMM(alpha=alpha, beta=beta,
                                trans_a=1, trans_b=0,
                                a=self.ip_lm_even[m][:,m:self.lmax+1], # (nth/2, n_l_even)
                                b=f_even[tuple(slc)])                  # (nth/2, ...)
                c_temp = np.reshape(c_temp, tuple(out_shape))          # (n_l_even, ...)

                # store output
                for j in range(1,self.n_l_even[m]+1):
                    l = self.lvalsi[m]['even'][j]
                    data_out[l,...] = c_temp[j,...]

            if (self.n_l_odd[m] > 0):
                out_shape = list(f_odd.shape)
                out_shape[transform_axis] = self.n_l_odd[m]
                c_temp = ZGEMM(alpha=alpha, beta=beta,
                                trans_a=1, trans_b=0,
                                a=self.ip_lm_odd[m][:,m:self.lmax+1], # (nth/2, n_l_odd)
                                b=f_odd[tuple(slc)])                  # (nth/2, ...)
                c_temp = np.reshape(c_temp, tuple(out_shape))         # (n_l_odd, ...)

                # store output
                for j in range(1,self.n_l_odd[m]+1):
                    l = self.lvalsi[m]['odd'][j]
                    data_out[l,...] = c_temp[j,...]

        # restore axis order
        data_out = swap_axis(data_out, transform_axis, axis)

        return data_out

    def to_physical(self, data_in, th_axis=0, phi_axis=1, fft=True):
        """
        SHT transform from spectral space (l,m) to physical space (theta,phi)

        Args
        ----
        data_in : ndarray
            Input data to be transformed
        th_axis : int, optional
            The axis over which the Legendre transform (theta/l) will take place
        phi_axis : int, optional
            The axis over which the Fourier transform (phi/m) will take place
        fft : bool, optional
            If fft=True, a full SHT transform from (l,m) to (th,phi) will occur,
            starting with a FFT in the phi direction. If False, the FFT will not
            occur and the Legendre transform will assume m=0.

        Returns
        -------
        data_out : ndarray
            Transformed data, same shape as input, except along the transformed axes
        """
        data_in = np.asarray(data_in)
        shp = np.shape(data_in); dim = len(shp)
        th_axis = pos_axis(th_axis, dim)
        phi_axis = pos_axis(phi_axis, dim)

        if (dim < 2 and fft):
            e = "Full SHT requires data to be at least 2d. Found {}d"
            raise ValueError(e.format(dim))

        if (shp[th_axis] != self.nth):
            e = "SHT: Legendre transform expected length={} along axis={}, found N={}"
            raise ValueError(e.format(self.nth, th_axis, shp[th_axis]))

        if (fft and (shp[phi_axis] != self.nphi)):
            e = "SHT: Fourier transform expected length={} along axis={}, found N={}"
            raise ValueError(e.format(self.nphi, phi_axis, shp[phi_axis]))

        if (fft): # first do FFT
            data_in = self.fourier.to_physical(data_in, axis=phi_axis)

        # make the transform axis first
        transform_axis = 0
        data_in = swap_axis(data_in, axis, transform_axis)

        # allocate space
        shape = list(data_in.shape); shape[transform_axis] = self.nth
        data_out = np.zeros(tuple(shape))

        shape[transform_axis] = self.nth/2
        out_shape = list(shape)

        for i in range(int(self.nth/2)):
            f_even[i,...] = data_in[i,...] + data_in[self.nth-1-i,...]
            f_odd[i,...]  = data_in[i,...] - data_in[self.nth-1-i,...]

        alpha = 1.0; beta = 0.0
        if (fft):
            _nm = self.nm
        else:
            _nm = [0] # assume only m=0 for "pure" Legendre transform

        slc = [slice(None)]*dim
        for m in range(_nm):
            if (fft):
                slc[phi_axis] = m

            if (self.n_l_even[m] > 0):
                shape[transform_axis] = self.n_l_even[m]
                f_even = np.zeros_like(tuple(shape))

                for j in range(1,self.n_l_even[m]+1):
                    l = self.lvalsi[m]['even'][j]
                    f_even[j,...] = data_in[l,...]

                c_temp = ZGEMM(alpha=alpha, beta=beta,
                                trans_a=1, trans_b=0,
                                a=self.p_lm_even[m][:,m:self.lmax+1], # (n_l_even, nth/2)
                                b=f_even[tuple(slc)])                 # (n_l_even, ...)
                c_temp = np.reshape(c_temp, tuple(out_shape))         # (nth/2, ...)

                # store first half of output
                data_out[:int(self.nth/2),...] = c_temp[:,...]

                # reflect evenly
                for j in range(int(self.nth/2)):
                    data_out[self.nth-j-1,...] = data_out[j,...]

            if (self.n_l_odd[m] > 0):
                shape[transform_axis] = self.n_l_odd[m]
                f_odd  = np.zeros_like(tuple(shape))

                for j in range(1,self.n_l_odd[m]+1):
                    l = self.lvalsi[m]['odd'][j]
                    f_odd[j,...] = data_in[l,...]

                c_temp = ZGEMM(alpha=alpha, beta=beta,
                                trans_a=1, trans_b=0,
                                a=self.p_lm_odd[m][:,m:self.lmax+1], # (n_l_odd, nth/2)
                                b=f_odd[tuple(slc)])                 # (n_l_odd, ...)
                c_temp = np.reshape(c_temp, tuple(out_shape))        # (nth/2, ...)

                # add output to already computed even stuff
                for j in range(int(self.nth/2)):
                    data_out[j,...] += c_temp[j,...]            # include first half
                    data_out[self.nth-j-1,...] -= c_temp[j,...] # odd-reflect

        # restore axis order
        data_out = swap_axis(data_out, transform_axis, axis)

        return data_out

