# X.  Modal Outputs
# ========
# 
# 
# **Summary:**    Spherical Harmonic Spectral Coefficients sampled at discrete radii and degree ell. 
# 
# **Subdirectory:**  SPH_Modes
# 
# **main_input prefix:** sph_mode
# 
# **Python Classes:** SPH_Modes
# 
# 
# 
# **Additional Namelist Variables:**  
# 
# * sph_mode_levels (indicial) : indices along radial grid at which to output spectral coefficients.
# 
# * sph_mode_levels_nrm (normalized) : normalized radial grid coordinates at which to output spectral coefficients.
# 
# * sph_mode_ell : Comma-separated list of spherical harmonic degree ell to output.
# 
# 
# The Modal output type allows us to output a restricted set of complex spherical harmonic coefficients at discrete radii.  For each specified ell-value, all associated azimuthal wavenumbers are output.
# 
# This output can be useful for storing high-time-cadence spectral data for a few select modes.  In the example below, we illustrate how to read in this output type, and we plot the temporal variation of the real and complex components of radial velocity for mode ell = 4, m = 4.
# 
# 
# Examining the *main_input* file, we see that the following output values have been denoted for the Shell Spectra (see *rayleigh_output_variables.pdf* for mathematical formulae):
# 
# 
# | Menu Code  | Description |
# |------------|-------------|
# | 1          | Radial Velocity |
# | 2          | Theta Velocity |
# | 3          | Phi Velocity  |
# 
# We also see that ell=2,4,8 have been selected in the *main_input* file, leading to power at the following modes:
# 
# |ell-value | m-values |
# |----------|---------------|
# | 2        | 0,1,2 |
# | 4        | 0,1,2,3,4 |
# | 8        | 0,1,2,3,4,5,6,7,8 |
# 
# 

# In[26]:

from rayleigh_diagnostics import SPH_Modes, build_file_list
import matplotlib.pyplot as plt
import numpy

qind = 1  # Radial velocity
rind = 0  # First radius stored in file


files = build_file_list(0,1000000,path='SPH_Modes')
nfiles = len(files)
for i in range(nfiles):
    spm = SPH_Modes(files[i],path='')
    if (i == 0):
        nell = spm.nell
        nr = spm.nr
        nq = spm.nq
        niter = spm.niter
        lvals = spm.lvals
        max_ell = numpy.max(lvals)
        nt = niter*nfiles
        vr = spm.lut[qind]
        vals=numpy.zeros( (max_ell+1,nell,nr,nq,nt),dtype='complex64')
        time=numpy.zeros(nt,dtype='float64')
    vals[:,:,:,:, i*niter:(i+1)*niter] = spm.vals
    time[i*niter:(i+1)*niter]=spm.time

#####################################################3
# Print some information regarding the bookkeeping
print('...........')
print(' Contents')
print('  nr = ', nr)
print('  nq = ', nq)
print('  nt = ', nt)
for i in range(nell):
    lstring=str(lvals[i])
    estring = 'Ell='+lstring+' Complex Amplitude : vals[0:'+lstring+','+str(i)+',0:nr-1,0:nq-1,0:nt-1]'
    print(estring)
print(' First dimension is m-value.')
print('...........')

######################################
# Create a plot of the ell=4, m=4 real and imaginary amplitudes
radius = spm.radius[rind]
lfour_mfour = vals[4,1,rind,vr,:]
fig, ax = plt.subplots()
ax.plot(time,numpy.real(lfour_mfour), label='real part')
ax.plot(time,numpy.imag(lfour_mfour), label='complex part')
ax.set_xlabel('Time')
ax.set_ylabel('Amplitude')
rstring = "{0:4.2f}".format(radius)
ax.set_title(r'Radial Velocity ( $\ell=4$ , m=4, radius='+rstring+' ) ')
ax.legend(shadow=True)
ax.set_xlim([0.5,4.0])
savefile = 'SPH_Modes.pdf'
print('Saving figure to: ', savefile)
plt.savefig(savefile)
