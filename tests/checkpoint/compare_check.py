from rayleigh_diagnostics import G_Avgs, SPH_Modes
import numpy as np

def chisq(a1,a2):
    b = (a1-a2)**2
    c = np.sum(b)/np.sum(a1*a1)
    return c**0.5

dirs = ['old_format', 'new_format', 'new_format_same_res', 'new_format_up_res', 'new_format_down_res']
iters = ['00000100', '00000150', '00000200','00000250']
tols = [1e-10, 1e-10, 1e-2, 1e-2]
ndirs = len(dirs)
for i in range(1,ndirs):
    #print('=========================')
    chisqs = []
    old = dirs[i-1]
    new = dirs[i]
    istring = iters[i-1]
    #print(old,new,istring)
    
    sphm_o = SPH_Modes(istring,path=old+'/SPH_Modes/')
    sphm_n = SPH_Modes(istring,path=new+'/SPH_Modes/')
    ga_o = G_Avgs(istring,path=old+'/G_Avgs/')
    ga_n = G_Avgs(istring,path=new+'/G_Avgs/')
    
    mx_ref = np.max([sphm_o.vals.real**2,sphm_o.vals.imag**2])**0.5
    for j, l in enumerate(sphm_o.lvals):
        for m in range(0,l+1):
            val_o = sphm_o.vals[m,j,0,0,:]
            val_n = sphm_n.vals[m,j,0,0,:]
            c = chisq(val_o.real,val_n.real)

            mx = np.max([val_o.real**2,val_n.imag**2])**0.5
            chisqs.append(c*mx/mx_ref)
            
            if (m > 0):  # m = 0 has no imaginary component
                c = chisq(val_o.imag,val_n.imag)
                chisqs.append(c*mx/mx_ref)
            if ((l == 4000) and (m ==4)):
                fig,ax = plt.subplots(ncols=2,figsize=(10,5))
                ax[0].plot(val_o.real)
                ax[0].plot(val_n.real)
                ax[1].plot(val_o.imag)
                ax[1].plot(val_n.imag)
                plt.show()
    sph_csq = np.max(chisqs)

    val_o = ga_o.vals[:,0]
    val_n = ga_n.vals[:,0]
    ga_csq = chisq(val_o,val_n)

    #print('sph modes: ', sph_csq)
    #print('G_Avgs: ', ga_csq) 
    if ( (ga_csq > tols[i-1]) or (sph_csq > tols[i-1]) ):
        print('Checkpoint Test Error:  Time series do not agree.')
        exit(1)
print('Checkpoint Test Passed')
exit(0)
