"""
Plot grid of Meridional Slices to verify helicity components

R. Orvedahl 9-19-2018
"""
from __future__ import print_function
from rayleigh_diagnostics import Meridional_Slices, build_file_list, plot_azav
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os

def main():

    helicity_type = 1

    saveplot = True
    tindex = -1
    phi_index = 0

    if (helicity_type == 1): # full
        vr_index = 1;   vt_index = 2;   vp_index = 3
        wr_index = 301; wt_index = 302; wp_index = 303
        hr_index = 324; ht_index = 325; hp_index = 326
        helicity = 339
        print("full helicity check")
    elif (helicity_type == 2): # u-prime, w-mean
        vr_index = 4;   vt_index = 5;   vp_index = 6
        wr_index = 307; wt_index = 308; wp_index = 309
        hr_index = 336; ht_index = 337; hp_index = 338
        helicity = 343
        print("prime-mean helicity check")
    elif (helicity_type == 3): # u-mean, w-prime
        vr_index = 7;   vt_index = 8;   vp_index = 9
        wr_index = 304; wt_index = 305; wp_index = 306
        hr_index = 333; ht_index = 334; hp_index = 335
        helicity = 342
        print("mean-prime helicity check")
    elif (helicity_type == 4): # u-mean, w-mean
        vr_index = 7;   vt_index = 8;   vp_index = 9
        wr_index = 307; wt_index = 308; wp_index = 309
        hr_index = 330; ht_index = 331; hp_index = 332
        helicity = 341
        print("mean-mean helicity check")
    elif (helicity_type == 5): # u-prime, w-prime
        vr_index = 4;   vt_index = 5;   vp_index = 6
        wr_index = 304; wt_index = 305; wp_index = 306
        hr_index = 327; ht_index = 328; hp_index = 329
        helicity = 340
        print("prime-prime helicity check")

    # build file list and get data
    files = build_file_list(0, 10000000, path="./Meridional_Slices/")
    F = Meridional_Slices(filename=files[-1], path='')

    radius = F.radius
    costh  = F.costheta
    sinth  = F.sintheta

    # extract data
    vr = F.vals[phi_index, :, :, F.lut[vr_index], tindex]
    vt = F.vals[phi_index, :, :, F.lut[vt_index], tindex]
    vp = F.vals[phi_index, :, :, F.lut[vp_index], tindex]
    wr = F.vals[phi_index, :, :, F.lut[wr_index], tindex]
    wt = F.vals[phi_index, :, :, F.lut[wt_index], tindex]
    wp = F.vals[phi_index, :, :, F.lut[wp_index], tindex]
    hr = F.vals[phi_index, :, :, F.lut[hr_index], tindex]
    ht = F.vals[phi_index, :, :, F.lut[ht_index], tindex]
    hp = F.vals[phi_index, :, :, F.lut[hp_index], tindex]
    hel = F.vals[phi_index,:, :, F.lut[helicity], tindex]

    # True value
    hr_T = vr*wr
    ht_T = vt*wt
    hp_T = vp*wp
    hel_T = hr_T + ht_T + hp_T

    # setup grid
    fig = plt.figure(1, dpi=100, figsize=(9,9))
    Grid = gridspec.GridSpec(ncols=3, nrows=4)

    ax = fig.add_subplot(Grid[0,0])
    plot_azav(fig, ax, hr_T, radius, costh, sinth, mycmap='RdYlBu_r',
              boundsfactor=4.5, boundstype='rms')
    ax.set_title("Radial True")

    ax = fig.add_subplot(Grid[1,0])
    plot_azav(fig, ax, ht_T, radius, costh, sinth, mycmap='RdYlBu_r',
              boundsfactor=4.5, boundstype='rms')
    ax.set_title("Theta True")

    ax = fig.add_subplot(Grid[2,0])
    plot_azav(fig, ax, hp_T, radius, costh, sinth, mycmap='RdYlBu_r',
              boundsfactor=4.5, boundstype='rms')
    ax.set_title("Phi True")

    ax = fig.add_subplot(Grid[3,0])
    plot_azav(fig, ax, hel_T, radius, costh, sinth, mycmap='RdYlBu_r',
              boundsfactor=4.5, boundstype='rms')
    ax.set_title("Total True")

    ax = fig.add_subplot(Grid[0,1])
    plot_azav(fig, ax, hr, radius, costh, sinth, mycmap='RdYlBu_r',
              boundsfactor=4.5, boundstype='rms')
    ax.set_title("Radial Diag")

    ax = fig.add_subplot(Grid[1,1])
    plot_azav(fig, ax, ht, radius, costh, sinth, mycmap='RdYlBu_r',
              boundsfactor=4.5, boundstype='rms')
    ax.set_title("Theta Diag")

    ax = fig.add_subplot(Grid[2,1])
    plot_azav(fig, ax, hp, radius, costh, sinth, mycmap='RdYlBu_r',
              boundsfactor=4.5, boundstype='rms')
    ax.set_title("Phi Diag")

    ax = fig.add_subplot(Grid[3,1])
    plot_azav(fig, ax, hel, radius, costh, sinth, mycmap='RdYlBu_r',
              boundsfactor=4.5, boundstype='rms')
    ax.set_title("Total Diag")

    ax = fig.add_subplot(Grid[0,2])
    plot_azav(fig, ax, np.abs(hr_T-hr), radius, costh, sinth, mycmap='RdYlBu_r',
              boundsfactor=4.5, boundstype='rms')
    ax.set_title("Radial Error")
    print("Max error, radial", np.amax(np.abs(hr_T-hr)))

    ax = fig.add_subplot(Grid[1,2])
    plot_azav(fig, ax, np.abs(ht_T-ht), radius, costh, sinth, mycmap='RdYlBu_r',
              boundsfactor=4.5, boundstype='rms')
    ax.set_title("Theta Error")
    print("Max error, theta", np.amax(np.abs(ht_T-ht)))

    ax = fig.add_subplot(Grid[2,2])
    plot_azav(fig, ax, np.abs(hp_T-hp), radius, costh, sinth, mycmap='RdYlBu_r',
              boundsfactor=4.5, boundstype='rms')
    ax.set_title("Phi Error")
    print("Max error, phi", np.amax(np.abs(hp_T-hp)))

    ax = fig.add_subplot(Grid[3,2])
    plot_azav(fig, ax, np.abs(hel_T-hel), radius, costh, sinth, mycmap='RdYlBu_r',
              boundsfactor=4.5, boundstype='rms')
    ax.set_title("Total Error")
    print("Max error, total", np.amax(np.abs(hel_T-hel)))

    fig.tight_layout()

    if (saveplot):
        output = "test_helicity_{}.png".format(helicity_type)
        plt.savefig(output, bbox_inches='tight', dpi=300)
    else:
        plt.show()

if __name__ == "__main__":

    main()

