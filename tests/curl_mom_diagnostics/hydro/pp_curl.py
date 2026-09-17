#!/home/pgrad1/2831664s/anaconda3/bin/python

from rayleigh_diagnostics import Point_Probes, GridInfo, build_file_list

import itertools

import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# SETTINGS
# ============================================================

start_file = 1
end_file   = 10000000

# Probe location(s) to examine, as a list of (pind, tind, rind) tuples.
# Set to None to plot every sampled probe location (all combinations of
# the phi, theta, and r grids read from the Point_Probes files).
probe_points = None

output_file = 'Force_Balance_Curls_Hydro.pdf'


# ============================================================
# READ POINT PROBE FILES
# ============================================================

files = build_file_list(
    start_file,
    end_file,
    path='Point_Probes'
)

nfiles = len(files)

print(f'Number of Point_Probes files: {nfiles}')

for i, f in enumerate(files):

    print(f'Reading file {i+1}/{nfiles}: {f}')

    pp = Point_Probes(f, path='')

    if i == 0:

        nphi   = pp.nphi
        ntheta = pp.ntheta
        nr     = pp.nr
        nq     = pp.nq
        niter  = pp.niter

        vals = np.zeros(
            (nphi, ntheta, nr, nq, niter * nfiles)
        )

        time = np.zeros(niter * nfiles)

    vals[:, :, :, :, i*niter:(i+1)*niter] = pp.vals

    time[i*niter:(i+1)*niter] = pp.time


lut = pp.lut

if probe_points is None:
    probe_points = list(
        itertools.product(range(nphi), range(ntheta), range(nr))
    )

nprobes = len(probe_points)

print(f'Number of probe locations to plot: {nprobes}')


# ============================================================
# FULL-GRID EXTENT (for normalized coordinates, matching
# Rayleigh's own r_nrm/theta_nrm/phi_nrm definition, i.e. a
# linear min-max normalization over the full simulation grid)
# ============================================================

grid = GridInfo(path='./')

r_min,     r_max     = grid.radius.min(), grid.radius.max()
theta_min, theta_max = grid.theta.min(),  grid.theta.max()
phi_min,   phi_max   = grid.phi.min(),    grid.phi.max()


# ============================================================
# FUNCTION TO EXTRACT A DIAGNOSTIC AT A GIVEN PROBE LOCATION
# ============================================================

def F(idx, pind, tind, rind):
    return vals[
        pind,
        tind,
        rind,
        lut[idx],
        :
    ]


# ============================================================
# FORCE COMPONENTS AT A GIVEN PROBE LOCATION
# ============================================================

def force_components(pind, tind, rind):

    # Inertia
    FI_r = F(1301, pind, tind, rind)
    FI_t = F(1302, pind, tind, rind)
    FI_p = F(1303, pind, tind, rind)

    # Coriolis
    FC_r = -F(1327, pind, tind, rind)
    FC_t = -F(1328, pind, tind, rind)
    FC_p = -F(1329, pind, tind, rind)

    # Buoyancy
    FA_t = -F(1319, pind, tind, rind)
    FA_p = -F(1320, pind, tind, rind)

    # Viscous
    FV_r = -F(1339, pind, tind, rind)
    FV_t = -F(1340, pind, tind, rind)
    FV_p = -F(1341, pind, tind, rind)

    # Pressure
    FP_t = -F(1358, pind, tind, rind)
    FP_p = -F(1359, pind, tind, rind)

    balance_r = FI_r + FC_r + FV_r
    balance_t = FI_t + FC_t + FA_t + FV_t + FP_t
    balance_p = FI_p + FC_p + FA_p + FV_p + FP_p

    return {
        'FI': (FI_r, FI_t, FI_p),
        'FC': (FC_r, FC_t, FC_p),
        'FA': (None, FA_t, FA_p),
        'FV': (FV_r, FV_t, FV_p),
        'FP': (None, FP_t, FP_p),
        'BAL': (balance_r, balance_t, balance_p),
    }


# ============================================================
# COLORS
# ============================================================

COLORS = {
    'FI': '#1f77b4',   # blue
    'FA': '#ff7f0e',   # orange
    'FC': '#2ca02c',   # green
    'FM': '#9467bd',   # purple
    'FV': '#d62728',   # red
    'FP': '#8c564b',   # brown
    'BAL': '#000000'   # black
}

LABELS = {
    'FI': r'$\nabla\times\mathbf{F}_I$',
    'FC': r'$\nabla\times\mathbf{F}_C$',
    'FA': r'$\nabla\times\mathbf{F}_A$',
    'FV': r'$\nabla\times\mathbf{F}_V$',
    'FP': r'$\nabla\times\mathbf{F}_P$',
}

BAL_LABELS = (r'$R_r$', r'$R_\theta$', r'$R_\phi$')

COMPONENT_TITLES = ('Radial component', r'$\theta$ component', r'$\phi$ component')


# ============================================================
# FIGURE
# ============================================================

fig, ax = plt.subplots(
    nprobes,
    3,
    figsize=(27, 6 * nprobes),
    sharex=True,
    squeeze=False
)


for row, (pind, tind, rind) in enumerate(probe_points):

    forces = force_components(pind, tind, rind)

    r_val     = pp.radius[rind]
    theta_val = np.arccos(pp.costheta[tind])
    phi_val   = pp.phi[pind]

    theta_deg = np.degrees(theta_val)
    phi_deg   = np.degrees(phi_val)

    r_nrm     = (r_val - r_min) / (r_max - r_min)
    theta_nrm = (theta_val - theta_min) / (theta_max - theta_min)
    phi_nrm   = (phi_val - phi_min) / (phi_max - phi_min)

    probe_desc = (
        f'Probe {row} '
        f'(r={r_val:.3f}, '
        r'$\theta$=' f'{theta_deg:.1f}' r'$^\circ$, '
        r'$\phi$=' f'{phi_deg:.1f}' r'$^\circ$)'
        '\n'
        r'($r_{nrm}$=' f'{r_nrm:.3f}, '
        r'$\theta_{nrm}$=' f'{theta_nrm:.3f}, '
        r'$\phi_{nrm}$=' f'{phi_nrm:.3f})'
    )

    for col in range(3):

        for key in ('FI', 'FC', 'FA', 'FV', 'FP'):

            comp = forces[key][col]

            if comp is None:
                continue

            ax[row, col].plot(
                time, comp,
                color=COLORS[key],
                lw=1.2,
                label=LABELS[key]
            )

        ax[row, col].plot(
            time, forces['BAL'][col],
            color=COLORS['BAL'],
            lw=1.5,
            label=BAL_LABELS[col]
        )

        ax[row, col].axhline(
            0,
            color='black',
            lw=0.8,
            linestyle='--'
        )

        if row == 0:
            ax[row, col].set_title(
                COMPONENT_TITLES[col],
                fontsize=16
            )

        ax[row, col].legend(
            loc='upper right',
            fontsize=10,
            frameon=False
        )

    ax[row, 0].set_ylabel(
        f'{probe_desc}\n\n' r'$\nabla \times$ Force',
        fontsize=12
    )


# ============================================================
# GENERAL FORMATTING
# ============================================================

for a in ax.flat:

    a.tick_params(
        axis='both',
        labelsize=12
    )

for a in ax[-1, :]:

    a.set_xlabel(
        r'$t$',
        fontsize=14
    )

    # Uncomment if you want grid lines
    # a.grid(True, alpha=0.3)


plt.tight_layout()

plt.savefig(
    output_file,
    dpi=300,
    bbox_inches='tight'
)

plt.close()

print(f'{output_file} Done')
