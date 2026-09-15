#!/home/pgrad1/2831664s/anaconda3/bin/python

from rayleigh_diagnostics import Meridional_Slices, GridInfo, Point_Probes, plot_azav, build_file_list, get_lims

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

from tqdm import tqdm


# ============================================================
# SETTINGS
# ============================================================

start_file = 1
end_file   = 10000000

phi_index = 0  # Only one Meridional Slice (phi=0) was output

COMPONENT_TITLES = (r'$r$ component', r'$\theta$ component', r'$\phi$ component')


# ============================================================
# QUANTITY CODE GROUPS (r, theta, phi) common to every case
# (case-specific groups, e.g. J_CROSS_B/CURL_J_CROSS_B for MHD,
# are defined in each case's own plot_meridional.py)
# ============================================================

VISCOUS_FORCE      = (1228, 1229, 1230)
CURL_VISCOUS_FORCE = (1339, 1340, 1341)


def field(ms, tindex, code):
    return ms.vals[phi_index, :, :, ms.lut[code], tindex]


def code_row(codes):
    """Build a row function (ms, tindex, col) -> field from a
    (r, theta, phi) quantity-code tuple.
    """
    def _row(ms, tindex, col):
        return field(ms, tindex, codes[col])
    return _row


VISCOUS_FORCE_ROW      = code_row(VISCOUS_FORCE)
CURL_VISCOUS_FORCE_ROW = code_row(CURL_VISCOUS_FORCE)


# ============================================================
# POINT PROBE LOCATIONS (read from the point probe output itself,
# so the markers reflect where Rayleigh actually sampled -- snapped
# to the nearest grid point -- rather than the requested main_input
# values) that fall on the phi=0 meridional slice, in plot_azav's
# (x, y) plane
# ============================================================

def point_probe_coords(path='Point_Probes'):
    """(x, y) plot-plane coordinates of every point probe that lies on
    the phi=0 meridional slice. Only the first Point_Probes file is
    read, since probe locations don't change over time.
    """

    files = build_file_list(start_file, end_file, path=path)

    if not files:
        return []

    pp = Point_Probes(files[0], path='')

    on_slice = np.isclose(pp.phi, 0.0, atol=1e-12) | np.isclose(pp.phi, 2*np.pi, atol=1e-12)

    if not on_slice.any():
        return []  # no probes on the phi=0 slice

    r_max = GridInfo(path='./').radius.max()

    coords = []

    for r_val in pp.radius:

        r_frac = r_val / r_max

        for costheta in pp.costheta:

            sintheta = (1.0 - costheta**2)**0.5

            coords.append((r_frac * sintheta, r_frac * costheta))

    return coords


POINT_PROBE_COORDS = point_probe_coords()


# ============================================================
# FIGURE/AXES FOR A 2-ROW x 3-COLUMN GRID OF MERIDIONAL PANELS
# (margins/spacing tuned to minimize whitespace around and between
# the half-disk panels, which plot_azav draws with axis('equal') --
# that letterboxes each panel inside its cell, so the usual
# tight_layout/default subplot spacing leaves large gaps)
# ============================================================

def make_frame_figure(figsize=(9, 9.5)):
    fig, axes = plt.subplots(2, 3, figsize=figsize, squeeze=False)
    fig.subplots_adjust(
        left=0.07, right=0.98, top=0.93, bottom=0.01,
        wspace=-0.2, hspace=-0.08
    )
    return fig, axes


# ============================================================
# FUNCTION TO PLOT ONE FRAME OF A 2-ROW x 3-COLUMN GRID
# ============================================================

def plot_2row_frame(fig, axes, ms, tindex, row_funcs, row_labels,
                     cbar_cache=None, mycmap='seismic',
                     boundsfactor=4.5, boundstype='rms', bounds=None):
    """Plot a 2-row (row_funcs) x 3-column (r, theta, phi) grid of
    meridional-slice fields for a single time snapshot. Reused as-is
    to draw each frame of an animation.

    row_funcs  : length-2 sequence of functions (ms, tindex, col) ->
                 2D field, e.g. code_row(...) or a residual function.
    row_labels : length-2 sequence of row labels.
    cbar_cache : optional dict, reused across calls, used to remove
                 each panel's previous colorbar before drawing a new
                 one -- required because plot_azav's colorbar axes
                 are not cleared automatically by ax.clear(), so
                 without this they would pile up frame after frame.
    bounds     : optional {(row, col): (mini, maxi)} dict of fixed
                 color limits, e.g. from animate()'s first frame. If a
                 panel is missing from bounds (or bounds is None),
                 that panel's limits are auto-computed from its own
                 field each call via boundsfactor/boundstype.
    """

    radius   = ms.radius
    costheta = ms.costheta
    sintheta = ms.sintheta

    for row, row_func in enumerate(row_funcs):

        for col in range(3):

            ax = axes[row, col]
            key = (row, col)

            if cbar_cache is not None and key in cbar_cache:
                cbar_cache[key].remove()

            ax.clear()

            field_vals = row_func(ms, tindex, col)

            mini, maxi = bounds.get(key, (-1, -1)) if bounds is not None else (-1, -1)

            img = plot_azav(
                fig, ax, field_vals, radius, costheta, sintheta,
                mini=mini, maxi=maxi,
                mycmap=mycmap,
                boundsfactor=boundsfactor,
                boundstype=boundstype,
                units='',
                fontsize=8
            )

            if cbar_cache is not None:
                cbar_cache[key] = img.colorbar

            if POINT_PROBE_COORDS:
                xs, ys = zip(*POINT_PROBE_COORDS)
                ax.plot(
                    xs, ys,
                    marker='x', color='black', linestyle='none',
                    markersize=6, markeredgewidth=1.5, zorder=5
                )

            if row == 0:
                ax.set_title(COMPONENT_TITLES[col], fontsize=14)

        axes[row, 0].text(
            -0.15, 0.5, row_labels[row],
            transform=axes[row, 0].transAxes,
            fontsize=14,
            ha='center', va='center',
            rotation=90
        )

    time_val = ms.time[tindex]
    iter_val = ms.iters[tindex]

    fig.suptitle(
        f't = {time_val:.5f}    iteration = {iter_val}',
        fontsize=16
    )


# ============================================================
# DRIVE A FULL ANIMATION FROM A PAIR OF ROW FUNCTIONS
# ============================================================

def animate(output_file, row_funcs, row_labels, fps=15, dpi=150, figsize=(9, 9.5),
            frame_stride=1, mycmap='seismic', boundsfactor=4.5, boundstype='rms'):
    """Read every available Meridional_Slices file and write an mp4
    animation of the 2x3 (row_funcs x r/theta/phi) grid, one frame per
    time snapshot, over the second half of the available snapshots
    (skipping the first half lets the early transient settle out of
    both the animated range and the fixed color limits below). Color
    limits are fixed for the whole animation, computed once from the
    first animated frame, rather than auto-rescaling every frame.

    frame_stride : keep only every frame_stride-th snapshot of the
                   (already second-half-only) range, e.g. to thin out
                   a long, slowly-varying run.
    """

    files = build_file_list(
        start_file,
        end_file,
        path='Meridional_Slices'
    )

    nfiles = len(files)

    print(f'Number of Meridional_Slices files: {nfiles}')

    frames = []  # list of (ms, local time-index) pairs, one per animation frame

    for f in tqdm(files, desc='Reading Meridional_Slices files'):

        ms = Meridional_Slices(f, path='')

        for tindex in range(ms.niter):
            frames.append((ms, tindex))

    nsnapshots = len(frames)

    half_start = int(nsnapshots / 2 + 1) - 1  # e.g. x available -> keep int(x/2+1)..x (1-indexed)
    frames = frames[half_start:]

    frames = frames[::frame_stride]

    nframes = len(frames)

    print(f'Number of animation frames: {nframes} (of {nsnapshots} available snapshots)')

    fig, axes = make_frame_figure(figsize=figsize)

    # Fixed color limits for the whole animation, computed once from
    # the first animated frame (the early transient is already
    # excluded by the second-half selection above).
    bounds = {}

    if frames:
        ms0, tindex0 = frames[0]

        for row, row_func in enumerate(row_funcs):
            for col in range(3):
                field_vals = row_func(ms0, tindex0, col)
                bounds[(row, col)] = (
                    get_lims(field_vals, boundsfactor=boundsfactor, boundstype=boundstype, themin=True),
                    get_lims(field_vals, boundsfactor=boundsfactor, boundstype=boundstype, themin=False)
                )

    cbar_cache = {}

    def update(frame_index):

        ms, tindex = frames[frame_index]

        plot_2row_frame(
            fig, axes, ms, tindex, row_funcs, row_labels,
            cbar_cache=cbar_cache, mycmap=mycmap,
            boundsfactor=boundsfactor, boundstype=boundstype,
            bounds=bounds
        )

        return axes.flat

    anim = animation.FuncAnimation(
        fig,
        update,
        frames=nframes,
        blit=False
    )

    writer = animation.FFMpegWriter(fps=fps, bitrate=-1)

    with tqdm(total=nframes, desc=f'Writing {output_file}') as pbar:

        anim.save(
            output_file, writer=writer, dpi=dpi,
            progress_callback=lambda i, n: pbar.update(1)
        )

    plt.close()

    print(f'{output_file} Done')
