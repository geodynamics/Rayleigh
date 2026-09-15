#!/home/pgrad1/2831664s/anaconda3/bin/python

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'common'))

from meridional_animation import (
    animate, code_row, field, make_frame_figure, plot_2row_frame,
    Meridional_Slices, build_file_list, plt,
    start_file, end_file,
    VISCOUS_FORCE_ROW, CURL_VISCOUS_FORCE_ROW,
)


# ============================================================
# QUANTITY CODE GROUPS (r, theta, phi) -- MHD-only, on top of the
# VISCOUS_FORCE/CURL_VISCOUS_FORCE groups shared via meridional_animation
# ============================================================

J_CROSS_B      = (1248, 1249, 1250)
CURL_J_CROSS_B = (1369, 1370, 1371)

J_CROSS_B_ROW      = code_row(J_CROSS_B)
CURL_J_CROSS_B_ROW = code_row(CURL_J_CROSS_B)


# ============================================================
# FORCE-BALANCE RESIDUALS
# (mirrors the balance_r/t/p sums in pp.py and pp_curl.py -- the
# forces should cancel, so a nonzero residual reveals imbalance/noise)
# ============================================================

def force_residual(ms, tindex, col):
    F = lambda code: field(ms, tindex, code)

    if col == 0:    # r: Inertia - Coriolis - Pressure - Buoyancy - Lorentz - Viscous
        return F(1201) - F(1219) - F(1237) - F(1216) - F(1248) - F(1228)
    elif col == 1:  # theta
        return F(1202) - F(1220) - F(1238) - F(1249) - F(1229)
    else:           # phi
        return F(1203) - F(1221) - F(1239) - F(1250) - F(1230)


def curl_force_residual(ms, tindex, col):
    F = lambda code: field(ms, tindex, code)

    if col == 0:    # r: curl(Inertia - Coriolis - Lorentz - Viscous)
        return F(1301) - F(1327) - F(1369) - F(1339)
    elif col == 1:  # theta
        return F(1302) - F(1328) - F(1319) - F(1370) - F(1358) - F(1340)
    else:           # phi
        return F(1303) - F(1329) - F(1320) - F(1371) - F(1359) - F(1341)


# ============================================================
# READ FIRST MERIDIONAL SLICE FILE AND PLOT A PREVIEW
# (force-balance residuals, as a sanity check of the grid/layout)
# ============================================================

if __name__ == '__main__':

    output_file = 'Meridional_Slices_frame.pdf'

    files = build_file_list(
        start_file,
        end_file,
        path='Meridional_Slices'
    )

    print(f'Number of Meridional_Slices files: {len(files)}')

    ms = Meridional_Slices(files[0], path='')

    fig, axes = make_frame_figure()

    plot_2row_frame(
        fig, axes, ms, tindex=0,
        row_funcs=(force_residual, curl_force_residual),
        row_labels=('Force residual', 'Curl-of-force residual')
    )

    plt.savefig(
        output_file,
        dpi=200,
        bbox_inches='tight'
    )

    plt.close()

    print(f'{output_file} Done')
