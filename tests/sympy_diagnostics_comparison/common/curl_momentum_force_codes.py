"""Tested quantity codes from src/Diagnostics/curl_momentum_equation_codes.F,

QUANTITY_CODES covers the non-magnetic forces.
MAGNETIC_QUANTITY_CODES (curl_j_cross_b_*) requires magnetism=.true. and a
B-field reconstruction (see magnetic_field.py).
"""

QUANTITY_CODES = {
    1301: 'curl_v_grad_v_r',
    1302: 'curl_v_grad_v_theta',
    1303: 'curl_v_grad_v_phi',

    1316: 'curl_buoyancy_force_theta',
    1317: 'curl_buoyancy_force_phi',

    1322: 'curl_coriolis_force_r',
    1323: 'curl_coriolis_force_theta',
    1324: 'curl_coriolis_force_phi',

    1331: 'curl_viscous_force_r',
    1332: 'curl_viscous_force_theta',
    1333: 'curl_viscous_force_phi',

    1341: 'curl_pressure_force_theta',
    1342: 'curl_pressure_force_phi',
}

# Curl of the fluctuating (non-axisymmetric) and mean (axisymmetric) viscous
# forces; only meaningful for a velocity field with both m=0 and m!=0 content
# (see free_slip).
MEAN_FLUCTUATING_CODES = {
    1334: 'curl_viscous_pforce_r',
    1335: 'curl_viscous_pforce_theta',
    1336: 'curl_viscous_pforce_phi',

    1337: 'curl_viscous_mforce_r',
    1338: 'curl_viscous_mforce_theta',
    1339: 'curl_viscous_mforce_phi',
}

MAGNETIC_QUANTITY_CODES = {
    1349: 'curl_j_cross_b_r',
    1350: 'curl_j_cross_b_theta',
    1351: 'curl_j_cross_b_phi',
}
