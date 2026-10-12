"""Tested quantity codes from src/Diagnostics/curl_momentum_equation_codes.F,

QUANTITY_CODES covers the non-magnetic forces.
MAGNETIC_QUANTITY_CODES (curl_j_cross_b_*) requires magnetism=.true. and a
B-field reconstruction (see magnetic_field.py).
"""

QUANTITY_CODES = {
    1301: 'curl_v_grad_v_r',
    1302: 'curl_v_grad_v_theta',
    1303: 'curl_v_grad_v_phi',

    1319: 'curl_buoyancy_force_theta',
    1320: 'curl_buoyancy_force_phi',

    1327: 'curl_coriolis_force_r',
    1328: 'curl_coriolis_force_theta',
    1329: 'curl_coriolis_force_phi',

    1339: 'curl_viscous_force_r',
    1340: 'curl_viscous_force_theta',
    1341: 'curl_viscous_force_phi',

    1358: 'curl_pressure_force_theta',
    1359: 'curl_pressure_force_phi',
}

# Curl of the fluctuating (non-axisymmetric) and mean (axisymmetric) viscous
# forces; only meaningful for a velocity field with both m=0 and m!=0 content
# (see free_slip).
MEAN_FLUCTUATING_CODES = {
    1345: 'curl_viscous_pforce_r',
    1346: 'curl_viscous_pforce_theta',
    1347: 'curl_viscous_pforce_phi',

    1351: 'curl_viscous_mforce_r',
    1352: 'curl_viscous_mforce_theta',
    1353: 'curl_viscous_mforce_phi',
}

# Squares and magnitude of the curl of the (full) viscous force.
CURL_VISCOUS_DERIVED_CODES = {
    1342: 'curl_viscous_force_r_squared',
    1343: 'curl_viscous_force_theta_squared',
    1344: 'curl_viscous_force_phi_squared',
    1391: 'curl_viscous_force_abs',
}

MAGNETIC_QUANTITY_CODES = {
    1369: 'curl_j_cross_b_r',
    1370: 'curl_j_cross_b_theta',
    1371: 'curl_j_cross_b_phi',
}
