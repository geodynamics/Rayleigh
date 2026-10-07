"""
Quantity codes from src/Diagnostics/momentum_equation_codes.F.

QUANTITY_CODES covers the non-magnetic forces.
MAGNETIC_QUANTITY_CODES (j_cross_b_*) requires magnetism=.true. and a
B-field reconstruction (see magnetic_field.py).
"""

QUANTITY_CODES = {
    1201: 'v_grad_v_r', 1202: 'v_grad_v_theta', 1203: 'v_grad_v_phi',
    1216: 'buoyancy_force',
    1219: 'Coriolis_Force_r', 1220: 'Coriolis_Force_theta', 1221: 'Coriolis_Force_phi',
    1228: 'viscous_Force_r', 1229: 'viscous_Force_theta', 1230: 'viscous_Force_phi',
    1237: 'pressure_Force_r', 1238: 'pressure_Force_theta', 1239: 'pressure_Force_phi',
}

MAGNETIC_QUANTITY_CODES = {
    1248: 'j_cross_b_r', 1249: 'j_cross_b_theta', 1250: 'j_cross_b_phi',
}
