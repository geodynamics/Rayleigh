"""
Quantity codes from src/Diagnostics/momentum_equation_codes.F.
"""

QUANTITY_CODES = {
    1201: 'v_grad_v_r', 1202: 'v_grad_v_theta', 1203: 'v_grad_v_phi',
    1216: 'buoyancy_force',
    1219: 'Coriolis_Force_r', 1220: 'Coriolis_Force_theta', 1221: 'Coriolis_Force_phi',
    1228: 'viscous_Force_r', 1229: 'viscous_Force_theta', 1230: 'viscous_Force_phi',
    1237: 'pressure_Force_r', 1238: 'pressure_Force_theta', 1239: 'pressure_Force_phi',
}
