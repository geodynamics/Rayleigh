"""
Velocity quantity codes from src/Diagnostics/velocity_field_codes.F.
"""

QUANTITY_CODES = {
    1: 'v_r', 2: 'v_theta', 3: 'v_phi',

    10: 'dv_r_dr', 11: 'dv_theta_dr', 12: 'dv_phi_dr',
    19: 'dv_r_dt', 20: 'dv_theta_dt', 21: 'dv_phi_dt',
    28: 'dv_r_dp', 29: 'dv_theta_dp', 30: 'dv_phi_dp',
    37: 'dv_r_dtr', 38: 'dv_theta_dtr', 39: 'dv_phi_dtr',
    46: 'dv_r_dprs', 47: 'dv_theta_dprs', 48: 'dv_phi_dprs',

    55: 'dv_r_d2r', 56: 'dv_theta_d2r', 57: 'dv_phi_d2r',
    64: 'dv_r_d2t', 65: 'dv_theta_d2t', 66: 'dv_phi_d2t',
    73: 'dv_r_d2p', 74: 'dv_theta_d2p', 75: 'dv_phi_d2p',
    82: 'dv_r_d2rt', 83: 'dv_theta_d2rt', 84: 'dv_phi_d2rt',
    91: 'dv_r_d2rp', 92: 'dv_theta_d2rp', 93: 'dv_phi_d2rp',
    100: 'dv_r_d2tp', 101: 'dv_theta_d2tp', 102: 'dv_phi_d2tp',
}
