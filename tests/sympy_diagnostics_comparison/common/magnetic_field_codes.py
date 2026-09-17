"""
Magnetic field quantity codes from src/Diagnostics/magnetic_field_codes.F.
"""

QUANTITY_CODES = {
    801: 'b_r', 802: 'b_theta', 803: 'b_phi',

    810: 'db_r_dr', 811: 'db_theta_dr', 812: 'db_phi_dr',
    819: 'db_r_dt', 820: 'db_theta_dt', 821: 'db_phi_dt',
    828: 'db_r_dp', 829: 'db_theta_dp', 830: 'db_phi_dp',
    837: 'db_r_dtr', 838: 'db_theta_dtr', 839: 'db_phi_dtr',
    846: 'db_r_dprs', 847: 'db_theta_dprs', 848: 'db_phi_dprs',

    855: 'db_r_d2r', 856: 'db_theta_d2r', 857: 'db_phi_d2r',
    864: 'db_r_d2t', 865: 'db_theta_d2t', 866: 'db_phi_d2t',
    873: 'db_r_d2p', 874: 'db_theta_d2p', 875: 'db_phi_d2p',
    882: 'db_r_d2rt', 883: 'db_theta_d2rt', 884: 'db_phi_d2rt',
    891: 'db_r_d2rp', 892: 'db_theta_d2rp', 893: 'db_phi_d2rp',
    900: 'db_r_d2tp', 901: 'db_theta_d2tp', 902: 'db_phi_d2tp',
}
