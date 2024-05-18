#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ====================================
# @File name        : Makela08_alloc_parameter.py
# @First created    : 2024/5/15 17:38
# @Author           : Yuzhi Zhu
# @E-mail           : edwardmashed@gmail.com
# ====================================

params_makela08_spruce = {
    "AvgLongevity_foliage":     8,      # Mean lifetime (yr)
    "AvgLongevity_wood":        33.3,   # Mean lifetime (yr)
    "AvgLongevity_root":        1.25,   # Mean lifetime (yr)

    "alpha_w":                  0.4,    # Sapwood weight per unit foliage and pipe length
    "c_H":                      3400,   # ‘Steady-state’ pipe length coefficient

    "Kr":                      2000,   # Amount of roots capturing 50% of available N (kg/ha)
    "Kf":                      8000,   # Amount of foliage capturing 50% of maximum C gain (kg/ha)

    "NResorbFrac_foliage":      0.3,    # Proportion N recycled
    "NResorbFrac_wood":         0.3,    # Proportion N recycled
    "NResorbFrac_root":         0.3,    # Proportion N recycled
    "Resp_Nspecific":           16,     # Specific rate of maintenance respiration (kg−1 N yr−1)
    "CtoDM_frac":               1.54,   # Growth efficiency, Yg (kg DW kg−1 C)
    "NrNf_ratio":               1,      # Ratio of fine root [N] to foliage [N]
    "NwNf_ratio":               0.07,   # Ratio of sapwood [N] to foliage [N]

    "Nconc_foliage_structural": 0.008,  # Structural (non-photosynthetic) N concentration in foliage (kg N/kg DW)
    "Nconc_ref":                0.002,  # Reference photosynthetic N at sigma_fM = sigma_fM0/2 (kg N/kg DW)
    "Photosyn_Nsat":            4,      # N-saturated rate of photosynthesis (kg C / kg DW-1 yr-1)
}

params_makela08_pine = {
    "AvgLongevity_foliage":     3.3,
    "AvgLongevity_wood":        40,
    "AvgLongevity_root":        1.25,

    "alpha_w":                  0.8,
    "c_H":                      2800,

    "Kr":                      2000,
    "Kf":                      2500,

    "NResorbFrac_foliage":      0.3,
    "NResorbFrac_wood":         0.3,
    "NResorbFrac_root":         0.3,
    "Resp_Nspecific":           16,
    "CtoDM_frac":               1.54,
    "NrNf_ratio":               1,
    "NwNf_ratio":               0.07,

    "Nconc_foliage_structural": 0.009,
    "Nconc_ref":                0.002,
    "Photosyn_Nsat":            8,
}


def get_param_dict(dict_name="params_makela08_spruce"):
    dict_name_list = ["params_makela08_spruce", "params_makela08_pine"]
    if str(dict_name) in dict_name_list:
        return eval(dict_name)
    else:
        raise ValueError(f"Parameter set '{dict_name}' not found.")

