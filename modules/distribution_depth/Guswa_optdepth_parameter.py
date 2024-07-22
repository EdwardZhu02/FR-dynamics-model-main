#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ====================================
# @File name        : Guswa_optdepth_parameter.py
# @First created    : 2024/6/4 17:52
# @Author           : Yuzhi Zhu
# @E-mail           : edwardmashed@gmail.com
# ====================================

params_guswa_yang2016 = {

    # TODO: derive `RootResp_20` from biome-specific lookup table
    # TODO: derive other variables

    "RootResp_20": 1.0,             # gamma_r20 | Root respiration rate @ 20 degrees, g C / (g root * day)
    "RootLengthDensity": 0.1,       # RLD, cm root / cm3 soil
    "SpecificRootLength": 1500,     # SRL, cm root / g root
    "WUE_photosynthesis": 0.5,      # WUE, fraction (need derivation)
    "FractionGrowingSeason": 0.5,   # f_GS, fraction (need derivation)
    "RainfallDepth_singevent_mean": 1.0,       # alpha (need derivation)
    "WaterContent_Avail": 0.5,      # theta (need derivation)
}


def get_param_dict(dict_name="params_default"):
    dict_name_list = ["params_default"]
    if str(dict_name) in dict_name_list:
        return eval(dict_name)
    else:
        raise ValueError(f"Parameter set '{dict_name}' not found.")
