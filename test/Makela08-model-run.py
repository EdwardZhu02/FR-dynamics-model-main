#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ====================================
# @File name        : Makela08-model-run.py
# @First created    : 2024/5/15 20:59
# @Author           : Yuzhi Zhu
# @E-mail           : edwardmashed@gmail.com
# ====================================

from modules.cn_allocation import Makela08_alloc_parameter
import numpy as np
import sympy as sp

params_dict = Makela08_alloc_parameter.get_param_dict(dict_name="params_makela08_pine")


def matrix_generator(x_start, x_end, y_start, y_end, step):
    """
    Generate 2-D matrix containing Light and N availability gradients

    X: (sigma_rM) N availability (kg N / kg FR / yr)
    `Nup_max_specific`

    Y: (sigma_fM) Light saturated foliar-specific photosynthetic rate (kg C / kg foliage / yr)
    `Photosyn_lightsat`

    :return: target matrix
    """
    num_x = int((x_end - x_start) / step) + 1
    num_y = int((y_end - y_start) / step) + 1
    matrix_n = np.zeros((num_y, num_x))
    matrix_p = np.zeros((num_y, num_x))

    for _i in range(num_y):
        for _j in range(num_x):
            matrix_n[_i][_j] = x_start + _j * step
            matrix_p[_i][_j] = y_start + _i * step
    return matrix_n, matrix_p


gradient_matrix_n, gradient_matrix_p = matrix_generator(0, 10, 0, 10, 0.01)

# Load parameters into separate values
AvgLongevity_foliage = params_dict["AvgLongevity_foliage"]
AvgLongevity_wood = params_dict["AvgLongevity_wood"]
AvgLongevity_root = params_dict["AvgLongevity_root"]

NResorbFrac_foliage = params_dict["NResorbFrac_foliage"]
NResorbFrac_wood = params_dict["NResorbFrac_wood"]
NResorbFrac_root = params_dict["NResorbFrac_root"]

alpha_w = params_dict["alpha_w"]
c_H = params_dict["c_H"]
NrNf_ratio = params_dict["NrNf_ratio"]
NwNf_ratio = params_dict["NwNf_ratio"]
Resp_Nspecific = params_dict["Resp_Nspecific"]
CtoDM_frac = params_dict["CtoDM_frac"]
Kr = params_dict["Kr"]
Kf = params_dict["Kf"]

# Define symbolic variables
Nconc_foliage, psi_r, Photosyn_lightsat, Nup_max_specific = sp.symbols(
    'Nconc_foliage psi_r Photosyn_lightsat Nup_max_specific'
)

exp_beta1 = CtoDM_frac * Photosyn_lightsat * Kf / (1 / AvgLongevity_root + Resp_Nspecific * Nconc_foliage * NrNf_ratio)

exp_beta2 = (1 / AvgLongevity_foliage + Nconc_foliage * (
            alpha_w * c_H / AvgLongevity_wood +
            CtoDM_frac * Resp_Nspecific * (1 + NwNf_ratio * alpha_w * c_H * Nconc_foliage)
    )) / (1 / AvgLongevity_root + CtoDM_frac * Resp_Nspecific * Nconc_foliage * NrNf_ratio)

exp_beta3 = (Nup_max_specific * Kr) / (
            Nconc_foliage * (1 - NResorbFrac_root) * NrNf_ratio / AvgLongevity_root
    )

exp_beta4 = (
                    (1 - NResorbFrac_foliage) / AvgLongevity_foliage +
                    (1 - NResorbFrac_wood) * NwNf_ratio * alpha_w * c_H * Nconc_foliage / AvgLongevity_wood
            ) / (
                    (1 - NResorbFrac_root) * NrNf_ratio / AvgLongevity_root
            )

exp_a1 = (exp_beta1 - exp_beta3 + Kr - Kf * (exp_beta2 + exp_beta4)) / -Kf
exp_a2 = (exp_beta1 * exp_beta4 - exp_beta2 * exp_beta3 + Kr * (exp_beta2 + exp_beta4) - Kf * exp_beta2 * exp_beta4) / -Kf
exp_a3 = Kr * exp_beta2 * exp_beta4 / -Kf

exp_cubiceq_psi_r = psi_r**3 + exp_a1*psi_r**2 + exp_a2*psi_r + exp_a3

# Symbolic solutions, containing Nconc_foliage, psi_r, Photosyn_lightsat, Nup_max_specific
# Photosyn_lightsat, Nup_max_specific is evaluated later so the final solution is psi_r given a Nconc_foliage.
exp_cubiceq_psi_r_solutions = sp.solve(exp_cubiceq_psi_r, psi_r)

# Enumerate the 2-D matrix
try:
    for i in range(len(gradient_matrix_n)):  # or gradient_matrix_p
        for j in range(len(gradient_matrix_n[i])):
            _Nup_max_specific = gradient_matrix_n[i, j]
            _Photosyn_lightsat = gradient_matrix_p[i, j]

            # Symbolic solutions of psi_r containing only Nconc_foliage as variable
            psi_r_solution_list = [
                i.evalf(subs={Nup_max_specific: _Nup_max_specific, Photosyn_lightsat: _Photosyn_lightsat}) for i in
                exp_cubiceq_psi_r_solutions]

            raise StopIteration
except StopIteration:
    pass

