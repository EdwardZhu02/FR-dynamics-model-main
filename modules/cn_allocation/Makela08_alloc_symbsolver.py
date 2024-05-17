#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ====================================
# @File name        : Makela08_alloc_symbsolver.py
# @First created    : 2024/5/16 16:39
# @Author           : Yuzhi Zhu
# @E-mail           : edwardmashed@gmail.com
# ====================================

# Input: psi_r = f(Nconc_foliage) - three possible roots for a given Navail and Photosyn conditions
# Output: G = f(Nconc_foliage)
import sympy as sp
import numpy as np
from scipy.optimize import minimize_scalar


class DryMassFoliageSolver:
    def __init__(self, Nup_max_specific, Photosyn_lightsat, params_dict):
        # Fixed parameters
        self.alpha_w = params_dict["alpha_w"]
        self.c_H = params_dict["c_H"]
        self.NrNf_ratio = params_dict["NrNf_ratio"]
        self.NwNf_ratio = params_dict["NwNf_ratio"]
        self.Resp_Nspecific = params_dict["Resp_Nspecific"]
        self.CtoDM_frac = params_dict["CtoDM_frac"]
        self.Kr = params_dict["Kr"]
        self.Kf = params_dict["Kf"]
        self.AvgLongevity_foliage = params_dict["AvgLongevity_foliage"]
        self.AvgLongevity_wood = params_dict["AvgLongevity_wood"]
        self.AvgLongevity_root = params_dict["AvgLongevity_root"]
        self.NResorbFrac_foliage = params_dict["NResorbFrac_foliage"]
        self.NResorbFrac_wood = params_dict["NResorbFrac_wood"]
        self.NResorbFrac_root = params_dict["NResorbFrac_root"]
        # Specified by environmental conditions
        self.Photosyn_lightsat = Photosyn_lightsat
        self.Nup_max_specific = Nup_max_specific
        # Variable used for symbolic calculation
        self.Nconc_foliage = sp.symbols("Nconc_foliage")

    def solve_carbon(self, symb_psi_r):
        """
        Eqn. S3, steady state W_f calculation: W_f(C) = f(Nconc_foliage)
        :param symb_psi_r: psi_r = f(Nconc_foliage)
        :return:
        """
        symb_beta1 = (self.CtoDM_frac * self.Photosyn_lightsat * self.Kf /
                      (1 / self.AvgLongevity_root + self.Resp_Nspecific * self.Nconc_foliage * self.NrNf_ratio))

        symb_beta2 = ((1 / self.AvgLongevity_foliage + self.Nconc_foliage * (
                self.alpha_w * self.c_H / self.AvgLongevity_wood +
                self.CtoDM_frac * self.Resp_Nspecific * (
                        1 + self.NwNf_ratio * self.alpha_w * self.c_H * self.Nconc_foliage)
        )) / (
                              1 / self.AvgLongevity_root + self.CtoDM_frac * self.Resp_Nspecific * self.Nconc_foliage * self.NrNf_ratio))

        symb_DM_foliage_C = symb_beta1 / (symb_beta2 + symb_psi_r) - self.Kf
        return [symb_DM_foliage_C, symb_beta1, symb_beta2]

    def solve_nitrogen(self, symb_psi_r):
        """
        Eqn. S7,  steady state W_f calculation: W_f(N) = f(Nconc_foliage)
        :param symb_psi_r: psi_r = f(Nconc_foliage)
        :return:
        """
        Nup_max = self.Nup_max_specific * self.Kr

        symb_beta3 = Nup_max / (
                self.Nconc_foliage * (1 - self.NResorbFrac_root) * self.NrNf_ratio / self.AvgLongevity_root)

        symb_beta4 = (
                ((1 - self.NResorbFrac_foliage) / self.AvgLongevity_foliage +
                 (
                         1 - self.NResorbFrac_wood) * self.NwNf_ratio * self.alpha_w * self.c_H * self.Nconc_foliage / self.AvgLongevity_wood
                 ) / ((1 - self.NResorbFrac_root) * self.NrNf_ratio / self.AvgLongevity_root)
        )

        symb_DM_foliage_N = symb_beta3 / (symb_beta4 + symb_psi_r) - self.Kr / symb_psi_r
        return [symb_DM_foliage_N, symb_beta3, symb_beta4]


class BiomassProductionSolver:
    def __init__(self, params_dict):
        self.Nconc_foliage = sp.symbols("Nconc_foliage")
        self.AvgLongevity_foliage = params_dict["AvgLongevity_foliage"]
        self.AvgLongevity_wood = params_dict["AvgLongevity_wood"]
        self.AvgLongevity_root = params_dict["AvgLongevity_root"]
        self.alpha_w = params_dict["alpha_w"]
        self.c_H = params_dict["c_H"]

    def solve_total_biomass_production(self, symb_DM_foliage, symb_psi_r):
        """

        :param symb_DM_foliage: f(Nconc_foliage)
        :param symb_psi_r: f(Nconc_foliage)
        :return:
        """
        symb_DM_production = symb_DM_foliage * (
                1 / self.AvgLongevity_foliage + symb_psi_r / self.AvgLongevity_root +
                self.alpha_w * self.c_H * self.Nconc_foliage / self.AvgLongevity_wood
        )
        return symb_DM_production


class BiomassProductionOptimizer:
    def __init__(self, symb_G_C, symb_G_N):
        """
        Init G functions in which symb_DM_foliage calculated using either C or N balance
        :param symb_G_C: symb_DM_foliage: f(Nconc_foliage) carbon balance
        :param symb_G_N: symb_DM_foliage: f(Nconc_foliage) nitrogen balance
        """
        self._symb_G_C = symb_G_C
        self._symb_G_N = symb_G_N

    def optimize_total_biomass_production(self, range_lower=0, range_upper=10, method="C"):
        """
        Optimize total G by conducting numeric search of [N]_f (Nconc_foliage)
        :param range_lower: upper limit of optimization range (should be positive)
        :param range_upper: lower limit of optimization range (should be positive)
        :param method: steady state W_f calculated from C or N balance (C/N)
        :return:
        """
        # Convert the symbolic equation to a numerical function
        Nconc_foliage = sp.symbols("Nconc_foliage")
        if method == "C":
            f_numeric = sp.lambdify(Nconc_foliage, self._symb_G_C, 'numpy')
        else:  # method == "N"
            f_numeric = sp.lambdify(Nconc_foliage, self._symb_G_N, 'numpy')

        def func_to_minimize(x_val):
            return -f_numeric(x_val)

        # Use scipy's minimize_scalar function to find the maximum
        f_max_result = minimize_scalar(func_to_minimize,
                                       bounds=(float(range_lower), float(range_upper)), method='bounded')
        # Get the maximum value and the corresponding x value
        Nconc_foliage_maxG = f_max_result.x
        maxG_value = -f_max_result.fun

        return Nconc_foliage_maxG, maxG_value


