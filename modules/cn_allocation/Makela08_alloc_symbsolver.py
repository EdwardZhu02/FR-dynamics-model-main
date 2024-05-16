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
                    self.CtoDM_frac * self.Resp_Nspecific * (1 + self.NwNf_ratio * self.alpha_w * self.c_H * self.Nconc_foliage)
                )) / (1 / self.AvgLongevity_root + self.CtoDM_frac * self.Resp_Nspecific * self.Nconc_foliage * self.NrNf_ratio))

        symb_DM_foliage_C = symb_beta1 / (symb_beta2 + symb_psi_r) - self.Kf
        return symb_DM_foliage_C, symb_beta1, symb_beta2

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
             (1 - self.NResorbFrac_wood) * self.NwNf_ratio * self.alpha_w * self.c_H * self.Nconc_foliage / self.AvgLongevity_wood
             ) / ((1 - self.NResorbFrac_root) * self.NrNf_ratio / self.AvgLongevity_root)
        )

        symb_DM_foliage_N = symb_beta3 / (symb_beta4 + symb_psi_r) - self.Kr / symb_psi_r
        return symb_DM_foliage_N, symb_beta3, symb_beta4


class BiomassProductionSolver:
    def __init__(self, psi_r, params_dict):
        self.Nconc_foliage = sp.symbols("Nconc_foliage")
        self.psi_r = psi_r
        self.AvgLongevity_foliage = params_dict["AvgLongevity_foliage"]
        self.AvgLongevity_wood = params_dict["AvgLongevity_wood"]
        self.AvgLongevity_root = params_dict["AvgLongevity_root"]
        self.alpha_w = params_dict["alpha_w"]
        self.c_H = params_dict["c_H"]

    def solve_total_biomass_production(self, symb_DM_foliage):
        symb_DM_production = symb_DM_foliage * (
            1 / self.AvgLongevity_foliage + self.psi_r / self.AvgLongevity_root +
            self.alpha_w * self.c_H * self.Nconc_foliage / self.AvgLongevity_wood
        )
        return symb_DM_production


if __name__ == "__main__":

    """
    # Example usage:
    solver = DryMassFoliageSolver()

    psi_r = 0.1
    DM_foliage_C, beta1, beta2 = solver.solve_carbon(psi_r)
    DM_foliage_N, beta3, beta4 = solver.solve_nitrogen(psi_r)

    biomass_solver = BiomassProductionSolver(
        Nconc_foliage=0.02, psi_r=psi_r, DM_foliage=DM_foliage_C,
        AvgLongevity_foliage=3, AvgLongevity_wood=10, AvgLongevity_root=5,
        alpha_w=0.3, c_H=0.5
    )

    DM_production = biomass_solver.solve_total_biomass_production()
    print(f"DM_foliage_C: {DM_foliage_C}, DM_foliage_N: {DM_foliage_N}, DM_production: {DM_production}")
    """
    pass
