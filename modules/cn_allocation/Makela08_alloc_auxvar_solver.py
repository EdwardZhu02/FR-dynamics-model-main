#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ====================================
# @File name        : Makela08_alloc_auxvar_solver.py
# @First created    : 2024/5/15 20:49
# @Author           : Yuzhi Zhu
# @E-mail           : edwardmashed@gmail.com
# ====================================
from modules.cn_allocation import Makela08_alloc_mainsolver


class DownstreamValueSolver(Makela08_alloc_mainsolver.BaseSolver):
    def __init__(self, params_dict, Nup_max_specific, Nconc_foliage_opt, DM_production_opt):
        """
        :param params_dict: inherited from Makela08_alloc_mainsolver.BaseSolver
        :param Nup_max_specific: sigma_rM
        :param Nconc_foliage_opt: [N]f
        :param DM_production_opt: G
        """
        super().__init__(params_dict)
        self.Nconc_foliage = Nconc_foliage_opt  # Overwrite the BaseSolver init values
        self.Nup_max_specific = Nup_max_specific
        self.DM_production = DM_production_opt

        self.params_dict = params_dict  # Init the whole dict for psi_r solving

    def solve_optimum_biomass(self, method="C"):

        # Step1: derive psi_r under the optimum Nconc_foliage
        psi_r_solver = Makela08_alloc_mainsolver.PsiRCubicEqnSolver(self.params_dict)

        Photosyn_lightsat = psi_r_solver.solve_photosyn_rate_lightsat_Ndep(self.Nconc_foliage)
        _, _, psi_r_realroots_conv_positive, _ = psi_r_solver.solve_cubic_eqn_numeric(
            self.Nup_max_specific, Photosyn_lightsat, self.Nconc_foliage)
        psi_r_real = psi_r_realroots_conv_positive[0]

        # Step2: derive Wf under the optimum Nconc_foliage
        Wf_solver = Makela08_alloc_mainsolver.DryMassFoliageSolver(
            self.Nup_max_specific, Photosyn_lightsat, self.params_dict,
            Nconc_foliage=self.Nconc_foliage, use_numeric_Nconc_foliage=True
        )

        if method == "C":
            DM_foliage = Wf_solver.solve_carbon(psi_r_real)[0]
        else:  # method == "N"
            DM_foliage = Wf_solver.solve_nitrogen(psi_r_real)[0]
        # print(DM_foliage)

        # Step3: calculate steady-state Wr and Ww based on optimum steady-state Wf
        DM_root = psi_r_real * DM_foliage
        # TODO: Check DM_wood calculation equation, now not so concordant with the manuscript
        DM_wood = self.alpha_w * DM_foliage * (self.c_H * self.Nconc_foliage)  # Eqn. 12 and 15

        return DM_foliage[0], DM_wood[0], DM_root[0]


