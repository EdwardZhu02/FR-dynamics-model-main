#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ====================================
# @File name        : Makela08_alloc_mainsolver.py
# @First created    : 2024/5/15 20:49
# @Author           : Yuzhi Zhu
# @E-mail           : edwardmashed@gmail.com
# ====================================
import numpy as np
from numpy.polynomial.polynomial import Polynomial
from scipy.optimize import minimize_scalar


# Function to calculate W_f (1)
def ss_drymass_foliage_solver_carbon(psi_r, Nconc_foliage,
                                     alpha_w, c_H, NrNf_ratio, NwNf_ratio, Resp_Nspecific, CtoDM_frac, Kf,
                                     AvgLongevity_foliage, AvgLongevity_wood, AvgLongevity_root,
                                     Photosyn_lightsat,
                                     ):
    """
    Eqn. S3, steady state W_f calculation
    `N_f`(Nconc_foliage) and `psi_r` provide variable `W_f`(DM_foliage)

    Pre-given parameters: alpha_w, c_H, NrNf_ratio, NwNf_ratio, Resp_Nspecific, CtoDM_frac, Kf, AvgLongevity_foliage, AvgLongevity_wood, AvgLongevity_root

    Photosyn_lightsat: (sigma_fM) light-saturated foliage-specific rate of photosynthesis
    - kg C / kg foliage / yr

    :param Kf:
    :param CtoDM_frac:
    :param Photosyn_lightsat:
    :param Resp_Nspecific:
    :param AvgLongevity_root:
    :param AvgLongevity_foliage:
    :param AvgLongevity_wood:
    :param psi_r:
    :param Nconc_foliage:
    :param alpha_w:
    :param c_H:
    :param NrNf_ratio: n_r
    :param NwNf_ratio: n_f
    :return: Wf under carbon balance


    """
    beta1 = CtoDM_frac * Photosyn_lightsat * Kf / (1 / AvgLongevity_root + Resp_Nspecific * Nconc_foliage * NrNf_ratio)
    # Kf (kg/ha): density of foliage DM that reduce rate to 50% of the Photosyn_lightsat rate

    beta2 = (1 / AvgLongevity_foliage + Nconc_foliage * (
            alpha_w * c_H / AvgLongevity_wood +
            CtoDM_frac * Resp_Nspecific * (1 + NwNf_ratio * alpha_w * c_H * Nconc_foliage)
    )) / (1 / AvgLongevity_root + CtoDM_frac * Resp_Nspecific * Nconc_foliage * NrNf_ratio)

    DM_foliage_C = beta1 / (beta2 + psi_r) - Kf  # Eqn. S3
    return DM_foliage_C, beta1, beta2


# Function to calculate W_f (2)
def ss_drymass_foliage_solver_nitrogen(psi_r, Nconc_foliage,
                                       alpha_w, c_H, NrNf_ratio, NwNf_ratio, Kr,
                                       AvgLongevity_foliage, AvgLongevity_wood, AvgLongevity_root,
                                       NResorbFrac_foliage, NResorbFrac_wood, NResorbFrac_root,
                                       Nup_max_specific
                                       ):
    """
    Eqn. S7,  steady state W_f calculation
    `N_f`(Nconc_foliage) and `psi_r` provide variable `W_f`(DM_foliage)

    Pre-given parameters: alpha_w, c_H, NrNf_ratio, NwNf_ratio, Kr, AvgLongevity_foliage, AvgLongevity_wood, AvgLongevity_root, NResorbFrac_foliage, NResorbFrac_wood, NResorbFrac_root,

    Umax vary between sites, reflecting an N gradient
    Nup_max_specific: (sigma_rM) maximum fine-root specific N uptake rate, depends on availability of N in the soil
    - kg N / kg FR / yr
    - Somewhat represents soil N availability

    :param Kr:
    :param Nup_max_specific:
    :param NResorbFrac_foliage:
    :param NResorbFrac_wood:
    :param NResorbFrac_root:
    :param AvgLongevity_root:
    :param AvgLongevity_wood:
    :param AvgLongevity_foliage:
    :param psi_r:
    :param Nconc_foliage:
    :param alpha_w:
    :param c_H:
    :param NrNf_ratio:
    :param NwNf_ratio:
    :return:
    """
    Nup_max = Nup_max_specific * Kr  # U_max, maximum rate of N uptake
    beta3 = Nup_max / (
            Nconc_foliage * (1 - NResorbFrac_root) * NrNf_ratio / AvgLongevity_root
    )

    beta4 = (
                    (1 - NResorbFrac_foliage) / AvgLongevity_foliage +
                    (1 - NResorbFrac_wood) * NwNf_ratio * alpha_w * c_H * Nconc_foliage / AvgLongevity_wood
            ) / (
                    (1 - NResorbFrac_root) * NrNf_ratio / AvgLongevity_root
            )

    DM_foliage_N = beta3 / (beta4 + psi_r) - Kr / psi_r  # Eqn. S7
    return DM_foliage_N, beta3, beta4


# Function to calculate G using N_f, psi_r, W_f (Nconc_foliage, psi_r, DM_foliage)
def ss_total_biomass_production_solver(Nconc_foliage, psi_r, DM_foliage,
                                       AvgLongevity_foliage, AvgLongevity_wood, AvgLongevity_root, alpha_w, c_H
                                       ):
    # Production rate of total biomass at steady state
    DM_production = DM_foliage * (
            1 / AvgLongevity_foliage + psi_r / AvgLongevity_root + alpha_w * c_H * Nconc_foliage / AvgLongevity_wood
    )  # Eqn. S1
    return DM_production


def ss_psi_r_solver(beta1, beta2, beta3, beta4, Kr, Kf):
    # a, b, c = (5, 7, 8)
    # Nf = random.randint(1, 10)
    a_1 = float((beta1 - beta3 + Kr - Kf * (beta2 + beta4)) / -Kf)
    a_2 = float((beta1 * beta4 - beta2 * beta3 + Kr * (beta2 + beta4) - Kf * beta2 * beta4) / -Kf)
    a_3 = float(Kr * beta2 * beta4 / -Kf)

    p = Polynomial([a_3, a_2, a_1, 0])  # x^3 + a1x^2 + a2x + a3
    psi_r_roots = p.roots()
    psi_r_real_roots = psi_r_roots[np.isreal(psi_r_roots)].real

    if len(psi_r_real_roots) == 0:
        # No real roots
        print("No real roots found for given Nf\n")
        return None
    else:
        # ---
        # Choose the root > 0, but not always the max.
        psi_r = max(psi_r_real_roots)  # Choose the root that you are interested in, e.g., the maximum root
        # ---
        print("psi_r = ", psi_r, "\n", sep=" ")
        return psi_r


def nconc_foliage_optimizer():
    # Objective function to minimize (negative of G to maximize G)
    def maximum_g_objective(Nconc_foliage,
                            var1, var2,
                            calc_method="C"):

        psi_r = ss_psi_r_solver(Nconc_foliage)
        if calc_method == "C":  # carbon balance
            W_f = ss_drymass_foliage_solver_carbon(Nconc_foliage, psi_r)
        else:  # elif calc_method == "N":  # nitrogen balance
            W_f = ss_drymass_foliage_solver_nitrogen(Nconc_foliage, psi_r)
        G = ss_total_biomass_production_solver(Nconc_foliage, psi_r, W_f)
        return -G  # Minimize negative G to maximize G

    # Create a lambda function to fix the other parameters and only pass N_f (Nconc_foliage) to the optimizer
    optimize_g_func = lambda Nconc_foliage: maximum_g_objective(Nconc_foliage, var1, var2)

    # Use scipy.optimize to find the N_f that maximizes G
    optimize_g_result = minimize_scalar(optimize_g_func, bounds=(0, 10), method='bounded')

    # Extract the optimal N_f and corresponding G
    optimal_Nconc_foliage = optimize_g_result.x
    optimal_psi_r = ss_psi_r_solver(optimal_Nconc_foliage)
    optimal_W_f = ss_drymass_foliage_solver_carbon(optimal_Nconc_foliage, optimal_psi_r)
    optimal_G = ss_total_biomass_production_solver(optimal_Nconc_foliage, optimal_psi_r, optimal_W_f)

    print(f'Optimal N_f: {optimal_Nconc_foliage}')
    print(f'Optimal psi_r: {optimal_psi_r}')
    print(f'Optimal W_f: {optimal_W_f}')
    print(f'Optimal G: {optimal_G}')
