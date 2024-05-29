#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ====================================
# @File name        : Schym08-watbal-main.py
# @First created    : 2024/5/29 10:18
# @Author           : Yuzhi Zhu
# @E-mail           : edwardmashed@gmail.com
# ====================================

class HourlyFluxSolver:
    def __init__(self, unsat_layer_count, k_unsat, suchead_matrix, layer_thickness):
        """

        :param unsat_layer_count:
        :param k_unsat:
        :param suchead_matrix:
        :param layer_thickness:
        """
        self.unsat_layer_count = unsat_layer_count  # int
        self.k_unsat = k_unsat  # list
        self.suchead_matrix = suchead_matrix  # list
        self.layer_thickness = layer_thickness  # list, assumed default all = 0.5m

    def solve_water_flux_unsat(self):
        """
        Solve water fluxes in unsaturated layers
        :return:
        """
        upwards_flow_bylayer = []  # Layer1: [0], Layer2: [1]...

        if self.unsat_layer_count > 1:
            # Runoff occurs only from the layer `_water_layer_count`(id),
            # therefore, no downward flow into the layers below is allowed.
            for layer_id in range(0, self.unsat_layer_count):
                # Eqn.2 in SI
                upwards_flow_bylayer[layer_id] = \
                    0.5 * (self.k_unsat[layer_id] + self.k_unsat[layer_id + 1]) * \
                    (
                            (self.suchead_matrix[layer_id] - self.suchead_matrix[layer_id + 1]) /
                            (0.5 * (self.layer_thickness[layer_id] + self.layer_thickness[layer_id + 1])) - 1
                    )

    def solve_evaporative_flux(self, soil_conductivity):
        # E_s: soil evaporation per unit horizontal catchment area (m/s)
        water_molar_weight = 0.018  # kg/mol, M_w
        water_density = 1000  # kg/m3, pho_w

        water_molar_frac_laminar = 0  # W_s
        water_molar_frac_atmosphere = 0  # W_a

        soil_evap_flux = ((water_molar_weight / water_density) * soil_conductivity *
                          (water_molar_frac_laminar - water_molar_frac_atmosphere))
        # esoil__ = (par_h(th_) / (srad2par_h * l_E_ * rho_wat)) * (
        #            1 - (1 - i_trans_vegcov) * (o_cait + caig_d(2))) * su__(1)


