#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ====================================
# @File name        : Schym08-watbal-main.py
# @First created    : 2024/5/29 10:18
# @Author           : Yuzhi Zhu
# @E-mail           : edwardmashed@gmail.com
# ====================================
import math
import numpy as np


class HourlyFluxSolver:
    def __init__(self, unsat_layer_count, k_unsat, suchead_matrix, layer_thickness):
        """

        :param unsat_layer_count:
        :param k_unsat:
        :param suchead_matrix:
        :param layer_thickness:
        """
        self.unsat_layer_count = unsat_layer_count  # int: wlayer_
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

    def solve_seepage_face_flow(self, elevation_water_table, elevation_channel_avg, elevation_soil_avg,
                                k_sat_saturated):
        """
        :param k_sat_saturated: s_ksat: Saturated hydraulic conductivity (m/s)
        :param elevation_water_table: zw_: Elevation of water table (m)
        :param elevation_channel_avg: i_zr: Average channel elevation (m)
        :param elevation_soil_avg: i_cz: Average soil elevation(m)
        :return: flow_seepage_face: Q_sf (m/s)
        """
        capital_gamma_s = 0  # i_cgs, parameter: Capital Gamma S (length scale for seepage outflow) (m)

        if elevation_water_table > elevation_channel_avg:
            flow_seepage_face = max(0.0, 0.5 * (
                    math.sqrt(elevation_soil_avg - elevation_channel_avg) -
                    math.sqrt(elevation_soil_avg - elevation_water_table)
            ) * (
                    elevation_water_table - elevation_channel_avg
            ) * k_sat_saturated[self.unsat_layer_count - 1] / (
                    math.sqrt(elevation_soil_avg - elevation_channel_avg) * capital_gamma_s * math.cos(i_go))
                                    )  # max
            # k_sat_saturated[self.unsat_layer_count - 1]: the bottom-most layer
        else:
            flow_seepage_face = None

        return flow_seepage_face


def solve_evaporative_flux(soil_conductivity):
    # E_s: soil evaporation per unit horizontal catchment area (m/s)
    water_molar_weight = 0.018  # kg/mol, M_w
    water_density = 1000  # kg/m3, pho_w

    water_molar_frac_laminar = 0  # W_s
    water_molar_frac_atmosphere = 0  # W_a

    soil_evap_flux = ((water_molar_weight / water_density) * soil_conductivity *
                      (water_molar_frac_laminar - water_molar_frac_atmosphere))
    # esoil__ = (par_h(th_) / (srad2par_h * l_E_ * rho_wat)) * (
    #            1 - (1 - i_trans_vegcov) * (o_cait + caig_d(2))) * su__(1)
