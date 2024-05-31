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


class WaterBalanceInitializer:

    def __init__(self, elevation_channel_avg, elevation_soil_avg, soil_layer_thickness):
        """
        Initialize conditions of water balance simulation
        :param elevation_water_table: zw_: Elevation of water table (m)
        :param elevation_channel_avg: i_zr: Average channel elevation (m)
        :param elevation_soil_avg: i_cz: Average soil elevation(m)
        :param soil_layer_thickness: s_delz: Thickness of each soil layer (m) - SIMPLIFIED as CONSTANT

        """
        # Subsequent calculation use zw_ only
        self.elevation_water_table = elevation_soil_avg  # zw_ = i_cz
        self.elevation_channel_avg = elevation_channel_avg
        self.soil_layer_thickness = soil_layer_thickness  # assumed default all = 0.5m

        # Calculate unsaturated soil layers
        self.unsat_layer_count = 0

        while self.elevation_water_table > self.elevation_channel_avg:
            self.unsat_layer_count += 1
            self.elevation_water_table -= self.soil_layer_thickness

    def solve_equilibrium_pressure_head(self, soil_layers_total):
        # phead_matric_bylayer = []  # pcap_: Matric pressure head in each layer (at current time step)
        # initial values
        phead_matric_bylayer = np.zeros(int(soil_layers_total))  # TODO: confirm initial value

        for layer_id in range(0, self.unsat_layer_count):
            phead_matric_bylayer[layer_id] = 0.5 * (
                    self.soil_layer_thickness[layer_id + 1] + self.soil_layer_thickness[layer_id]
            ) + phead_matric_bylayer[layer_id + 1]

    """
    pcap_ = np.zeros(s_maxlayer)
    pcap_[wlayer_:] = 0.0

    # Equilibrium pressure head
    for jj in range(wlayer_, 0, -1):
        pcap_[jj - 1] = 0.5 * s_delz[jj] + 0.5 * s_delz[jj - 1] + pcap_[jj]

    sueq = np.zeros(s_maxlayer)

    for jj in range(1, s_maxlayer):
        sueq[jj - 1] = (1.0 / ((0.5 * (s_delz[jj] + s_delz[jj - 1]) * s_avg[jj - 1]) ** s_nvg[jj - 1] + 1.0)) ** c_mvg[
            jj - 1]

    sueq[s_maxlayer - 1] = (1.0 / (
                (0.5 * s_delz[s_maxlayer - 1] * s_avg[s_maxlayer - 1]) ** s_nvg[s_maxlayer - 1] + 1.0)) ** c_mvg[
                               s_maxlayer - 1]

    su__ = (1.0 / ((pcap_ * s_avg) ** s_nvg + 1.0)) ** c_mvg

    # Unsaturated hydraulic conductivity as a function of su
    kunsat_ = ((-su__ ** (1.0 / c_mvg) + 1.0) ** c_mvg - 1.0) ** 2.0 * s_ksat * np.sqrt(su__)
    cH2Ol_s = (-su__ * s_thetar + su__ * s_thetas + s_thetar) * s_delz

    zwnew = zw_
    wlayernew = wlayer_
    pcapnew = pcap_
    sunew = su__
    kunsatnew = kunsat_  # Unsaturated hydraulic conductivity (m s−1)

    return {
        'zwnew': zwnew,
        'wlayernew': wlayernew,
        'pcapnew': pcapnew,
        'sunew': sunew,
        'kunsatnew': kunsatnew,
        'cH2Ol_s': cH2Ol_s
    }
    """


class HourlyFluxSolver:
    def __init__(self, unsat_layer_count, k_unsat, suchead_matrix, layer_thickness):
        """

        :param unsat_layer_count: wlayer_
        :param k_unsat:
        :param suchead_matrix:
        :param layer_thickness: s_delz: Thickness of each soil layer (m)
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
                                            math.sqrt(
                                                elevation_soil_avg - elevation_channel_avg) * capital_gamma_s * math.cos(
                                        i_go))
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
