#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ====================================
# @File name        : Guswa_optdepth_mainsolver.py
# @First created    : 2024/6/4 17:27
# @Author           : Yuzhi Zhu
# @E-mail           : edwardmashed@gmail.com
# ====================================
import math
import scipy.special as sp  # For gamma functions
import numpy as np


class OptimalRootingDepthSolver:
    def __init__(self, params_dict, observed_temp):
        """
        Initialize parameters for each calculation
        :param params_dict:
        :param observed_temp:
        """
        self.RootResp_20 = params_dict["RootResp_20"]
        self.RootLengthDensity = params_dict["RootLengthDensity"]
        self.SpecificRootLength = params_dict["SpecificRootLength"]
        self.WUE_photosynthesis = params_dict["WUE_photosynthesis"]
        self.FractionGrowingSeason = params_dict["FractionGrowingSeason"]
        self.RainfallDepth_singevent_mean = params_dict["RainfallDepth_singevent_mean"]
        self.WaterContent_Avail = params_dict["WaterContent_Avail"]

        def convert_root_resp_dynamicQ10(_RootResp_20, _obs_temp):
            """
            Calculate root respiration under observed temperature in each grid using reference respiration @ 20 degrees
            Ref: Tjoelker et al. 2001 - temperature-dependent Q10
            :param _RootResp_20:
            :param _obs_temp:
            :return: root respiration @ observed temperature
            """

            itm_Q10 = 3.22 - 0.046 * _obs_temp
            itm_RootResp_amb = _RootResp_20 * itm_Q10 ** ((_obs_temp - 20) / 10)
            # specific_resp_25 = specific_resp_raw * (resp_q10 ^ ((25 - temp_measure) / 10))

            return itm_RootResp_amb

        self.RootResp_amb = convert_root_resp_dynamicQ10(self.RootResp_20, observed_temp)

    def solve_g08_method(self, Transp_pot, RainfallDepth_tot):
        """

        :param RainfallDepth_tot: P, annual total rainfall depth
            - Ratio P/alpha = number of rainfall events
        :param Transp_pot: PT, potential transpiration rate
            - Calculated using Wilson 2009 model in Yang 2016 manuscript
        :return: list containing:
        1: optimal rooting depth [mm], as per Guswa (2008)
        2: optimal rooting zone storage capacity [mm, i.e. l/m2]
        3: optimal rooting zone storage capacity, normalized by the average depth of a rainfall event (see e.g. Porporato et al. 2004)
        """

        # Eqn. 5, calculate A
        # itm_ stands for "intermediate values"
        itm_A = ((self.RootResp_amb * self.RootLengthDensity) /
                 (self.SpecificRootLength * self.WUE_photosynthesis * Transp_pot * self.FractionGrowingSeason))

        # Eqn. 4, calculate X
        itm_W = RainfallDepth_tot / Transp_pot  # W = mean annual P / PT, unit-less
        if itm_W >= 1:
            itm_X = itm_W * (
                    (1 + (self.WaterContent_Avail / self.RainfallDepth_singevent_mean) * (((1 - itm_W) ** 2) / (2 * itm_A))) -
                    math.sqrt(
                        (self.WaterContent_Avail / self.RainfallDepth_singevent_mean) * (((1 - itm_W) ** 2) / itm_A) +
                        (((self.WaterContent_Avail / self.RainfallDepth_singevent_mean) * (((1 - itm_W) ** 2) / (2 * itm_A))) ** 2)
                    )
            )
        else:
            itm_X = itm_W * (
                    (1 + (self.WaterContent_Avail / self.RainfallDepth_singevent_mean) * (((1 - itm_W) ** 2) / (2 * itm_A))) +
                    math.sqrt(
                        (self.WaterContent_Avail / self.RainfallDepth_singevent_mean) * (((1 - itm_W) ** 2) / itm_A) +
                        (((self.WaterContent_Avail / self.RainfallDepth_singevent_mean) * (((1 - itm_W) ** 2) / (2 * itm_A))) ** 2)
                    )
            )

        # Optimal rooting depth under g08 calculation (Eqn.3, solving for optimal Zr)
        RootingDepth = (self.RainfallDepth_singevent_mean / (self.WaterContent_Avail * (1 - itm_W))) * math.log(itm_X)
        # optimal rooting zone storage capacity
        RootingZoneStorageCap = RootingDepth * self.WaterContent_Avail
        # optimal rooting zone storage capacity, normalized by the average depth of a rainfall event
        RootingZoneStorageCap_norm = RootingDepth / (self.RainfallDepth_singevent_mean / self.WaterContent_Avail)

        return [RootingDepth, RootingZoneStorageCap, RootingZoneStorageCap_norm]

    def solve_g10_method(self, Transp_pot, RainfallDepth_tot,
                         RootingDepth_iter_lower=0, RootingDepth_iter_upper=3001, RootingDepth_iter_interval=0.1):
        """
        Approximation, ref. Speich 2018 manuscript
        Epn.1: Marginal cost = marginal benefit
        Optimal rooting depth achieved when dT/dZe equals (RootResp * RLD) / (SRL * WUE * FracGS)

        Optimal Ze is found by applying Eqn.6 to increasing values of Ze, until the difference to the previous iteration
            is less than or equal to that ratio
        :return:
        """

        def calc_single_transpiration(_Rootingdepth_indv):

            # Calculation of Zn, Eqn.7, speich 2018 manuscript
            RootingDepth_norm_precevents = (self.WaterContent_Avail / self.RainfallDepth_singevent_mean) * _Rootingdepth_indv
            itm_W = RainfallDepth_tot / Transp_pot  # W = mean annual P / PT, unit-less

            Transp_daily_mean = Transp_pot * itm_W - (
                (math.exp(-RootingDepth_norm_precevents) * RootingDepth_norm_precevents ** (itm_W * RootingDepth_norm_precevents - 1)) /
                (sp.gammainc(itm_W * RootingDepth_norm_precevents, RootingDepth_norm_precevents))
            )  # <T>, mm/day, ref. Speich 2018 manuscript function (Eqn.6)

            return Transp_daily_mean

        itm_dTdZe = (self.RootResp_amb * self.RootLengthDensity) / (
                self.SpecificRootLength * self.WUE_photosynthesis * self.FractionGrowingSeason)  # benchmark, Eqn.1

        for RootingDepth_indv in np.arange(RootingDepth_iter_lower, RootingDepth_iter_upper, RootingDepth_iter_interval):
            if RootingDepth_indv == RootingDepth_iter_lower:  # neglect the first calculation
                continue

            Transp_daily_mean_curr = calc_single_transpiration(RootingDepth_indv)
            Transp_daily_mean_last = calc_single_transpiration(RootingDepth_indv - RootingDepth_iter_interval)
            itm_iteration_diff = Transp_daily_mean_curr - Transp_daily_mean_last

            if itm_iteration_diff <= itm_dTdZe:
                return RootingDepth_indv
        return None  # Optimal rooting depth not found

