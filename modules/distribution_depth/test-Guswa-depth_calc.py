# The Guswa 2008, 2010 water-limited max rooting depth model
# Re-written from Speich et al. 2018 R code.

import pandas as pd
import numpy as np


def netrad(flx_day, elv):
    """
    Auxiliary function for the estimation of net radiation

    Net radiation, required for the estimation of potential evapotranspiration, is provided in the FLUXNET dataset.
    However, there are many data gaps for this variable. This function estimates net radiation (and its components,
    net shortwave and net longwave radiation) from air temperature, global radiation and Vapor Pressure Deficit, for
    which temporal coverage is complete in many cases.

    :param flx_day: dataframe containing daily meteorological measurements, format corresponding to FLUXNET2015 dataset
    :param elv: elevation in m asl
    :return: net radiation
    """
    # Constants
    cp = 1004
    at = 0.95
    sigma = 5.6703E-8
    zeroK = 273.15
    hundredK = 373.15
    pZero = 1013.25
    gammaT = 0.005
    g = 9.81
    rAir = 287.04
    rWat = 461.5
    rhoWat = 999.941
    H = [13.3185, -1.976, -0.6445, -0.1299]
    pp = [-137.0, -75.0, 30.0, 167.0, 236.0, 252.0, 213.0, 69.0, -85.0, -206.0, -256.0, -206.0]

    # Select required variables from FLUXNET data.frame
    rg = flx_day['SW_IN_F']  # Global radiation
    rg_pot = flx_day['SW_IN_POT']  # Potential incoming shortwave radiation
    ta = flx_day['TA_F']  # Air temperature
    vpd = flx_day['VPD_F']  # Vapor pressure deficit

    # Net shortwave radiation
    albedo = 0.1
    rs = (1 - albedo) * rg

    # Auxiliary calculations for net LW: Vapor pressure and relative SSD
    tk = ta + zeroK
    tr = 1 - (hundredK / tk)
    tp = tk + gammaT * elv / 2
    PT = pZero * np.exp(-g * elv / rAir / tp)
    l = 2500800 * (zeroK / tk)
    gamma = PT * cp / l / (rAir / rWat)
    rhoAir = 100 * PT / rAir / tk
    x = 0
    dx = 0
    for iter in range(1, 5):
        x += H[iter - 1] * tr ** iter
        dx += iter * H[iter - 1] * tr ** (iter - 1)
    elx = pZero * np.exp(x)
    del_val = hundredK * elx / tk / tk * dx
    el = np.maximum((elx - vpd), 0.1)

    ssd = np.minimum(1, (((rg / rg_pot) - 0.25) / 0.5))

    # Net longwave radiation - following Schulla (1997)
    rl = -sigma * ((ta + zeroK) ** 4) * (0.52 - (0.065 * np.sqrt(el))) * (0.23 + (0.77 * ssd))  # Schulla 1997

    rn = rs + rl

    return rn


def climate_vars_calculation(flx_day, phen, lai, elv, ksimax, kl):
    """
    Calculates climate variables required by the optimal rooting depth models of Guswa (2008, 2010)
    :param flx_day: dataframe containing daily meteorological measurements, format corresponding to FLUXNET 2015 dataset
    :param phen: dataframe with information on phenology, each row corresponding to one day.
        Must contain a column named "active", with a value greater than 0.0 for days in the growing season, and
        equal to 0.0 for days outside the growing season.
    :param lai: leaf area index during the growing season (numeric)
    :param elv: elevation in m asl (numeric)
    :param ksimax: parameter relating LAI to size of interception reservoir
        (see Menzel (1997) and Vegas-Galdos et al. (2012))
    :param kl: canopy light extinction coefficient [-]

    :return:Returns a list containing the following elements:
#	- pet: mean daily Penman PET for the growing season [mm/day]
#	- prec: mean daily precipitation for the growing season [mm/day]
#	- pfreq: frequency of rainfall events [events/day]
#	- pmdepth: average depth of a rainfall event [mm/event]
#	- ts: mean soil temperature [°C]
#	- mean daily *effective* precipitation (i.e. reaching the ground) for the growing season [mm/day]
#	- fgs: growing season length, expressed as a fraction of the year
    """

    # Set constants
    a = 86400 * 1000  # convert seconds to milliseconds
    rAir = 287.04
    rWat = 461.5
    rhoWat = 999.941
    zeroK = 273.15
    hundredK = 373.15
    gammaT = 0.005
    pZero = 10 ** 13.25
    g = 9.81
    cp = 1004
    H = np.array([13.3185, -1.976, -0.6445, -0.1299])
    pp = np.array([-137.0, -75.0, 30.0, 167.0, 236.0, 252.0, 213.0, 69.0, -85.0, -206.0, -256.0, -206.0])
    k = 0.5

    flx_day['TIMESTAMP'] = pd.to_datetime(flx_day['TIMESTAMP'])
    flx_day[['year', 'month', 'day']] = flx_day['TIMESTAMP'].dt.year, flx_day['TIMESTAMP'].dt.month, flx_day[
        'TIMESTAMP'].dt.day
    flx_day = flx_day.replace(-9999, np.nan)

    # Keep only years of selected period (contained in phen file)
    years = phen['year'].unique()
    flx_day = flx_day[flx_day['year'].isin(years)].reset_index(drop=True)

    # Variables from FLUXNET dataset
    t = flx_day['TA_F']  # Air temperature
    tsoil = flx_day['TS_F_MDS_1']  # Soil temperature
    u = flx_day['WS_F']  # Wind speed
    rn = netrad(flx_day, elv)  # Net radiation (calculated with the auxiliary function netrad)
    # also contained in this file
    vpd = flx_day['VPD_F']  # Vapor pressure deficit

    # Penman potential evaporation
    tk = t + zeroK
    tr = 1 - (hundredK / tk)
    tp = tk + gammaT * elv / 2
    PT = pZero * np.exp(-g * elv / rAir / tp)  # Air pressure
    l = 2500800 * (zeroK / tk)  # Latent heat of vaporization
    gamma = PT * cp / l / (rAir / rWat)  # Psychrometric "constant"
    rhoAir = 100 * PT / rAir / tk  # Air density

    x, dx = 0, 0
    for iter in range(4):
        x += H[iter] * tr ** iter
        dx += iter * H[iter] * tr ** (iter - 1)
    elx = pZero * np.exp(x)  # Atmospheric water vapor pressure at saturation
    # "delta" variable of Penman equation (slope of saturation vapor pressure against temperature)
    penman_delta = hundredK * elx / tk / tk * dx

    # Aerodynamic term of Penman equation, following Penman 1954
    alpha_L = gamma * rhoWat * l / a * 0.263 * (0.5 + 0.537 * u)
    energy = ((penman_delta * rn) + alpha_L * vpd) / (gamma + penman_delta)
    epen = energy * a / rhoWat / l  # Penman evaporation

    # Calculate sub-canopy fluxes to estimate soil evaporation
    asoil = np.maximum(0, rn * np.exp(-k * lai))  # Available energy below the canopy
    pesoil = ((penman_delta * 0.8 * asoil) / (penman_delta + gamma)) * a / rhoWat / l

    # Initialize array with zeros
    pet = np.zeros(len(years))
    tpot = np.zeros(len(years))
    prec = np.zeros(len(years))
    pes = np.zeros(len(years))
    temp = np.zeros(len(years))
    ts = np.zeros(len(years))
    for i in range(len(years)):
        phen_year = phen[phen['year'] == years[i]]
        epen_year = epen[phen['year'] == years[i]]
        pet[i] = np.mean(epen_year[phen_year['active'] > 0 & ~np.isnan(phen_year['active'])])
        p_year = flx_day[flx_day['P_ F']][phen['year'] == years[i]][
            phen_year['active'] > 0 & ~np.isnan(phen_year['active'])]
        prec[i] = np.mean(p_year)
        temp[i] = np.mean(t[phen_year['active'] > 0 & ~np.isnan(phen_year['active'])])
        ts[i] = np.mean(tsoil[phen_year['active'] > 0 & ~np.isnan(phen_year['active'])])

    # Temporal distribution of precipitation
    pmdepth = np.zeros(len(years))
    pfreq = np.zeros(len(years))
    for i in range(len(years)):
        # Get precipitation of the current growing season
        gsprec = flx_day[(phen['year'] == years[i]) & (phen['active'] > 0)]['P_F']
        prevp = False
        nevts = 0
        for j in range(len(gsprec)):
            if gsprec.iloc[j] >= 0.5 and not np.isnan(gsprec.iloc[j]):
                if not prevp:
                    nevts += 1
                prevp = True
            else:
                prevp = False
        pfreq[i] = nevts / len(gsprec)
        pmdepth[i] = np.sum(gsprec) / nevts

    climate_df = pd.DataFrame({'years': years, 'pet': pet, 'prec': prec, 'pfreq': pfreq, 'pmdepth': pmdepth, 'ts': ts})
    clim_mean = climate_df.mean()  # Get mean values over all available years

    # Effective precipitation - see Guswa (2008), Eq. 6
    # Assumption: interception evaporative depth corresponds to canopy interception capacity
    simax = ksimax * np.log10(1 + lai)
    peff = max(clim_mean['prec'] * np.exp(-simax / clim_mean['pmdepth']))
    clim_mean['peff'] = peff

    # Growing season length (fraction of a year)
    fgs = sum(phen['active']) / len(phen)
    clim_mean['fgs'] = fgs

    # Potential transpiration
    tpot = clim_mean['pet'] * (1 - np.exp(-kl * lai)) * 0.75
    clim_mean['tpot'] = tpot

    # Potential soil/understory transpiration
    pes = clim_mean['pet'] * np.exp(-kl * lai)
    clim_mean['pes'] = pes

    clim_mean = clim_mean.iloc[1:]
    return clim_mean


def ord_grass_guswa08(w, tpot, prec, avdepth, whc, rresp, rld, srl, wue, fgs):
    """
    Implementation of the optimal rooting depth function of Guswa (2008)
    :param w: the ratio of effective precipitation `Peff` and potential transpiration `Tpot` (eqn. 3)
    :param tpot: mean daily potential transpiration in the growing season [mm/day]
    :param prec: mean daily precipitation in the growing season [mm/day]
    :param avdepth: Average depth of a precipitation event [mm/day]
    :param whc: Soil water holding capacity [mm/mm]
    :param rresp: Root respiration rate
    :param rld: Root length density
    :param srl: specific rooting length
    :param wue: Water-use efficiency
    :param fgs: fraction of the year occupied by GS (growth season)
    :return: Returns a vector containing 3 values:
        1: optimal rooting zone storage capacity, normalized by the average depth of a rainfall event
            (see e.g. Porporato et al. 2004)
        2: optimal rooting depth [mm], as per Guswa (2008)
        3: optimal rooting zone storage capacity [mm, i.e. l/m2]
    """
    out = np.zeros(3, dtype=float)

    # Eqn. 5, calculate A
    a = (rresp * rld) / (srl * wue * tpot * fgs)

    # Eqn. 4, calculate X
    if w >= 1:
        x = w * ((1 + (whc / avdepth) * (((1 - w) ** 2) / (2 * a))) - np.sqrt(
            (whc / avdepth) * (((1 - w) ** 2) / a) + (((whc / avdepth) * (((1 - w) ** 2) / (2 * a))) ** 2)))
    else:
        x = w * ((1 + (whc / avdepth) * (((1 - w) ** 2) / (2 * a))) + np.sqrt(
            (whc / avdepth) * (((1 - w) ** 2) / a) + (((whc / avdepth) * (((1 - w) ** 2) / (2 * a))) ** 2)))

    # Eqn.3, solving for optimal Ze
    zr = (avdepth / (whc * (1 - w))) * np.log(x)

    out[0] = zr / (avdepth / whc)
    out[1] = zr
    out[2] = zr * whc

    return out


def ord_wood_guswa10(w, tpot, prec, avdepth, whc, fgs, wue, gr, srl, rld):
    """
    New implementation of Guswa's 2010 model.
    Older version (provided as supplement of the HESS discussion paper) gave incorrect results

    Inputs:
    - rresp: Root respiration rate

    :param w: the ratio of effective precipitation `Peff` and potential transpiration `Tpot` (eqn. 3)
    :param tpot: mean daily potential transpiration in the growing season [mm/day]
    :param prec: mean daily precipitation in the growing season [mm/day]
    :param avdepth: Average depth of a precipitation event [mm/day]
    :param whc: Soil water holding capacity [mm/mm]
    :param fgs: fraction of the year occupied by GS (growth season)
    :param wue: Water-use efficiency
    :param gr:
    :param srl: specific rooting length
    :param rld: Root length density
    :return: Returns a vector containing 3 values:
        1: optimal rooting zone storage capacity, normalized by the average depth of a rainfall event
            (see e.g. Porporato et al. 2004)
        2: optimal rooting depth [mm], as per Guswa (2008)
        3: optimal rooting zone storage capacity [mm, i.e. l/m2]
    """

    # As differentiating and rearranging the model of Porporato et al. (2004) (Eq. 6) leads to rather cumbersome
    # expressions, an approximation was used here for the G10 model. It follows from Eq. (1) that the optimal rooting
    # depth is the value of Ze for which dT = dZe equals the ratio Gamma_r * Dr = Lr * wph * fseas.
    # Therefore, the optimal Ze is found by applying Eq. (6) to increasing values of Ze, until the difference to the
    # previous iteration is less than or equal to that ratio.

    # Calculate para
    para = (gr * rld) / (srl * wue * fgs)

    # Create xvals and calculate porp.pred, porp.diff
    xvals = np.arange(0, 3001, 0.1)

    #
    porp_pred = []  # porp_ze(xvals, avdepth, whc, w, tpot)  # assume porp_ze is a function defined elsewhere
    porp_diff = np.diff(porp_pred) * 10

    # Find optimal Ze (rooting depth)
    z = xvals[np.where(porp_diff <= para)[0][0]]

    out = [z / (avdepth / whc), z, z * whc]

    return out
