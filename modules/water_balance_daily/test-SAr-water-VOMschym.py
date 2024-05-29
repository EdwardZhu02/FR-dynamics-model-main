# The Schymanski 2008, 2009 VOM-ROOT model
# Re-written from Fortran code from GitHub (https://github.com/schymans/VOM).
import random
import numpy as np
import pandas as pd

def waterbalance_fluxes():

    # Use variables declared in `vom_vegwat_mod`
    wlayer_ = random.randint(1,10) # Ramdom value to avoid warnings
    th_ = 0

    # Define variables
    dummy = 0.0
    ii, jj = 0, 0
    inf__ = 0.0  # REAL*8  :: inf__ ! Infiltration rate (m/s)
    infx__ = 0.0  # REAL*8  :: infx__ ! Infiltration excess runoff rate (m/s)
    qbl = np.zeros((wlayer_,))  # INTEGER :: wlayer_ ! Number of soil layers in unsaturated zone
    esoil__ = 0  # REAL*8  :: esoil__ ! Bare soil evaporation rate (m/s)
    spgfcf__ = 0.0  # REAL*8  :: spgfcf__ ! Seepage face flow rate (m/s)

    # Infiltration
    # REAL*8, ALLOCATABLE :: rain_h(:)  ! Hourly rainfall rate (m/s)
    # INTEGER :: th_  ! Hour since start of run
    def rain_h(hour_since_start):
        return random.randint(1, hour_since_start)

    if rain_h(th_) > 0:  # if rainfall rate in a certain hour >0
        if wlayer_ >= 1:  # if there are at least 1 unsaturated soil layers
            # calculate infiltration rate m/s
            inf__ = min(((s_ksat(1) + kunsat_(1)) / 2) * (1 + (2 * pcap_(1)) / s_delz(1)), rain_h(th_))
            # calculate infiltration excess runoff rate (m/s)
            # rainfall - infiltration
            infx__ = rain_h(th_) - inf__
        else:
            # all saturated layers
            inf__, infx__ = 0.0, 0.0
    else:
        # no rainfall
        inf__, infx__ = 0.0, 0.0

    # Unsaturated flow
    # Eqn. 2 in Schymanski 2015 SI?
    if wlayer_ > 1:
        for jj in range(1, wlayer_ - 1):
            qbl[jj] = -0.5 * (2 * (pcap_(jj + 1) - pcap_(jj))) / ((s_delz(jj + 1) + s_delz(jj)) + 1) * (
                        kunsat_(jj + 1) + kunsat_(jj))

    # Soil evaporation
    esoil__ = (par_h(th_) / (srad2par_h * l_E_ * rho_wat)) * (1 - (1 - i_trans_vegcov) * (o_cait + caig_d(2))) * su__(1)

    # Seepage face flow
    if zw_ > i_zr:
        spgfcf__ = max(0, 0.5 * (np.sqrt(i_cz - i_zr) - np.sqrt(i_cz - zw_)) * (zw_ - i_zr) * s_ksat(wlayer_) / (
                    np.sqrt(i_cz - i_zr) * i_cgs * np.cos(i_go)))

    # Ensure no sub-layer overflows
    if max(su__(1: wlayer_)) >= 1:
    if wlayer_ > 1:
        if su__(1) >= 0.99:
            dummy = esoil__ - inf__ + ruptkt__(1) + ruptkg__(1)
            qbl[0] = min(qbl[0], dummy - 1e-16)

    for ii in range(2, wlayer_ - 1):
        if su__(ii) >= 0.99:
            dummy = qbl[ii - 1] + ruptkt__(ii) + ruptkg__(ii)
            qbl[ii] = min(qbl[ii], dummy - 1e-16)

    if wlayer_ > 1:
        if su__(wlayer_) >= 1:
            dummy = -qbl[wlayer_ - 1] - ruptkt__(wlayer_) - ruptkg__(wlayer_)
            spgfcf__ = max(spgfcf__, dummy + 1e-16)
    else:
        if su__(wlayer_) >= 1:
            dummy = inf__ - esoil__ - ruptkt__(1) - ruptkg__(1)
            spgfcf__ = max(spgfcf__, dummy + 1e-16)


    return 0


def water_balance():
    """

    :return:
    """

    # CURRENT SOIL WATER CONTENT
    su_ = 0         # REAL*8, ALLOCATABLE :: su__(:)    ! Soil saturation degree in each layer (-)
    s_thetar = 0    # REAL*8, ALLOCATABLE :: s_thetar(:)  ! Residual soil water content (-)
    s_thetas = 0    # REAL*8, ALLOCATABLE :: s_thetas(:)  ! Saturated soil water content (-)
    s_delz = 0  # REAL*8, ALLOCATABLE :: s_delz(:)  ! Thickness of each soil layer (m)

    cH2Ol_s = (-su_ * s_thetar + su_ * s_thetas + s_thetar) * s_delz
    wc = np.sum(cH2Ol_s)

    # FLUXES (inf, infx, qbl, esoil__, spgfcf__)
    waterbalance_fluxes()

    # CHANGES IN WATER STORAGE (WATERBALANCE)
    if i_no_veg == 0:
        io__ = inf__ - esoil__ - spgfcf__ - np.sum(ruptkt_) - np.sum(ruptkg_)
    else:
        io__ = inf__ - esoil__ - spgfcf__

    iovec = np.zeros(wlayer_)
    iovec[0] = qbl[0] + inf__ - esoil__ - ruptkt_[0] - ruptkg_[0]
    if wlayer_ == 1:
        iovec[0] -= spgfcf__
    for jj in range(1, wlayer_ - 1):
        iovec[jj] = qbl[jj] - qbl[jj - 1] - ruptkt_[jj] - ruptkg_[jj]
    if wlayer_ > 1:
        iovec[wlayer_] = qbl[wlayer_] - qbl[wlayer_ - 1] - ruptkt_[wlayer_] - ruptkg_[wlayer_] - spgfcf__

    # CHANGE IN SATURATION DEGREE
    dsu = -iovec / ((s_thetar - s_thetas) * s_delz)

    # CALCULATION OF MAXIMAL TIME STEP SIZE
    waterbalance_timestep()

    # CALCULATING STATE VARIABLES AT NEXT TIME STEP
    # (sunew, wlayernew, pcapnew, kunsatnew, zwnew)
    waterbalance_update_state()

    waterbalance_diag()