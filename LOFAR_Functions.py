# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 09:32:51 2020

@author: kelli
"""

import numpy as np
import math
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM

# Set cosmology according to Garon et al. paper
cosmo = FlatLambdaCDM(H0=68, Om0=0.315)


def x_calc(RA_1, RA_0, Dec_0):
    
    """
    Parameters
    ----------
    RA_1 : float
        RA of object 1 (deg)
    RA_0 : float
        RA of object 0 (deg)
    Dec_0 : float
        Dec of object 0 (deg)

    Returns
    -------
    x_val : float
        x offset of objects 0 and 1 in arcseconds

    """
    x_val = (RA_1 - RA_0) * np.cos(math.radians(Dec_0)) * 3600
    return x_val

def y_calc(Dec_1, Dec_0):
    """

    Parameters
    ----------
    Dec_1 : float
        Dec of object 1 (deg)
    Dec_0 : float
        Dec of object 0 (deg)

    Returns
    -------
    y_val : float
        y offset of objects 0 and 1 in arcseconds

    """
    y_val = (Dec_1 - Dec_0) * 3600
    return y_val

def sepError(RA0, DEC0, RA1, DEC1, RA0Err, DEC0Err, RA1Err, DEC1Err):
    """

    Parameters
    ----------
    RA0 : float
        RA of object 0 (deg)
    Dec0 : float
        DEC of object 0 ((deg)
    RA1 : float
        RA of object 1 (deg)
    Dec1 : float
        DEC of object 1 (deg)
    RA0Err : float
        Error on RA of object 0 (arcsec)
    Dec0Err : float
        Error on DEC of object 0 (arcsec)
    RA1Err : float
        Error on RA of object 1 (arcsec)
    Dec1Err : float
        Error on DEC of object 1 (arcsec)

    Returns
    -------
    sepErr : float
        Error on the separation between objects 0 and 1 in arcseconds

    """
    x = (RA1 - RA0) * np.cos(np.radians(DEC0)) * 3600
    y = (DEC1 - DEC0) * 3600
    xErr = np.sqrt(np.cos(np.radians(DEC0))**2 *(RA1Err**2 + RA0Err**2) \
                     + np.sin(np.radians(DEC0))**2 *(RA0 - RA1)**2 * DEC0Err**2)
    yErr = np.sqrt(DEC1Err**2 + DEC0Err**2)
    sepErr = np.sqrt((x**2 * xErr**2 + y**2 * yErr**2) / (x**2 + y**2))
    return sepErr

def z_diff_calc(z_gal, z_clus):
    """

    Parameters
    ----------
    z_gal : float
        Redshift of galaxy
    z_clus : float
        Redshift of cluster centre

    Returns
    -------
    z_diff : float
        Delta z of galaxy and cluster centre

    """
    z_diff = (z_gal - z_clus) / (1 + z_gal)
    return z_diff

def offset(RA_0, DEC_0, RA_1, DEC_1):
    """

    Parameters
    ----------
    RA0 : float
        RA of object 0 (deg)
    DEC0 : float
        DEC of object 0 ((deg)
    RA1 : float
        RA of object 1 (deg)
    DEC1 : float
        DEC of object 1 (deg)

    Returns
    -------
    sepErr : float
        Separation between objects 0 and 1 in arcseconds

    """
    x = (RA_1 - RA_0) * np.cos(np.radians(DEC_0)) * 3600
    y = (DEC_1 - DEC_0) * 3600
    Offset = np.sqrt(x**2 + y**2)
    return Offset

def NAT_cut(Table):
    NATcut = Table['NAT'] == True
    return Table[NATcut]

def WAT_cut(Table):
    WATcut = Table['WAT'] == True
    return Table[WATcut]

def ext_cut(Table):
    extcut = (Table['NAT'] == True) | (Table['WAT'] == True)
    return Table[extcut]

def AGN_cut(Table):
    AGNcut = (Table['specAGN'] == 1.0) |\
             (Table['mqcAGN'] == True) |\
             (Table['XrayClass'] == 1.0)
    return Table[AGNcut]

def cut_Mpc(Table, num):
    cutCond = (Table['2D Distance'].data <= num)
    return Table[cutCond]

def cut_delZ(Table, num):
    cutCond = (abs(Table['Delta z'].data) <= num)
    return Table[cutCond]
              
# Define BCG region - exclusive to WHL matches
def BCG_cut(Table):
    BCGcut = Table['2D Distance'].data <= 0.01*Table['r500']
    return Table[BCGcut]

def Lum_calc(S_obs, z):
    # Take flux value in Wm^-2Hz^-1
    # Luminosity distance in Mpc, conversion to metres
    D_L = np.array(cosmo.luminosity_distance(z)*3.0857*10**22/u.Mpc)
    D_L = D_L.astype(float)
    # Spectral index for LOFAR data
    alpha = 0.7
    # Luminosity calculation
    L = (S_obs*4*np.pi*D_L**2)/((1+z)**(1+alpha))
    return L

def Vega_conv(Table):
    W1vega = Table['w1Mag'] - 2.699
    W2vega = Table['w2Mag'] - 3.339
    W3vega = Table['w3Mag'] - 5.174
    W4vega = Table['w4Mag'] - 6.620
    Table.add_columns([W1vega, W2vega, W3vega, W4vega], names=['w1Vega', 'w2Vega', 'w3Vega', 'w4Vega'])










