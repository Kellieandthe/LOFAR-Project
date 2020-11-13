# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 09:32:51 2020

@author: kelli
"""

import numpy as np
import math


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
        x offset of objects 0 and 1

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
        y offset of objects 0 and 1

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
        Error on the separation between objects 0 and 1

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




