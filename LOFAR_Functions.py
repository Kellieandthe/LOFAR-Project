# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 09:32:51 2020

@author: kelli
"""

import math


def x_calc(RA_1, RA_0, Dec_0):
    
    """
    Parameters
    ----------
    RA_1 : float
        RA of object 1
    RA_0 : float
        RA of object 0
    Dec_0 : float
        Dec of object 0

    Returns
    -------
    x_val : float
        x offset of objects 0 and 1

    """
    x_val = (RA_1 - RA_0) * math.cos(math.radians(Dec_0)) * 3600
    return x_val

def y_calc(Dec_1, Dec_0):
    """

    Parameters
    ----------
    Dec_1 : float
        Dec of object 1
    Dec_0 : float
        Dec of object 0

    Returns
    -------
    y_val : float
        y offset of objects 0 and 1

    """
    y_val = (Dec_1 - Dec_0) * 3600
    return y_val

def offset(RA_1, RA_0, Dec_1, Dec_0):
    """

    Parameters
    ----------
    RA_1 : float
        RA of object 1
    RA_0 : float
        RA of object 0
    Dec_1 : float
        Dec of object 1
    Dec_0 : float
        Dec of object 0

    Returns
    -------
    Offset : float
        Pythagoras distance from offsets x and y

    """
    x = (RA_1 - RA_0) * math.cos(math.radians(Dec_0)) * 3600
    y = (Dec_1 - Dec_0) * 3600
    Offset = math.sqrt(x**2 + y**2)
    return Offset

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

c = 3*10**5 #km/s   -   speed of light
H_0 = 68 #(km/s)/Mpc   -   Hubble constant

def distfE_calc(z):
    """

    Parameters
    ----------
    z : float
        Redshift of celestial object

    Returns
    -------
    d_fE : float
        Distance from Earth in Mpc

    """
    d_fE = (c*z)/H_0
    return d_fE

def distGC_calc(z, Offset):
    """
    
    Parameters
    ----------
    z : float
        Redshift of celestial objects
    Offset : float
        Offset between galaxy and cluster in arcseconds

    Returns
    -------
    d_GC : float
        Distance between galaxy and cluster in Mpc

    """
    d_GC = (c*z*(Offset/206265))/H_0
    return d_GC

def vertAngle(x, y):
    """

    Parameters
    ----------
    x : float
    y : float

    Returns
    -------
    theta : float
        Angle in degrees from vertical axis

    """
    theta = math.degrees(math.atan(x/y))
    return theta


