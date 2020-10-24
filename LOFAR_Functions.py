# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 09:32:51 2020

@author: kelli
"""

import math


def x_calc(RA_r, RA_o, Dec_o):
    x_val = (RA_r - RA_o) * math.cos(math.radians(Dec_o)) * 3600
    return x_val

def y_calc(Dec_r, Dec_o):
    y_val = (Dec_r - Dec_o) * 3600
    return y_val

def pythag(x, y):
    Offset = math.sqrt(x**2 + y**2)
    return Offset

def offset(RA_r, RA_o, Dec_r, Dec_o):
    x = (RA_r - RA_o) * math.cos(math.radians(Dec_o)) * 3600
    y = (Dec_r - Dec_o) * 3600
    Offset = math.sqrt(x**2 + y**2)
    return Offset

def z_diff_calc(z_gal, z_clus):
    z_diff = (z_gal - z_clus) / (1 + z_gal)
    return z_diff

c = 3*10**5 #km/s
H_0 = 68 #(km/s)/Mpc

def distfE_calc(z):
    d_fE = (c*z)/H_0
    return d_fE

def distGC_calc(z, Offset):
    # d_GC = abs((c*z*math.tan(Offset/206265))/H_0)
    d_GC = (c*z*(Offset/206265))/H_0
    return d_GC

def angleCalc(d_oc, d_or, d_rc):
    theta = math.acos((d_oc**2 + d_or**2 - d_rc**2)/(2 * d_oc * d_or))
    return theta

def vertAngle(x, y):
    theta = math.degrees(math.atan(x/y))
    return theta

def horAngle(x, y):
    theta = math.degrees(math.atan(y/x))
    return theta

