# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:02:52 2020

@author: -
"""

from astropy.table import Table, join
# import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord

# Import Cluster match data FITS file created in LOFAR_Optical_Cluster_match.py
dat = Table.read('Data\WHL match data', format='fits')

# Put all radio, optical and cluster RA and DEC into SkyCoords
rad = SkyCoord(dat['Radio RA'], dat['Radio DEC'], unit='deg')
opt = SkyCoord(dat['Optical RA'], dat['Optical DEC'], unit='deg')
clus = SkyCoord(dat['Cluster RA'], dat['Cluster DEC'], unit='deg')

# Use Astropy function to determine angle between optical and radio, and optical and cluster

theta_or = opt.position_angle(rad).degree
theta_oc = opt.position_angle(clus).degree

# Find difference between angles
theta_diff = abs(theta_or - theta_oc)

# If difference is > 180, take 360 - theta_diff
np.putmask(theta_diff, theta_diff > 180, (360 - theta_diff))
theta_diff = 180 - theta_diff

# Adding to the FITS table to include the angle between ROC
dat.add_column(theta_diff, name='Angle ROC')
dat['Angle ROC'].unit = 'deg'

# Import morphological data of sample of radio sources
MorphDat = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/Catalogue of morphological classifications.fits', format='fits')
MorphDat.remove_columns(['RA', 'DEC', 'LM_dec_size', 'LM_Flux'])
MorphDat.rename_column('Source_Name', 'Radio Source')

# Find radio sources in current tables that have morphological data and add those columns onto the end
MorphDat = join(dat, MorphDat, keys='Radio Source', join_type='left')

# Cut to only keep sources within 1 and 100arcsecond radio-optical offset, and z > 0.05
def cut_cond(Table):
    cutCond = (Table['Radio-Optical Offset'].data >= 1) &\
              (Table['Radio-Optical Offset'].data < 100) &\
              (Table['z_best'].data > 0.05)
    return Table[cutCond]

cut_cond(MorphDat.filled(False)).write('WHL angle data', format = 'fits')
    
    
    
    
    