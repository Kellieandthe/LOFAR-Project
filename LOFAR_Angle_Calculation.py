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
RMdat = Table.read('Data\RM match data (Garon)', format='fits')
WHLdat= Table.read('Data\WHL match data (Garon)', format='fits')

# Put all radio, optical and cluster RA and DEC into SkyCoords
RMrad = SkyCoord(RMdat['Radio RA'], RMdat['Radio DEC'], unit='deg')
RMopt = SkyCoord(RMdat['Optical RA'], RMdat['Optical DEC'], unit='deg')
RMclus = SkyCoord(RMdat['Cluster RA'], RMdat['Cluster DEC'], unit='deg')

WHLrad = SkyCoord(WHLdat['Radio RA'], WHLdat['Radio DEC'], unit='deg')
WHLopt = SkyCoord(WHLdat['Optical RA'], WHLdat['Optical DEC'], unit='deg')
WHLclus = SkyCoord(WHLdat['Cluster RA'], WHLdat['Cluster DEC'], unit='deg')

# Use Astropy function to determine angle between optical and radio, and optical and cluster
RMtheta_or = RMopt.position_angle(RMrad).degree
RMtheta_oc = RMopt.position_angle(RMclus).degree
WHLtheta_or = WHLopt.position_angle(WHLrad).degree
WHLtheta_oc = WHLopt.position_angle(WHLclus).degree

# Find difference between angles
RMtheta_diff = abs(RMtheta_or - RMtheta_oc)
WHLtheta_diff = abs(WHLtheta_or - WHLtheta_oc)

# If difference is > 180, take 360 - theta_diff
np.putmask(RMtheta_diff, RMtheta_diff > 180, (360 - RMtheta_diff))
np.putmask(WHLtheta_diff, WHLtheta_diff > 180, (360 - WHLtheta_diff))

# Adding to the FITS table to include the angle between ROC
RMdat.add_column(RMtheta_diff, name='Angle ROC')
RMdat['Angle ROC'].unit = 'deg'
WHLdat.add_column(WHLtheta_diff, name='Angle ROC')
WHLdat['Angle ROC'].unit = 'deg'

# Import morphological data of sample of radio sources
MorphDat = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/Catalogue of morphological classifications.fits', format='fits')
MorphDat.remove_columns(['RA', 'DEC', 'LM_dec_size', 'LM_Flux'])
MorphDat.rename_column('Source_Name', 'Radio Source')

# Find radio sources in current tables that have morphological data and add those columns onto the end
RMMorphDat = join(RMdat, MorphDat, keys='Radio Source', join_type='left')
WHLMorphDat = join(WHLdat, MorphDat, keys='Radio Source', join_type='left')

# Cut to only keep sources within 1 and 100arcsecond radio-optical offset, and z > 0.05
def cut_cond(Table):
    cutCond = (Table['Radio-Optical Offset'].data >= 1) &\
              (Table['Radio-Optical Offset'].data < 100) &\
              (Table['z_best'].data > 0.05)
    return Table[cutCond]

cut_cond(RMMorphDat.filled(False)).write('RM Angle data (Garon)', format = 'fits')
cut_cond(WHLMorphDat.filled(False)).write('WHL Angle data (Garon)', format = 'fits')
    
    
    
    
    
    