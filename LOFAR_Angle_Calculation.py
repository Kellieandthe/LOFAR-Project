# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:02:52 2020

@author: -
"""

from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord

# Import Cluster match data FITS file created in LOFAR_Optical_Cluster_match.py
RMdat = Table.read('RM Cluster match data')
WHLdat = Table.read('WHL Cluster match data')

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

#%%

plt.close('all')

# Create subplot showing angle distribution for both catalogues
plt.subplot(2, 1, 1)
plt.hist(WHLtheta_diff, bins=36)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 1, 2)
plt.hist(RMtheta_diff, bins=36)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between Cluster-Optical-Radio sources')
plt.tight_layout()

# Create subplot showing distance of galaxy from cluster centre for both catalogues
plt.figure()
plt.subplot(2, 1, 1)
plt.hist(WHLdat['3D Distance'], bins=30)
plt.xlabel('3D Distance (Mpc)')
plt.ylabel('Number of sources')
plt.title('WHL15')
plt.xticks(np.arange(0, 16, 1))

plt.subplot(2, 1, 2)
plt.hist(RMdat['3D Distance'], bins=30)
plt.xlabel('3D Distance (Mpc)')
plt.ylabel('Number of sources')
plt.title('redMaPPer')
plt.xticks(np.arange(0, 16, 1))
plt.suptitle('Distance between cluster centre and optical source')
plt.tight_layout()

plt.figure()
plt.subplot(2, 1, 1)
plt.hist(WHLdat['2D Distance'], bins=30)
plt.xlabel('2D Distance (Mpc)')
plt.ylabel('Number of sources')
plt.title('WHL15')
plt.xticks(np.arange(0, 16, 1))

plt.subplot(2, 1, 2)
plt.hist(RMdat['2D Distance'], bins=30)
plt.xlabel('2D Distance (Mpc)')
plt.ylabel('Number of sources')
plt.title('redMaPPer')
plt.xticks(np.arange(0, 16, 1))
plt.suptitle('Distance between cluster centre and optical source')
plt.tight_layout()

# Create subplot delta z between galaxy and cluster centre for both catalogues
plt.figure()
plt.subplot(2, 1, 1)
plt.hist(WHLdat['Delta z'], bins=15)
plt.xlabel(r'$\Delta$z')
plt.ylabel('Number of sources')
plt.title('WHL15')

plt.subplot(2, 1, 2)
plt.hist(RMdat['Delta z'], bins=15)
plt.xlabel(r'$\Delta$z')
plt.ylabel('Number of sources')
plt.title('redMaPPer')
plt.suptitle(r'$\Delta$z between matched galaxies and cluster centres')
plt.tight_layout()

# Adding to the FITS table to include the angle between ROC
RMdat.add_column(RMtheta_diff, name='Angle ROC')
RMdat['Angle ROC'].unit = 'deg'

WHLdat.add_column(WHLtheta_diff, name='Angle ROC')
WHLdat['Angle ROC'].unit = 'deg'

# RMdat.write('RM Angle data', format = 'fits')
# WHLdat.write('WHL Angle data', format = 'fits')
    
    
    
    
    
    