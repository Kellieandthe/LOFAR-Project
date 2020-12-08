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
RMdatG = Table.read('Data\RM match data (Garon)', format='fits')
WHLdatG = Table.read('Data\WHL match data (Garon)', format='fits')
RMdat3D = Table.read('Data\RM match data (3D)', format='fits')
WHLdat3D = Table.read('Data\WHL match data (3D)', format='fits')

# Put all radio, optical and cluster RA and DEC into SkyCoords
RMradG = SkyCoord(RMdatG['Radio RA'], RMdatG['Radio DEC'], unit='deg')
RMoptG = SkyCoord(RMdatG['Optical RA'], RMdatG['Optical DEC'], unit='deg')
RMclusG = SkyCoord(RMdatG['Cluster RA'], RMdatG['Cluster DEC'], unit='deg')

WHLradG = SkyCoord(WHLdatG['Radio RA'], WHLdatG['Radio DEC'], unit='deg')
WHLoptG = SkyCoord(WHLdatG['Optical RA'], WHLdatG['Optical DEC'], unit='deg')
WHLclusG = SkyCoord(WHLdatG['Cluster RA'], WHLdatG['Cluster DEC'], unit='deg')

RMrad3D = SkyCoord(RMdat3D['Radio RA'], RMdat3D['Radio DEC'], unit='deg')
RMopt3D = SkyCoord(RMdat3D['Optical RA'], RMdat3D['Optical DEC'], unit='deg')
RMclus3D = SkyCoord(RMdat3D['Cluster RA'], RMdat3D['Cluster DEC'], unit='deg')

WHLrad3D = SkyCoord(WHLdat3D['Radio RA'], WHLdat3D['Radio DEC'], unit='deg')
WHLopt3D = SkyCoord(WHLdat3D['Optical RA'], WHLdat3D['Optical DEC'], unit='deg')
WHLclus3D = SkyCoord(WHLdat3D['Cluster RA'], WHLdat3D['Cluster DEC'], unit='deg')

# Use Astropy function to determine angle between optical and radio, and optical and cluster
RMGtheta_or = RMoptG.position_angle(RMradG).degree
RMGtheta_oc = RMoptG.position_angle(RMclusG).degree
WHLGtheta_or = WHLoptG.position_angle(WHLradG).degree
WHLGtheta_oc = WHLoptG.position_angle(WHLclusG).degree

RM3Dtheta_or = RMopt3D.position_angle(RMrad3D).degree
RM3Dtheta_oc = RMopt3D.position_angle(RMclus3D).degree
WHL3Dtheta_or = WHLopt3D.position_angle(WHLrad3D).degree
WHL3Dtheta_oc = WHLopt3D.position_angle(WHLclus3D).degree

# Find difference between angles
RMGtheta_diff = abs(RMGtheta_or - RMGtheta_oc)
WHLGtheta_diff = abs(WHLGtheta_or - WHLGtheta_oc)
RM3Dtheta_diff = abs(RM3Dtheta_or - RM3Dtheta_oc)
WHL3Dtheta_diff = abs(WHL3Dtheta_or - WHL3Dtheta_oc)

# If difference is > 180, take 360 - theta_diff
np.putmask(RMGtheta_diff, RMGtheta_diff > 180, (360 - RMGtheta_diff))
np.putmask(WHLGtheta_diff, WHLGtheta_diff > 180, (360 - WHLGtheta_diff))
np.putmask(RM3Dtheta_diff, RM3Dtheta_diff > 180, (360 - RM3Dtheta_diff))
np.putmask(WHL3Dtheta_diff, WHL3Dtheta_diff > 180, (360 - WHL3Dtheta_diff))

#%%

plt.close('all')

# Create subplot showing angle distribution for both catalogues
plt.subplot(2, 2, 1)
plt.hist(WHL3Dtheta_diff, bins=36)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15 (3D match)')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 2, 2)
plt.hist(RM3Dtheta_diff, bins=36)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer (3D match)')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 2, 3)
plt.hist(WHLGtheta_diff, bins=36)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15 (Garon et al. match)')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 2, 4)
plt.hist(RMGtheta_diff, bins=36)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer (Garon et al. match)')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between Cluster-Optical-Radio sources')
plt.tight_layout()

# Create subplot showing distance of galaxy from cluster centre for both catalogues
plt.figure()
plt.subplot(2, 2, 1)
plt.hist(WHLdat3D['2D Distance'], bins=30)
plt.xlabel('2D Distance (Mpc)')
plt.ylabel('Number of sources')
plt.title('WHL15 (3D match)')
plt.xticks(np.arange(0, 16, 1))

plt.subplot(2, 2, 2)
plt.hist(RMdat3D['2D Distance'], bins=30)
plt.xlabel('2D Distance (Mpc)')
plt.ylabel('Number of sources')
plt.title('redMaPPer (3D match)')
plt.xticks(np.arange(0, 16, 1))

plt.subplot(2, 2, 3)
plt.hist(WHLdatG['2D Distance'], bins=30)
plt.xlabel('2D Distance (Mpc)')
plt.ylabel('Number of sources')
plt.title('WHL15 (Garon et al. match)')
plt.xticks(np.arange(0, 16, 1))

plt.subplot(2, 2, 4)
plt.hist(RMdatG['2D Distance'], bins=30)
plt.xlabel('2D Distance (Mpc)')
plt.ylabel('Number of sources')
plt.title('redMaPPer (Garon et al. match)')
plt.xticks(np.arange(0, 16, 1))
plt.suptitle('Distance between cluster centre and optical source')
plt.tight_layout()

# Create subplot delta z between galaxy and cluster centre for both catalogues
plt.figure()
plt.subplot(2, 2, 1)
plt.hist(WHLdat3D['Delta z'], bins=15)
plt.xlabel(r'$\Delta$z')
plt.ylabel('Number of sources')
plt.title('WHL15 (3D match)')

plt.subplot(2, 2, 2)
plt.hist(RMdat3D['Delta z'], bins=15)
plt.xlabel(r'$\Delta$z')
plt.ylabel('Number of sources')
plt.title('redMaPPer (3D match)')

plt.subplot(2, 2, 3)
plt.hist(WHLdatG['Delta z'], bins=15)
plt.xlabel(r'$\Delta$z')
plt.ylabel('Number of sources')
plt.title('WHL15 (Garon et al. match)')

plt.subplot(2, 2, 4)
plt.hist(RMdatG['Delta z'], bins=15)
plt.xlabel(r'$\Delta$z')
plt.ylabel('Number of sources')
plt.title('redMaPPer (Garon et al. match)')
plt.suptitle(r'$\Delta$z between matched galaxies and cluster centres')
plt.tight_layout()

# Adding to the FITS table to include the angle between ROC
RMdatG.add_column(RMGtheta_diff, name='Angle ROC')
RMdatG['Angle ROC'].unit = 'deg'
RMdat3D.add_column(RM3Dtheta_diff, name='Angle ROC')
RMdat3D['Angle ROC'].unit = 'deg'

WHLdatG.add_column(WHLGtheta_diff, name='Angle ROC')
WHLdatG['Angle ROC'].unit = 'deg'
WHLdat3D.add_column(WHL3Dtheta_diff, name='Angle ROC')
WHLdat3D['Angle ROC'].unit = 'deg'

#RMdatG.write('RM Angle data (Garon)', format = 'fits')
#WHLdatG.write('WHL Angle data (Garon)', format = 'fits')
#RMdat3D.write('RM Angle data (3D)', format = 'fits')
#WHLdat3D.write('WHL Angle data (3D)', format = 'fits')
    
    
    
    
    
    