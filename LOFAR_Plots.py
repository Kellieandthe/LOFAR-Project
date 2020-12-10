# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 11:10:17 2020

@author: ppykd1
"""

from astropy.table import Table
import LOFAR_Functions as lf
import matplotlib.pyplot as plt
import numpy as np

RMdat = Table.read('Data\RM Angle data (Garon)', format = 'fits')
WHLdat = Table.read('Data\WHL Angle data (Garon)', format = 'fits')

plt.close('all')

# Create subplot showing distance of galaxy from cluster centre for both catalogues
plt.figure()
plt.hist(WHLdat['2D Distance'], bins=30)
plt.xlabel('2D Distance (Mpc)')
plt.ylabel('Number of sources')
plt.title('Distance between cluster centre and optical source (WHL15)')
plt.xticks(np.arange(0, 16, 1))
plt.tight_layout()

# Create subplot of delta z between galaxy and cluster centre for various 2D distances
plt.figure()
plt.subplot(3, 1, 1)
plt.hist(WHLdat['Delta z'], bins=20)
plt.xlabel(r'$\Delta$z')
plt.ylabel('Number of sources')
plt.title('2D Distance < 15Mpc')

plt.subplot(3, 1, 2)
plt.hist(lf.cut_Mpc(WHLdat, 5)['Delta z'], bins=20)
plt.xlabel(r'$\Delta$z')
plt.ylabel('Number of sources')
plt.title('2D Distance < 5Mpc')

plt.subplot(3, 1, 3)
plt.hist(lf.cut_Mpc(WHLdat, 2)['Delta z'], bins=20)
plt.xlabel(r'$\Delta$z')
plt.ylabel('Number of sources')
plt.title('2D Distance < 2Mpc')
plt.suptitle(r'$\Delta$z between matched galaxies and cluster centres (WHL15)')
plt.tight_layout()

# Create plot of WAT and NAT delta z values
plt.figure()
plt.hist(lf.ext_cut(WHLdat)['Delta z'], bins=20)
plt.xlabel(r'$\Delta$z')
plt.ylabel('Number of sources')
plt.suptitle(r'$\Delta$z between matched WAT and NAT galaxies and cluster centres (WHL15)')
plt.tight_layout()

# Plot ROC angles with new constraints for NATs
plt.figure()
plt.subplot(3, 1, 1)
plt.hist(lf.NAT_cut(WHLdat)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('2D Distance < 15Mpc')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(3, 1, 2)
plt.hist(lf.cut_Mpc(lf.NAT_cut(WHLdat), 5)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('2D Distance < 5Mpc')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(3, 1, 3)
plt.hist(lf.cut_Mpc(lf.NAT_cut(WHLdat), 2)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('2D Distance < 2Mpc')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between cut Radio-Optical-Cluster sources for NATs (WHL15)')
plt.tight_layout()

# Plot ROC angles with new constraints for WATs
plt.figure()
plt.subplot(3, 1, 1)
plt.hist(lf.WAT_cut(WHLdat)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('2D Distance < 15Mpc')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(3, 1, 2)
plt.hist(lf.cut_Mpc(lf.WAT_cut(WHLdat), 5)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('2D Distance < 5Mpc')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(3, 1, 3)
plt.hist(lf.cut_Mpc(lf.WAT_cut(WHLdat), 2)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('2D Distance < 2Mpc')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between cut Radio-Optical-Cluster sources for WATs (WHL15)')
plt.tight_layout()

# Plot ROC angles with new constraints for AGN
plt.figure()
plt.subplot(3, 1, 1)
plt.hist(lf.AGN_cut(WHLdat)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('2D Distance < 15Mpc')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(3, 1, 2)
plt.hist(lf.cut_Mpc(lf.AGN_cut(WHLdat), 5)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('2D Distance < 5Mpc')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(3, 1, 3)
plt.hist(lf.cut_Mpc(lf.AGN_cut(WHLdat), 2)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('2D Distance < 2Mpc')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between cut Radio-Optical-Cluster sources for AGN (WHL15)')
plt.tight_layout()

#Plot richness histograms
plt.figure()
plt.subplot(2, 1, 1)
plt.hist(WHLdat['Richness'], bins=18)
plt.xlabel('Richness')
plt.ylabel('Number of sources')
plt.title('WHL15')

plt.subplot(2, 1, 2)
plt.hist(RMdat['Richness'], bins=18)
plt.xlabel('Richness')
plt.ylabel('Number of sources')
plt.title('redMaPPer')
plt.suptitle('Richness of cut matched galaxy clusters')
plt.tight_layout()