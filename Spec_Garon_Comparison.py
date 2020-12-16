# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 15:51:42 2020

@author: ppykd1
"""

from astropy.table import Table
import LOFAR_Functions as lf
import matplotlib.pyplot as plt
import numpy as np
# from astropy.coordinates import SkyCoord

WAll = Table.read('Data\WHL cut morphological data (Garon)', format = 'fits')
WSpec = Table.read('Data\WHL cut morphological data (Garon Spectroscopic)', format = 'fits')

plt.close('all')

# Plot of ROC angles
plt.subplot(2, 1, 1)
plt.hist(WAll['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('All Garon et al.')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 1, 2)
plt.hist(WSpec['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('Spectroscopic Garon et al.')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between cut Radio-Optical-Cluster sources')
plt.tight_layout()

# Plot of NATs
plt.figure()
plt.subplot(2, 1, 1)
plt.hist(lf.NAT_cut(WAll)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('All Garon et al.')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 1, 2)
plt.hist(lf.NAT_cut(WSpec)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('Spectroscopic Garon et al.')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between cut Radio-Optical-Cluster sources for NATs')
plt.tight_layout()

# Plot of WATs
plt.figure()
plt.subplot(2, 1, 1)
plt.hist(lf.WAT_cut(WAll)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('All Garon et al.')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 1, 2)
plt.hist(lf.WAT_cut(WSpec)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('Spectroscopic Garon et al.')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between cut Radio-Optical-Cluster sources for WATs')
plt.tight_layout()

# Plot of AGNs
plt.figure()
plt.subplot(2, 1, 1)
plt.hist(lf.AGN_cut(WAll)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('All Garon et al.')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 1, 2)
plt.hist(lf.AGN_cut(WSpec)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('Spectroscopic Garon et al.')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between cut Radio-Optical-Cluster sources for AGN')
plt.tight_layout()

# Plot of Delta z
plt.figure()
plt.subplot(2, 1, 1)
plt.hist(WAll['Delta z'], bins=15)
plt.xlabel(r'$\Delta$z')
plt.ylabel('Number of sources')
plt.title('All Garon et al.')

plt.subplot(2, 1, 2)
plt.hist(WSpec['Delta z'], bins=15)
plt.xlabel(r'$\Delta$z')
plt.ylabel('Number of sources')
plt.title('Spectroscopic Garon et al.')
plt.suptitle(r'$\Delta$z between matched galaxies and cluster centres')
plt.tight_layout()








