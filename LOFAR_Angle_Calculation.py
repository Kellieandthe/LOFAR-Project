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
dat = Table.read('Cluster match data')

# Remove sources that do not have a cluster match
condition = dat['Cluster RA'].data != b'None'
# I have edited my code in the meantime so this will be applicable when I'm able to run it in full:
# condition = np.isnan(dat['Cluster RA'].data) == False
dat = dat[condition]

# Put all radio, optical and cluster RA and DEC into SkyCoords
rad = SkyCoord(dat['Rad RA'], dat['Rad Dec'], unit='deg')
opt = SkyCoord(dat['Opt RA'], dat['Opt Dec'], unit='deg')
clus = SkyCoord(dat['Cluster RA'], dat['Cluster Dec'], unit='deg')

# Use Astropy function to determine angle between optical and radio, and optical and cluster
theta_or = opt.position_angle(rad).degree
theta_oc = opt.position_angle(clus).degree
# Find difference between angles
theta_diff = abs(theta_or - theta_oc)
# If difference is > 180, take 360 - theta_diff
np.putmask(theta_diff, theta_diff > 180, (360 - theta_diff))

#%%

plt.close('all')

WHLCond = dat['Cluster Catalogue'].data == b'WHL15'
rMCond = dat['Cluster Catalogue'].data == b'redMaPPer'

# Create subplot showing angle distribution for both catalogues together and individually
plt.subplot(3, 1, 1)
plt.hist(theta_diff, bins=36)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15 & redMaPPer')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(3, 1, 2)
plt.hist(theta_diff[WHLCond], bins=36)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(3, 1, 3)
plt.hist(theta_diff[rMCond], bins=36)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between Cluster-Optical-Radio sources')
plt.tight_layout()

dat['Dist_oc'] = dat['Dist_oc'].astype(float)

# Create subplot showing distance of galaxy from cluster centre for both catalogues together and individually
plt.figure()
plt.subplot(3, 1, 1)
plt.hist(dat['Dist_oc'], bins=30)
plt.xlabel('Distance (Mpc)')
plt.ylabel('Number of sources')
plt.title('WHL15 & redMaPPer')
plt.xticks(np.arange(0, 16, 1))

plt.subplot(3, 1, 2)
plt.hist(dat['Dist_oc'][WHLCond], bins=30)
plt.xlabel('Distance (Mpc)')
plt.ylabel('Number of sources')
plt.title('WHL15')
plt.xticks(np.arange(0, 16, 1))

plt.subplot(3, 1, 3)
plt.hist(dat['Dist_oc'][rMCond], bins=30)
plt.xlabel('Distance (Mpc)')
plt.ylabel('Number of sources')
plt.title('redMaPPer')
plt.xticks(np.arange(0, 16, 1))
plt.suptitle('Distance between cluster centre and optical source')
plt.tight_layout()

# Adding to the FITS table to include the angle between ROC
dat.add_column(theta_diff, name='Angle ROC')
dat['Angle ROC'].unit = 'deg'

# dat.write('Angle data', format = 'fits')
    
    
    
    
    
    