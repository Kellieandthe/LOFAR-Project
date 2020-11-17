# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 11:52:33 2020

@author: ppykd1
"""

import numpy as np
from astropy.table import Table, vstack
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt

# Set cosmology according to Garon et al. paper
cosmo = FlatLambdaCDM(H0=68, Om0=0.315)

# Import previously created FITS table with radio/optical data
RadDat = Table.read('Radio source data', format='fits')
# Import WHL15 cluster catalogues - old, updated catalogue and new smaller cluster catalogue
# No longer using the cluster duplicates code as we decided that keeping the cluster catalogues separate would be a better idea.
WHL15 = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/The WHL12 cluster catalog with updated parameters.txt', format='cds')
WHL15New = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/Newly identified rich clusters at high redshifts.txt', format='cds')

# Create column containing z source as either 'Photometric' or 'Spectroscopic'
WHL15zSource = np.full(len(WHL15['zspec']), 'Spectroscopic')
np.putmask(WHL15zSource, WHL15['zspec'] < -0.9, 'Photometric')
# Replace any photometric z values that have a spectroscopic value instead
np.putmask(WHL15['zphot'], WHL15['zspec'] > -0.9, WHL15['zspec'])

# Only keep (and rename) useful columns
WHL15.keep_columns(['Name', 'RAdeg', 'DEdeg', 'zphot', 'RL*500'])
WHL15.rename_column('zphot', 'z')
WHL15New.keep_columns(['Name', 'RAdeg', 'DEdeg', 'zspec', 'RL*500'])
WHL15New.rename_column('zspec', 'z')

# Combine old and new WHL catalogue data (and rename richness column)
WHLDat = vstack([WHL15, WHL15New])
WHLDat.rename_column('RL*500', 'Richness')

# Combine z Source arrays and add to WHL master table
WHL15NewzSource = np.full(len(WHL15New['z']), 'Spectroscopic') # All new WHL data has spectroscoptic redshift
WHLzSource = np.concatenate([WHL15zSource, WHL15NewzSource])
WHLDat.add_column(WHLzSource, name='z Source')

plt.hist(WHLDat['z'], bins=20)

z = np.array([0.2, 0.6])
num = len(WHLDat['z'][(WHLDat['z'] >= z[0]) & (WHLDat['z'] <= z[1])])
area = 11400
sphArea = 360**2/np.pi

d = cosmo.comoving_distance(z)

sphvol = (4 * np.pi * d**3) / 3
voldiff = sphvol[1] - sphvol[0]
volfrac = (area * voldiff) / sphArea
smVol = volfrac/num
clusDist = smVol**(1/3)

