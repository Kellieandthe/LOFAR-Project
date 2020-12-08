# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 10:41:22 2020

@author: ppykd1
"""

from astropy.table import Table, vstack
import numpy as np


# Import WHL15 cluster catalogues - old, updated catalogue and new smaller cluster catalogue
# No longer using the cluster duplicates code as we decided that keeping the cluster catalogues separate would be a better idea.
WHL15 = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/The WHL12 cluster catalog with updated parameters.txt', format='cds')
WHL15New = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/Newly identified rich clusters at high redshifts.txt', format='cds')

# Create column containing z source as either 'Photometric' or 'Spectroscopic'
WHL15zBestSource = np.full(len(WHL15['zspec']), 'Spectroscopic')
np.putmask(WHL15zBestSource, WHL15['zspec'] < -0.9, 'Photometric')
# Make z best column
WHL15z_best = WHL15['zphot']
# Replace any photometric z values in z best that have a spectroscopic value instead
np.putmask(WHL15z_best, WHL15['zspec'] > -0.9, WHL15['zspec'])

# Combine old and new WHL catalogue data (and rename richness column)
WHLDat = vstack([WHL15, WHL15New])
WHLDat.rename_columns(['Name', 'RAdeg', 'DEdeg', 'RL*500'], ['Cluster ID',
                        'Cluster RA', 'Cluster DEC', 'Richness'])

# Combine z best arrays
WHLzBest = np.concatenate([WHL15z_best, WHL15New['zspec']]) # All new WHL data has spectroscopic redshift
# Combine z Source arrays and add to WHL master table
WHL15NewzSource = np.full(len(WHL15New['zspec']), 'Spectroscopic') # All new WHL data has spectroscopic redshift
WHLzSource = np.concatenate([WHL15zBestSource, WHL15NewzSource])
WHLDat.add_columns([WHLzBest, WHLzSource], names=['Cluster z', 'z Source'])

# Calculate optical mass proxy from richness and add to WHL master table
M500 = 10**(1.08*np.log10(WHLDat['Richness']) - 1.37)
WHLDat.add_column(M500, name='M500')
# Dud rows that don't feed into SkyCoord nicely
WHLDat.remove_rows([132684, -1])
# Write to its own FITS file
WHLDat.write('Full WHL Cluster Catalogue', format = 'fits')

#%%

# Import redMaPPer cluster catalogue
RMDat = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/Redmapper DR8 Public v6.3 Cluster Catalog.fits.gz', format='fits')

# Create column containing z best source as either 'Photometric' or 'Spectroscopic'
RMzSource = np.full(len(RMDat['Z_SPEC']), 'Spectroscopic')
np.putmask(RMzSource, RMDat['Z_SPEC'] < -0.9, 'Photometric')
# Create z best column
RMzBest = RMDat['Z_LAMBDA']
np.putmask(RMzBest, RMDat['Z_SPEC'] > -0.9, RMDat['Z_SPEC'])

# Add z best and z Source column, and rename columns
RMDat.add_columns([RMzBest, RMzSource], names=['Cluster z', 'z Source'])
RMDat.rename_columns(['NAME', 'RA', 'DEC', 'LAMBDA', 'LAMBDA_ERR', 'S'],
                     ['Cluster ID', 'Cluster RA', 'Cluster DEC', 'Richness',
                      'E_Richness', 'Richness Scale Factor'])

# Write to its own FITS file
RMDat.write('Full RM Cluster Catalogue', format = 'fits')