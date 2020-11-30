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
WHL15zSource = np.full(len(WHL15['zspec']), 'Spectroscopic')
np.putmask(WHL15zSource, WHL15['zspec'] < -0.9, 'Photometric')
# Replace any photometric z values that have a spectroscopic value instead
np.putmask(WHL15['zphot'], WHL15['zspec'] > -0.9, WHL15['zspec'])

# Only keep (and rename) useful columns
WHL15.keep_columns(['Name', 'RAdeg', 'DEdeg', 'zphot', 'r500', 'RL*500', 'N500'])
WHL15.rename_column('zphot', 'z')
WHL15New.keep_columns(['Name', 'RAdeg', 'DEdeg', 'zspec', 'r500', 'RL*500', 'N500'])
WHL15New.rename_column('zspec', 'z')

# Combine old and new WHL catalogue data (and rename richness column)
WHLDat = vstack([WHL15, WHL15New])
WHLDat.rename_columns(['Name', 'RAdeg', 'DEdeg', 'z', 'RL*500'], ['Cluster ID',
                        'Cluster RA', 'Cluster DEC', 'Cluster z', 'Richness'])

# Combine z Source arrays and add to WHL master table
WHL15NewzSource = np.full(len(WHL15New['z']), 'Spectroscopic') # All new WHL data has spectroscoptic redshift
WHLzSource = np.concatenate([WHL15zSource, WHL15NewzSource])
WHLDat.add_column(WHLzSource, name='z Source')

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

# Create column containing z source as either 'Photometric' or 'Spectroscopic'
RMzSource = np.full(len(RMDat['Z_SPEC']), 'Spectroscopic')
np.putmask(RMzSource, RMDat['Z_SPEC'] < -0.9, 'Photometric')
# Replace any photometric z values that have a spectroscopic value instead
np.putmask(RMDat['Z_LAMBDA'], RMDat['Z_SPEC'] > -0.9, RMDat['Z_SPEC'])

# Add z Source column, and only keep (and rename) useful columns
RMDat.add_column(RMzSource, name='z Source')
RMDat.keep_columns(['NAME', 'RA', 'DEC', 'Z_LAMBDA', 'Z_LAMBDA_ERR', 'LAMBDA',
                    'LAMBDA_ERR', 'S', 'z Source'])
RMDat.rename_columns(['NAME', 'RA', 'DEC', 'Z_LAMBDA', 'Z_LAMBDA_ERR', 'LAMBDA',
                    'LAMBDA_ERR', 'S'], ['Cluster ID', 'Cluster RA', 'Cluster DEC',
                    'Cluster z', 'Cluster E_z', 'Richness', 'E_Richness',
                    'Richness Scale Factor'])

# Write to its own FITS file
RMDat.write('Full RM Cluster Catalogue', format = 'fits')