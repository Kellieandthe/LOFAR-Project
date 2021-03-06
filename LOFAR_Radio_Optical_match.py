# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 15:06:53 2020

@author: kelli
"""

from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
import LOFAR_Functions as lf

# Import FITS table of LoTTS radio & optical data
# LOFAR = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/HETDEX associations and optical IDs.fits', format='fits')
LOFAR = Table.read('Data\RLAGN')

# Rename some columns for ease of use
LOFAR.rename_columns(['Source_Name', 'RA', 'E_RA', 'DEC', 'E_DEC', 'ID_name', 
                      'ID_ra', 'ID_dec'],
                     ['Radio Source', 'Radio RA', 'Radio E_RA', 'Radio DEC',
                      'Radio E_DEC', 'Optical Source', 'Optical RA',
                      'Optical DEC'])
#%%
# Remove sources that have no associated redshift value and remove redshift values that are < 0.05
zCond = (np.isnan((LOFAR['z_best']).data) == False) & (LOFAR['z_best'] > 0.05)
LOFAR = LOFAR[zCond]

# Replace values in Optical z Source column with 'Spectroscopic' or 'Photometric'
LOFAR['z_best_source'] = (LOFAR['z_best_source']).astype(bytes) # convert to bytes from float otherwise str won't be accepted
specCond = LOFAR['z_best_source'] == '1.0' # correspond to spectroscopic redshifts
photCond = LOFAR['z_best_source'] == '0.0' # correspond to photometric redshifts
np.putmask(LOFAR['z_best_source'], specCond, 'Spectroscopic')
np.putmask(LOFAR['z_best_source'], photCond, 'Photometric')
LOFAR['z_best_source'] = (LOFAR['z_best_source']).astype(str) # convert to str as all entries are now strings

# Put radio and optical RA and DEC into SkyCoords
rad = SkyCoord(LOFAR['Radio RA'], LOFAR['Radio DEC'], unit='deg')
opt = SkyCoord(LOFAR['Optical RA'], LOFAR['Optical DEC'], unit='deg')
# Find separation between them in arcseconds
offset = (rad.separation(opt)).arcsecond

# Calculate the error in the offset from the error in the radio RA and DEC
offsetErr = lf.sepError(LOFAR['Radio RA'], LOFAR['Radio DEC'],
                        LOFAR['Optical RA'], LOFAR['Optical DEC'],
                        LOFAR['Radio E_RA'], LOFAR['Radio E_DEC'],
                        np.zeros(len(LOFAR['Radio RA'])), # There is no given error for optical RA and DEC so these are given as 0
                        np.zeros(len(LOFAR['Radio RA'])))

# Add offset and errors to table
LOFAR.add_columns([offset*u.arcsec, offsetErr*u.arcsec],
                  names=['Radio-Optical Offset', 'Offset Error'])

LOFAR.write('Radio source data', format='fits')

#%%

plt.close('all')    

# Remove any anomalous offset data
noOutliers = offset < 100

# Plot histogram with logarithmic y axis
plt.hist(offset, bins=50, log=True)
plt.xlabel('Offset (arcseconds)')
plt.ylabel('Number of sources')
plt.title('Offset from radio source to optical source')

# Plot histogram without anomalies
plt.figure()
plt.hist(offset[noOutliers], bins=50, log=True)
plt.xlabel('Offset (arcseconds)')
plt.ylabel('Number of sources')
plt.title('Offset from radio source to optical source (offset > 100arcsec removed)')

# Plot histogram with logarithmic x and y axes
# def plot_loghist(x, bins):
#     hist, bins = np.histogram(x, bins=bins)
#     logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
#     plt.hist(x, bins=logbins, log=True)
#     plt.xscale('log')

# plt.figure()    
# plot_loghist(d_or, 50)
# plt.xlabel('Offset (arcseconds)')
# plt.ylabel('Number of sources')
# plt.title('Offset from radio source to optical source')






    

