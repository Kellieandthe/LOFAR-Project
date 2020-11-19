# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 09:32:48 2020

@author: kelli
"""

from astropy.table import Table, vstack, hstack
import numpy as np
import LOFAR_Functions as lf
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
# import matplotlib.pyplot as plt


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
WHLDat.rename_columns(['Name', 'RAdeg', 'DEdeg', 'z', 'RL*500'], ['Cluster ID',
                        'Cluster RA', 'Cluster DEC', 'Cluster z', 'Richness'])

# Combine z Source arrays and add to WHL master table
WHL15NewzSource = np.full(len(WHL15New['z']), 'Spectroscopic') # All new WHL data has spectroscoptic redshift
WHLzSource = np.concatenate([WHL15zSource, WHL15NewzSource])
WHLDat.add_column(WHLzSource, name='z Source')

# Calculate distances to sources in Mpc using redshift values
optDist = cosmo.comoving_distance(RadDat['Optical z'])
ClusDist = cosmo.comoving_distance(WHLDat['Cluster z'])

# Put optical and cluster RA and DEC into SkyCoords
opt = SkyCoord(RadDat['Optical RA'], RadDat['Optical DEC'], distance=optDist, unit=('deg', 'deg', 'Mpc'))
clus = SkyCoord(WHLDat['Cluster RA'],WHLDat['Cluster DEC'], distance=ClusDist, unit=('deg', 'deg', 'Mpc'))

# Find 3D matches between opt and clus
idx, d2d, d3d = opt.match_to_catalog_3d(clus)
d2d = d2d.arcsecond

# Constrain matches to within 0.25Mpc and 15Mpc
cutCond = (d3d >= 0.25*u.Mpc) & (d3d <= 15*u.Mpc)
d3dcut = d3d[cutCond]
d2dcut = d2d[cutCond]
idxcut = idx[cutCond]
# Calculate 2D distance on sky in Mpc
Dist2d = optDist[cutCond]*np.radians(d2dcut/3600)
# Calculate Delta z between matched galaxies and cluster centres
delta_z = lf.z_diff_calc(RadDat['Optical z'][cutCond], WHLDat['Cluster z'][idxcut])

# Combine radio/optical data with matched cluster data
ClusMatchData = hstack([RadDat[cutCond], WHLDat[idxcut]])
ClusMatchData.add_columns([delta_z, d3dcut, d2dcut, Dist2d], names=['Delta z',
                          '3D Distance', '2D Separation', '2D Distance'])

ClusMatchData['Cluster RA'].unit = 'deg'
ClusMatchData['Cluster DEC'].unit = 'deg'
ClusMatchData['3D Distance'].unit = 'Mpc'
ClusMatchData['2D Separation'].unit = 'arcsec'
ClusMatchData['2D Distance'].unit = 'Mpc'

# Create new FITS file with all match data in
# ClusMatchData.write('WHL Cluster match data', format = 'fits')




