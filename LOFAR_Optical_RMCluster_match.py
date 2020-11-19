# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 09:32:48 2020

@author: kelli
"""

from astropy.table import Table, hstack
import numpy as np
import LOFAR_Functions as lf
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
# import matplotlib.pyplot as plt


# Set cosmology according to Garon et al. paper
cosmo = FlatLambdaCDM(H0=68, Om0=0.315)

# Import previously created FITS table with radio/optical data and redMaPPer cluster catalogue data
RadDat = Table.read('Radio source data', format='fits')
# No longer using the cluster duplicates code as we decided that keeping the cluster catalogues separate would be a better idea.
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

# Calculate distances to sources in Mpc using redshift values
optDist = cosmo.comoving_distance(RadDat['Optical z'])
ClusDist = cosmo.comoving_distance(RMDat['Cluster z'])

# Put optical and cluster RA and DEC into SkyCoords
opt = SkyCoord(RadDat['Optical RA'], RadDat['Optical DEC'], distance=optDist,
               unit=('deg', 'deg', 'Mpc'))
clus = SkyCoord(RMDat['Cluster RA'],RMDat['Cluster DEC'], distance=ClusDist,
                unit=('deg', 'deg', 'Mpc'))

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
delta_z = lf.z_diff_calc(RadDat['Optical z'][cutCond], RMDat['Cluster z'][idxcut])

# Combine radio/optical data with matched cluster data
ClusMatchData = hstack([RadDat[cutCond], RMDat[idxcut]])
ClusMatchData.add_columns([delta_z, d3dcut, d2dcut, Dist2d], names=['Delta z',
                          '3D Distance', '2D Separation', '2D Distance'])

ClusMatchData['Cluster RA'].unit = 'deg'
ClusMatchData['Cluster DEC'].unit = 'deg'
ClusMatchData['3D Distance'].unit = 'Mpc'
ClusMatchData['2D Separation'].unit = 'arcsec'
ClusMatchData['2D Distance'].unit = 'Mpc'

# Create new FITS file with all match data in
# ClusMatchData.write('RM Cluster match data', format = 'fits')






