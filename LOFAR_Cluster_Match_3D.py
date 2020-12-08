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

# Import previously created FITS tables with radio/optical data and cluster data
RadDat = Table.read('Data\Radio source data', format='fits')
ClusDat = Table.read('Data\Full RM Cluster Catalogue', format='fits') # This input can be changed for either WHL15 or RM catalogues

# Calculate distances to sources in Mpc using redshift values
optDist = cosmo.comoving_distance(RadDat['z_best'])
ClusDist = cosmo.comoving_distance(ClusDat['Cluster z'])

# Put optical and cluster RA and DEC into SkyCoords
opt = SkyCoord(RadDat['Optical RA'], RadDat['Optical DEC'], distance=optDist,
               unit=('deg', 'deg', 'Mpc'))
clus = SkyCoord(ClusDat['Cluster RA'],ClusDat['Cluster DEC'], distance=ClusDist,
                unit=('deg', 'deg', 'Mpc'))

# Find 3D matches between opt and clus
idx, d2d, d3d = opt.match_to_catalog_3d(clus)
d2d = d2d.arcsecond

# Calculate Delta z between matched galaxies and cluster centres
delta_z = lf.z_diff_calc(RadDat['z_best'], ClusDat['Cluster z'][idx])
# Calculate 2D distance on sky in Mpc
Dist2d = optDist*np.radians(d2d/3600)
# Constrain matches to within abs(Delta z) < 0.04 and within 15Mpc
cutCond = (abs(delta_z) < 0.04) & (Dist2d <= 15*u.Mpc)
d3dcut = d3d[cutCond]
d2dcut = d2d[cutCond]
idxcut = idx[cutCond]


# Combine radio/optical data with matched cluster data
ClusMatchData = hstack([RadDat[cutCond], ClusDat[idxcut]])
ClusMatchData.add_columns([delta_z[cutCond], d3dcut, d2dcut, np.array(Dist2d[cutCond])],
                          names=['Delta z', '3D Distance', '2D Separation', '2D Distance'])

ClusMatchData['Cluster RA'].unit = 'deg'
ClusMatchData['Cluster DEC'].unit = 'deg'
ClusMatchData['3D Distance'].unit = 'Mpc'
ClusMatchData['2D Separation'].unit = 'arcsec'
ClusMatchData['2D Distance'].unit = 'Mpc'

# Create new FITS file with all match data in
# ClusMatchData.write('WHL Cluster match data', format = 'fits')






