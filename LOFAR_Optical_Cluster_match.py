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

'''
This is currently my preferred method of matching a cluster to a radio-optical source.
However, my supervisors would prefer me to follow the same method used in a paper
I am following, so that is what the module Slow_Cluster_Match is. Please look
at that module for where I am having issues.
'''

# Set cosmology according to Garon et al. paper
cosmo = FlatLambdaCDM(H0=68, Om0=0.315)

# Import previously created FITS tables with radio/optical data and cluster data
RadDat = Table.read('Radio source data', format='fits')
ClusDat = Table.read('Full WHL Cluster Catalogue') # This input can be changed for either WHL15 or RM catalogues

# Calculate distances to sources in Mpc using redshift values
optDist = cosmo.comoving_distance(RadDat['Optical z'])
ClusDist = cosmo.comoving_distance(ClusDat['Cluster z'])

# Put optical and cluster RA and DEC into SkyCoords
opt = SkyCoord(RadDat['Optical RA'], RadDat['Optical DEC'], distance=optDist,
               unit=('deg', 'deg', 'Mpc'))
clus = SkyCoord(ClusDat['Cluster RA'],ClusDat['Cluster DEC'], distance=ClusDist,
                unit=('deg', 'deg', 'Mpc'))

# Find 3D matches between opt and clus
idx, d2d, d3d = opt.match_to_catalog_3d(clus)
d2d = d2d.arcsecond

# Constrain matches to within 15Mpc
cutCond = d3d <= 15*u.Mpc
d3dcut = d3d[cutCond]
d2dcut = d2d[cutCond]
idxcut = idx[cutCond]
# Calculate 2D distance on sky in Mpc
Dist2d = optDist[cutCond]*np.radians(d2dcut/3600)
# Calculate Delta z between matched galaxies and cluster centres
delta_z = lf.z_diff_calc(RadDat['Optical z'][cutCond], ClusDat['Cluster z'][idxcut])

# Combine radio/optical data with matched cluster data
ClusMatchData = hstack([RadDat[cutCond], ClusDat[idxcut]])
ClusMatchData.add_columns([delta_z, d3dcut, d2dcut, Dist2d], names=['Delta z',
                          '3D Distance', '2D Separation', '2D Distance'])

ClusMatchData['Cluster RA'].unit = 'deg'
ClusMatchData['Cluster DEC'].unit = 'deg'
ClusMatchData['3D Distance'].unit = 'Mpc'
ClusMatchData['2D Separation'].unit = 'arcsec'
ClusMatchData['2D Distance'].unit = 'Mpc'

# Create new FITS file with all match data in
# ClusMatchData.write('WHL Cluster match data', format = 'fits')






