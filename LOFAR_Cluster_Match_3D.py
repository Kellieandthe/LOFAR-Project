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
Rad = Table.read('Data\Radio source data', format='fits')
Clus = Table.read('Data\Full RM Cluster Catalogue', format='fits') # This input can be changed for either WHL15 or RM catalogues

# Spectroscopic z value cuts
# ClusDat = ClusDat[ClusDat['z Source'] == 'Spectroscopic']
# RadDat = RadDat[RadDat['z_best_source'] == 'Spectroscopic']

# Calculate distances to sources in Mpc using redshift values
optDist = cosmo.angular_diameter_distance(Rad['z_best'])
ClusDist = cosmo.angular_diameter_distance(Clus['Cluster z'])

# Put optical and cluster RA and DEC into SkyCoords
opt = SkyCoord(Rad['Optical RA'], Rad['Optical DEC'], distance=optDist,
               unit=('deg', 'deg', 'Mpc'))
clus = SkyCoord(Clus['Cluster RA'],Clus['Cluster DEC'], distance=ClusDist,
                unit=('deg', 'deg', 'Mpc'))

#%% Single match method

# Find 3D matches between opt and clus
idx, d2d, d3d = opt.match_to_catalog_3d(clus)
d2d = d2d.arcsecond

# Calculate Delta z between matched galaxies and cluster centres
delta_z = lf.z_diff_calc(Rad['z_best'], Clus['Cluster z'][idx])
# Calculate 2D distance on sky in Mpc
Dist2d = optDist*np.radians(d2d/3600)
# Constrain matches to within abs(Delta z) < 0.04 and within 15Mpc
cutCond = (abs(delta_z) < 0.04) & (Dist2d <= 15*u.Mpc)
d3dcut = d3d[cutCond]
d2dcut = d2d[cutCond]
idxcut = idx[cutCond]

# Combine radio/optical data with matched cluster data
ClusMatchData = hstack([Rad[cutCond], Clus[idxcut]])
ClusMatchData.add_columns([delta_z[cutCond], d3dcut, d2dcut, np.array(Dist2d[cutCond])],
                          names=['Delta z', '3D Distance', '2D Separation', '2D Distance'])

ClusMatchData['Cluster RA'].unit = 'deg'
ClusMatchData['Cluster DEC'].unit = 'deg'
ClusMatchData['3D Distance'].unit = 'Mpc'
ClusMatchData['2D Separation'].unit = 'arcsec'
ClusMatchData['2D Distance'].unit = 'Mpc'

# Create new FITS file with all match data in
ClusMatchData.write('RM match data (3D)', format = 'fits')

#%% Multiple match method

# Find 3D matches between opt and clus
idxclus, idxopt, d2d, d3d = opt.search_around_3d(clus, 180*u.Mpc)
d2d = d2d.arcsecond

# Calculate Delta z between matched galaxies and cluster centres
delta_z = lf.z_diff_calc(Rad['z_best'][idxopt], Clus['Cluster z'][idxclus])
# Calculate 2D distance on sky in Mpc
Dist2d = optDist[idxopt]*np.radians(d2d/3600)
# Constrain matches to within abs(Delta z) < 0.04 and within 15Mpc
cutCond = (abs(delta_z) < 0.04) & (Dist2d <= 15*u.Mpc)
d3dcut = d3d[cutCond]
d2dcut = d2d[cutCond]
idxcluscut = idxclus[cutCond]
idxoptcut = idxopt[cutCond]

# Combine radio/optical data with matched cluster data
ClusMatchData = hstack([Rad[idxoptcut], Clus[idxcluscut]])
ClusMatchData.add_columns([delta_z[cutCond], d3dcut, d2dcut, np.array(Dist2d[cutCond])],
                          names=['Delta z', '3D Distance', '2D Separation', '2D Distance'])

ClusMatchData['Cluster RA'].unit = 'deg'
ClusMatchData['Cluster DEC'].unit = 'deg'
ClusMatchData['3D Distance'].unit = 'Mpc'
ClusMatchData['2D Separation'].unit = 'arcsec'
ClusMatchData['2D Distance'].unit = 'Mpc'

#%%
match = ClusMatchData.copy()

for i in range(0, len(Rad)):
    if i%500 == 0: # Checking how fast the code is running
        print(i)
    # Look for all sources based on ID (looking for sources with >1 match)
    keepCond = match['Radio Source'] == Rad['Radio Source'][i]
    keepDist = match['2D Distance'][keepCond]
    if len(keepDist) > 1:
        allInd = np.where(keepCond == True)[0]
        keepInd = allInd[np.where(keepDist == min(keepDist))[0][0]]
        dropInd = allInd[allInd != keepInd]
        match.remove_rows(list(dropInd))
        
        
#%%

# Create new FITS file with all match data in
match.write('RM match data (new 3D)', format = 'fits')








