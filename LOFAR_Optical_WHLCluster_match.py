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
import time

# Currently this code takes forever so that's why I'm measuring how long it takes
start_time = time.time()

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

# Find 3D matches between opt and clus that are within 5Mpc
idx, d2d, d3d = opt.match_to_catalog_3d(clus)
d2d = d2d.arcsecond

cutCond = (d3d >= 0.25*u.Mpc) & (d3d <= 5*u.Mpc)
d3dcut = d3d[cutCond]
d2dcut = d2d[cutCond]
idxcut = idx[cutCond]
Dist2d = optDist[cutCond]*np.radians(d2dcut/3600)
delta_z = lf.z_diff_calc(RadDat['Optical z'][cutCond], WHLDat['Cluster z'][idxcut])

ClusMatchData = hstack([RadDat[cutCond], WHLDat[idxcut]])


ClusMatchData.add_columns([delta_z, d3dcut, d2dcut, Dist2d], names=['Delta z',
                          '3D Distance', '2D Separation', '2D Distance'])

ClusMatchData['Cluster RA'].unit = 'deg'
ClusMatchData['Cluster DEC'].unit = 'deg'
ClusMatchData['3D Distance'].unit = 'Mpc'
ClusMatchData['2D Separation'].unit = 'arcsec'
ClusMatchData['2D Distance'].unit = 'Mpc'

# Create new FITS file with all match data in
ClusMatchData.write('WHL Cluster match data', format = 'fits')

'''
This question is duplicated in both cluster match modules:

I am still in the process of determing which method I should use to match
optical sources to clusters, as this method seems perfectly good, however I am
trying to reproduce results from the Garon et al. paper I am following, and the
method he has used is to find galaxies that are within |Delta z| < 0.04 (see
LOFAR_Functions for definition of Delta z function) of a cluster, and then determine the closest galaxy
within that bound that is within a 15Mpc radius of the cluster centre. This is why
I have kept the below (extremely slow!!) code. I can't seem to find functionality
in Astropy as to how I could make it faster, as the Delta z function I am using
seems quite unique, and any 2D search around the sky functions in Astropy require
a limit in angular separation on the sky, rather than in physical distance (15Mpc).
Currently my thought process is to try and alter the Astropy code for the
search_around_sky function so that it takes z values, uses the Delta z function,
and matches galaxies to a cluster within a boundary (0.04 in my case).

If you have any suggestions they would be greatly appreciated!!
'''

#%%

# Create empty lists to store cluster match data in
delta_z = np.array([])
dist = np.array([])*u.Mpc

MatchID = np.array([])
MatchRA = np.array([])
MatchDec = np.array([])
Matchz = np.array([])
MatchzSource = np.array([])
MatchRich = np.array([])

   
for i in range(0, len(RadDat['Optical z'])): # running through every optical source
    min_d_val = 15*u.Mpc # set a value for the minimum distance value that is out of range of the constraints
    for j in range(0, len(WHLDat['z'])): # running through every WHL cluster
        Dz_val = lf.z_diff_calc(RadDat['Optical z'][i], WHLDat['z'][j]) # calculating delta z from z values of optical source and cluster
        if abs(Dz_val) < 0.04: # constraint on delta z
            d_val = (opt[i].separation(clus[j])).radian * cosmo.comoving_distance(WHLDat['z'][j])
            if d_val < 15*u.Mpc and d_val < min_d_val: # constraint on distance radius and only keeping smallest value
                min_d_val = d_val # set minimum d value to smaller d value
                corr_Dz_val = Dz_val # remember corresponding delta z
                ind = j # remember corresponding index in cluster array
    if min_d_val < 15*u.Mpc:
        dist = np.append(dist, min_d_val) # store minimum cluster distance
        delta_z = np.append(delta_z, corr_Dz_val) # store corresponding delta z
        MatchID = np.append(MatchID, WHLDat['Name'][ind]) # store corresponding cluster ID
        MatchRA = np.append(MatchRA, WHLDat['RAdeg'][ind]) # store corresponding cluster RA
        MatchDec = np.append(MatchDec, WHLDat['DEdeg'][ind]) # store corresponding cluster Dec
        Matchz = np.append(Matchz, WHLDat['z'][ind]) # store corresponding cluster z
        MatchzSource = np.append(MatchzSource, WHLDat['z Source'][ind]) # store corresponding cluster z source
        MatchRich = np.append(MatchRich, WHLDat['Richness'][ind]) # store corresponding cluster richness
    else:
        dist = np.append(dist, np.nan) # if no match is found, store NaN for the optical source
        delta_z = np.append(delta_z, np.nan)
        MatchID = np.append(MatchID, np.nan)
        MatchRA = np.append(MatchRA, np.nan)
        MatchDec = np.append(MatchDec, np.nan)
        Matchz = np.append(Matchz, np.nan)
        MatchzSource = np.append(MatchzSource, np.nan)
        MatchRich = np.append(MatchRich, np.nan)

print("This program took", (time.time() - start_time)/3600, "hours to run")















