# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 09:32:48 2020

@author: kelli
"""

from astropy.table import Table
import numpy as np
import LOFAR_Functions as lf

# Import previously created FITS tables with radio/optical data and cluster data
RadDat = Table.read('Radio source data', format='fits')
ClusDat = Table.read('Cluster data', format='fits')

# Extract relevant arrays to be able to match optical sources to clusters
gal_z = RadDat['z'].data
galRA = RadDat['Opt RA'].data
galDec = RadDat['Opt Dec'].data

clusID = ClusDat['Cluster ID'].data
clusRA = ClusDat['RA'].data
clusDec = ClusDat['Dec'].data
clus_z = ClusDat['z'].data
clus_zSource = ClusDat['z Source'].data
clusCat = ClusDat['Catalogue'].data

# Create empty lists to store cluster match data in
delta_z = []
dist = []

MatchID = []
MatchRA = []
MatchDec = []
Matchz = []
MatchzSource = []
MatchCat = []
   
for i in range(0, len(gal_z)): # running through every optical source
    min_d_val = 15 # set a value for the minimum distance value that is out of range of the constraints
    for j in range(0, len(clus_z)): # running through every cluster
        z_val = lf.z_diff_calc(gal_z[i], clus_z[j]) # calculating delta z from z values of optical source and cluster
        if abs(z_val) < 0.04: # constraint on delta z
            offset = lf.offset(clusRA[j], galRA[i], clusDec[j], galDec[i]) # calculating offset between optical source and cluster
            d_val = lf.distGC_calc(clus_z[j], offset) # using offset to calculate distance in Mpc between optical source and cluster
            if d_val < 15 and d_val < min_d_val: # constraint on distance radius and only keeping smallest value
                min_d_val = d_val # set minimum d value to smaller d value
                corr_z_val = z_val # remember corresponding delta z
                ind = j # remember corresponding index in cluster array
    if min_d_val < 15:
        dist.append(min_d_val) # store minimum cluster distance
        delta_z.append(corr_z_val) # store corresponding delta z
        MatchID.append(clusID[ind]) # store corresponding cluster ID
        MatchRA.append(clusRA[ind]) # store corresponding cluster RA
        MatchDec.append(clusDec[ind]) # store corresponding cluster Dec
        Matchz.append(clus_z[ind]) # store corresponding cluster z
        MatchzSource.append(clus_zSource[ind]) # store corresponding cluster z source
        MatchCat.append(clusCat[ind]) # store corresponding cluster catalogue
    else:
        dist.append(np.nan) # if no match is found, store NaN for the optical source
        delta_z.append(np.nan)
        MatchID.append(np.nan)
        MatchRA.append(np.nan)
        MatchDec.append(np.nan)
        Matchz.append(np.nan)
        MatchzSource.append(np.nan)
        MatchCat.append(np.nan)

# Import Radio Source data FITS file created from LOFAR_Radio_Optical_match.py to append matching cluster data to
ClusMatchData = Table.read('Radio source data', format='fits')

ClusMatchData.add_column(MatchID, name='Cluster ID')
ClusMatchData.add_column(MatchRA, name='Cluster RA')
ClusMatchData.add_column(MatchDec, name='Cluster Dec')
ClusMatchData.add_column(Matchz, name='Cluster z')
ClusMatchData.add_column(MatchzSource, name='Cluster z Source')
ClusMatchData.add_column(MatchCat, name='Cluster Catalogue')
ClusMatchData.add_column(delta_z, name='Delta z')
ClusMatchData.add_column(dist, name='Dist_oc')

ClusMatchData['Cluster RA'].unit = 'deg'
ClusMatchData['Cluster Dec'].unit = 'deg'
ClusMatchData['Dist_oc'].unit = 'Mpc'

# Create new FITS file with all match data in
# ClusMatchData.write('Cluster match data', format = 'fits')















