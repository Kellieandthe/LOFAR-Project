# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 12:40:31 2020

@author: ppykd1
"""

from astropy.table import Table, join
import numpy as np
import LOFAR_Functions as lf
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM

# Set cosmology according to Garon et al. paper
cosmo = FlatLambdaCDM(H0=68, Om0=0.315)

# Read cluster catalogues and radio-optical data
Clus = Table.read('Data\Full WHL Cluster Catalogue', format='fits')
Rad = Table.read('Data\Radio source data', format='fits')

# Put RA and DECs into SkyCoord co ordinates
Clusco = SkyCoord(Clus['Cluster RA'], Clus['Cluster DEC'], unit='deg')
Radco = SkyCoord(Rad['Optical RA'], Rad['Optical DEC'], unit='deg')

# Calculate angular diameter distance of all z values in radio-optical data
ADDist = np.array(cosmo.angular_diameter_distance(Rad['z_best'])/u.Mpc)

# Create empty arrays to hold cluster IDs, Delta z values, and dist values
ClusIDs = np.array([])
ClusAllDz = np.array([])
ClusAllDist = np.array([])

# Run through all radio-optical values
for i in range(0, len(Rad)):
    if i%500 == 0: # Checking how fast the code is running
        print(i)
    # Calculate Delta z between current optical source and all clusters
    Dz_vals = lf.z_diff_calc(Rad['z_best'][i], Clus['Cluster z'])
    # Find indices for all Delta z values where absolute value is < 0.04
    ind = np.where(abs(Dz_vals) < 0.04)[0]
    # Check if any Delta z values where abs < 0.04 have been found
    if len(ind) != 0:
        # Find offset between RA and DEC of current optical source and all clusters where ind applies
        offset = np.radians((Radco[i].separation(Clusco[ind])).degree)
        # Calulate distance on sky in Mpc using angular diameter distance to optical source and offset
        d_vals = offset*ADDist[i]
        # Keep the minimum distance value
        d_keep = min(d_vals)
        # Find the corresponding index of cluster with min distance in full cluster catalogue
        ind2 = ind[np.where(d_vals == d_keep)]
        # Check to see if minimum distance values are < 5*r200
        if d_keep < 7*Clus['r500'][ind2]:
            # Append corresponding Cluster IDs, Delta zs, and distances to relevant arrays
            ClusIDs = np.append(ClusIDs, Clus['Cluster ID'][ind2][0])
            ClusAllDz = np.append(ClusAllDz, Dz_vals[ind2][0])
            ClusAllDist = np.append(ClusAllDist, d_keep)
        else: # If no min dist values < 5*r200 have been found, append nan
            ClusIDs = np.append(ClusIDs, np.nan)
            ClusAllDz = np.append(ClusAllDz, np.nan)
            ClusAllDist = np.append(ClusAllDist, np.nan)
    else: # If no Delta z values where abs < 0.04 have been found, append nan
        ClusIDs = np.append(ClusIDs, np.nan)
        ClusAllDz = np.append(ClusAllDz, np.nan)
        ClusAllDist = np.append(ClusAllDist, np.nan)
        
        
#%%
# Append new arrays as columns to Rad
Rad.add_columns([ClusIDs, ClusAllDz, ClusAllDist], names=['Cluster ID', 'Delta z', '2D Distance'])
# Join with Cluster catalogue to get all other relevant information
clusMatch = join(Rad, Clus, keys='Cluster ID', join_type='left')
# Remove radio/optical sources with no cluster match
clusMatch2 = clusMatch[np.isnan(clusMatch['Delta z']) == False]
# Write to FITS file
clusMatch2.write('WHL match data', format = 'fits')

            
                
                