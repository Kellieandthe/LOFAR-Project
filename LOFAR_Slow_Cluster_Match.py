# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 12:40:31 2020

@author: ppykd1
"""

from astropy.table import Table
import numpy as np
import LOFAR_Functions as lf
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM

# Set cosmology according to Garon et al. paper
cosmo = FlatLambdaCDM(H0=68, Om0=0.315)

# Read cluster catalogues and radio-optical data
Clus = Table.read('Full WHL Cluster Catalogue') # Input can be changed depending on cluster catalogue
Rad = Table.read('Radio source data')

# Put RA and DECs into SkyCoord co ordinates
Clusco = SkyCoord(Clus['Cluster RA'], Clus['Cluster DEC'], unit='deg')
Radco = SkyCoord(Rad['Optical RA'], Rad['Optical DEC'], unit='deg')

# Calculate co-moving distance of all z values in radio-optical data
comDist = np.array(cosmo.comoving_distance(Rad['Optical z'])/u.Mpc)

'''
What I would like to do in the below code is:
    - Find Delta z = (z_gal - z_clus) / (1 + z_gal) for every cluster for each
    radio-optical source
    - When |Delta z| < 0.04, Find the physical distance on the sky in Mpc, at
    the z of the radio-optical source
    - Find the radio-optical/cluster match with the minimum physical distance
    on the sky that is < 15Mpc and within |Delta z| < 0.04, for each
    radio-optical source
    
I would prefer not to do this with a nested for loop but I'm not sure how else
to do it! I have a faster alternative for cluster matching which is in the
Optical_Cluster_match module, but my supervisors would prefer for me to use
this method as it follows the same method used in a paper I am looking at.

Please let me know if you know how to make this faster!
'''
#%%

# Create empty arrays to hold cluster IDs, Delta z values, and dist values
ClusIDs = np.array([])
ClusAllDz = np.array([])
ClusAllDist = np.array([])

# Run through all radio-optical values
for i in range(0, len(Rad['Optical RA'])):
    # Set dist values that are outside of the constraint of being < 15Mpc
    Clusdist = 15
    # Run through length of all WHL values
    for j in range(0, len(Clus['Cluster RA'])):
        # Calculate Delta z values between current optical source and cluster
        Clusdz = lf.z_diff_calc(Rad['Optical z'][i], Clus['Cluster z'][j])
        # See if absolute value of Delta z is less than 0.04
        if abs(Clusdz) < 0.04:
            # If so, calculate distance on sky between cluster and optical source
            # Clusoffset = np.radians((Radco[i].separation(Clusco[j])).degree) # Find offset between RA and DEC
            Clusoffset = np.radians(lf.offset(Clus['Cluster RA'][j], Clus['Cluster DEC'][j],
                                             Rad['Optical RA'][i], Rad['Optical DEC'][i])/3600)
            ClusdistNew = Clusoffset*comDist[i] # Find distance on sky at optical z
            # Check if this new distance is less than 15Mpc or the previous value
            if ClusdistNew < Clusdist:
                # If so, save it and its corresponding Delta z and index
                Clusdist = ClusdistNew
                ClusdzKeep = Clusdz
                Clusind = j
    # Once all j indices are run through for single i, check to see if any matches
    # were found that fit the constraints, and if so append it, but if not append nan
    if Clusdist < 15:
        ClusIDs = np.append(ClusIDs, Clus['Cluster ID'][Clusind])
        ClusallDz = np.append(ClusAllDz, ClusdzKeep)
        ClusallDist = np.append(ClusAllDist, Clusdist)
    else:
        ClusIDs = np.append(ClusIDs, np.nan)
        ClusallDz = np.append(ClusallDz, np.nan)
        ClusallDist = np.append(ClusallDist, np.nan)
            
                
                