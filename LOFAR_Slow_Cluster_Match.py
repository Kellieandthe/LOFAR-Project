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
WHL = Table.read('Full WHL Cluster Catalogue')
Rad = Table.read('Radio source data')

# Only keep radio sources that have a reasonable radio-optical offset
# Rad = Rad[(Rad['Radio-Optical Offset'] >= 1) & (Rad['Radio-Optical Offset'] <= 100)]

# Put RA and DECs into SkyCoord co ordinates
WHLco = SkyCoord(WHL['Cluster RA'], WHL['Cluster DEC'], unit='deg')
Radco = SkyCoord(Rad['Optical RA'], Rad['Optical DEC'], unit='deg')

# Calculate co-moving distance of all z values in radio-optical data
comDist = np.array(cosmo.comoving_distance(Rad['Optical z'])/u.Mpc)


#%%

# Create empty arrays to hold cluster IDs, Delta z values, and dist values
WHLIDs = np.array([])
WHLallDz = np.array([])
WHLallDist = np.array([])

# Run through all radio-optical values
for i in range(0, len(Rad['Optical RA'])):
    # Set dist values that are outside of the constraint of being < 15Mpc
    WHLdist = 15
    # Run through length of all WHL values
    for j in range(0, len(WHL['Cluster RA'])):
        # Calculate Delta z values between current optical source and cluster
        WHLdz = lf.z_diff_calc(Rad['Optical z'][i], WHL['Cluster z'][j])
        # See if absolute value of Delta z is less than 0.04
        if abs(WHLdz) < 0.04:
            # If so, calculate distance on sky between cluster and optical source
            # WHLoffset = np.radians((Radco[i].separation(WHLco[j])).degree) # Find offset between RA and DEC
            WHLoffset = np.radians(lf.offset(WHL['Cluster RA'][j], WHL['Cluster DEC'][j],
                                             Rad['Optical RA'][i], Rad['Optical DEC'][i])/3600)
            WHLdistNew = WHLoffset*comDist[i] # Find distance on sky at optical z
            # Check if this new distance is less than 15Mpc or the previous value
            if WHLdistNew < WHLdist:
                # If so, save it and its corresponding Delta z and index
                WHLdist = WHLdistNew
                WHLdzKeep = WHLdz
                WHLind = j
    # Once all j indices are run through for single i, check to see if any matches
    # were found that fit the constraints, and if so append it, but if not append nan
    if WHLdist < 15:
        WHLIDs = np.append(WHLIDs, WHL['Cluster ID'][WHLind])
        WHLallDz = np.append(WHLallDz, WHLdzKeep)
        WHLallDist = np.append(WHLallDist, WHLdist)
    else:
        WHLIDs = np.append(WHLIDs, np.nan)
        WHLallDz = np.append(WHLallDz, np.nan)
        WHLallDist = np.append(WHLallDist, np.nan)
            
                
                