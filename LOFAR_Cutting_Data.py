# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 12:16:13 2020

@author: ppykd1
"""

from astropy.table import Table
# import LOFAR_Functions as lf
import matplotlib.pyplot as plt
import numpy as np
# from astropy.coordinates import SkyCoord


RMdat = Table.read('RM Angle data')
WHLdat = Table.read('WHL Angle data')

# Constrain radio-optical doffset to within 1 and 100 arcseconds
# Constrain redshifts to be > 0.05
# Constrain 2D distance on sky to be < 5Mpc
RMcutCond = (RMdat['Radio-Optical Offset'].data >= 1) &\
            (RMdat['Radio-Optical Offset'].data < 100) &\
            (RMdat['Optical z'].data > 0.05) &\
            (RMdat['2D Distance'].data < 5)
            
WHLcutCond = (WHLdat['Radio-Optical Offset'].data >= 1) &\
             (WHLdat['Radio-Optical Offset'].data < 100) &\
             (WHLdat['Optical z'].data > 0.05) &\
             (WHLdat['2D Distance'].data < 5)

# Apply constraints
RMdatCut = RMdat[RMcutCond]
WHLdatCut = WHLdat[WHLcutCond]

plt.close('all')

# Plot ROC angles with new constraints
plt.subplot(2, 1, 1)
plt.hist(WHLdatCut['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 1, 2)
plt.hist(RMdatCut['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between Cluster-Optical-Radio sources')
plt.tight_layout()

# Import morphological data of sample of radio sources
MorphDat = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/Catalogue of morphological classifications.fits', format='fits')

"""
I'm sure there is a much better way to go about this than a nested for loop, but
I can't figure out what it is! I need to find common Radio Source IDs between
my now cut down tables and the morphological data table, so that I can then add
boolean columns onto the cut tables to show whether the sources are FR1, FR2, etc.
Hence, that's what the code below is doing. Please let me know if you know how
to do this more efficiently!
"""

RM_FR1 = np.full(len(RMdatCut), False)
RM_FR2 = np.full(len(RMdatCut), False)
RM_Indet = np.full(len(RMdatCut), False)
RM_Small = np.full(len(RMdatCut), False)
RM_NAT = np.full(len(RMdatCut), False)
RM_WAT = np.full(len(RMdatCut), False)
RM_DD = np.full(len(RMdatCut), False)

for i in np.arange(0, len(RMdatCut)):
    for j in np.arange(0, len(MorphDat)):
        if RMdatCut['Radio Source'][i] == MorphDat['Source_Name'][j]:
            np.put(RM_FR1, i, MorphDat['FR1'][j])
            np.put(RM_FR2, i, MorphDat['FR2'][j])
            np.put(RM_Indet, i, MorphDat['Indeterminate'][j])
            np.put(RM_Small, i, MorphDat['Small'][j])
            np.put(RM_NAT, i, MorphDat['NAT'][j])
            np.put(RM_WAT, i, MorphDat['WAT'][j])
            np.put(RM_DD, i, MorphDat['D-D'][j])
            
WHL_FR1 = np.full(len(WHLdatCut), False)
WHL_FR2 = np.full(len(WHLdatCut), False)
WHL_Indet = np.full(len(WHLdatCut), False)
WHL_Small = np.full(len(WHLdatCut), False)
WHL_NAT = np.full(len(WHLdatCut), False)
WHL_WAT = np.full(len(WHLdatCut), False)
WHL_DD = np.full(len(WHLdatCut), False)

for i in np.arange(0, len(WHLdatCut)):
    for j in np.arange(0, len(MorphDat)):
        if WHLdatCut['Radio Source'][i] == MorphDat['Source_Name'][j]:
            np.put(WHL_FR1, i, MorphDat['FR1'][j])
            np.put(WHL_FR2, i, MorphDat['FR2'][j])
            np.put(WHL_Indet, i, MorphDat['Indeterminate'][j])
            np.put(WHL_Small, i, MorphDat['Small'][j])
            np.put(WHL_NAT, i, MorphDat['NAT'][j])
            np.put(WHL_WAT, i, MorphDat['WAT'][j])
            np.put(WHL_DD, i, MorphDat['D-D'][j])

# Add new columns onto cut data tables
RMdatCut.add_columns([RM_FR1, RM_FR2, RM_Indet, RM_Small, RM_NAT, RM_WAT, RM_DD],
                     names=['FR1', 'FR2', 'Indeterminate', 'Small', 'NAT', 'WAT', 'D-D'])

WHLdatCut.add_columns([WHL_FR1, WHL_FR2, WHL_Indet, WHL_Small, WHL_NAT, WHL_WAT, WHL_DD],
                     names=['FR1', 'FR2', 'Indeterminate', 'Small', 'NAT', 'WAT', 'D-D'])

#%%

WHLWATCond = WHLdatCut['WAT'] == True
RMWATCond = RMdatCut['WAT'] == True

# Plot ROC angles with new constraints
plt.subplot(2, 1, 1)
plt.hist(WHLdatCut['Angle ROC'][WHLWATCond], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 1, 2)
plt.hist(RMdatCut['Angle ROC'][RMWATCond], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between Cluster-Optical-Radio sources')
plt.tight_layout()






