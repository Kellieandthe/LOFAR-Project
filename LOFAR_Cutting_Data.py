# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 12:16:13 2020

@author: ppykd1
"""

from astropy.table import Table, join
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
            (RMdat['3D Distance'].data <= 5)
            
WHLcutCond = (WHLdat['Radio-Optical Offset'].data >= 1) &\
             (WHLdat['Radio-Optical Offset'].data < 100) &\
             (WHLdat['Optical z'].data > 0.05) &\
             (WHLdat['3D Distance'].data <= 5) &\
             (WHLdat['3D Distance'].data > 0.01*WHLdat['r500'])
             
WHLBCGcut = WHLdat['3D Distance'].data <= 0.01*WHLdat['r500']

# Apply constraints
RMdatCut = RMdat[RMcutCond]
WHLdatCut = WHLdat[WHLcutCond]
WHLBCGs = WHLdat[WHLBCGcut]


# Import morphological data of sample of radio sources
MorphDat = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/Catalogue of morphological classifications.fits', format='fits')
MorphDat.remove_columns(['RA', 'DEC', 'LM_dec_size', 'LM_Flux'])
MorphDat.rename_column('Source_Name', 'Radio Source')
# Import LOFAR data to take AGN columns
AllRad = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/HETDEX associations and optical IDs.fits', format='fits')
AllRad.keep_columns(['Source_Name', 'specAGN', 'mqcAGN', 'XrayClass'])
AllRad.rename_column('Source_Name', 'Radio Source')

# Find radio sources in current tables that have morphological data and add those columns onto the end
RMMorphDat = join(RMdatCut, MorphDat, keys='Radio Source', join_type='left')
WHLMorphDat = join(WHLdatCut, MorphDat, keys='Radio Source', join_type='left')

# Check that the matched sources in the tables look right (they do!)
RMmatch = set(RMdatCut['Radio Source']) & set(MorphDat['Radio Source'])
RMmatch = sorted(RMmatch)
WHLmatch = set(WHLdatCut['Radio Source']) & set(MorphDat['Radio Source'])
WHLmatch = sorted(WHLmatch)

RMMorphDat = join(RMMorphDat, AllRad, keys='Radio Source', join_type='left')
WHLMorphDat = join(WHLMorphDat, AllRad, keys='Radio Source', join_type='left')

# Since added columns are boolean, set any with -- to False, and write to fits tables
# (RMMorphDat.filled(False)).write('RM morphological data', format = 'fits')
# (WHLMorphDat.filled(False)).write('WHL morphological data', format = 'fits')

#%%
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
plt.suptitle('Angle between Radio-Optical-Cluster sources')
plt.tight_layout()


WHLNATCond = WHLMorphDat['NAT'] == True
RMNATCond = RMMorphDat['NAT'] == True

# Plot ROC angles with new constraints for NATs
plt.figure()
plt.subplot(2, 1, 1)
plt.hist(WHLdatCut['Angle ROC'][WHLNATCond], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 1, 2)
plt.hist(RMdatCut['Angle ROC'][RMNATCond], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between Radio-Optical-Cluster sources for NATs')
plt.tight_layout()

WHLAGNCond = (WHLMorphDat['specAGN'] == 1.0) | (WHLMorphDat['mqcAGN'] == True) |\
             (WHLMorphDat['XrayClass'] == 1.0)
RMAGNCond = (RMMorphDat['specAGN'] == 1.0) | (RMMorphDat['mqcAGN'] == True) |\
             (RMMorphDat['XrayClass'] == 1.0)

# Plot ROC angles with new constraints for AGN
plt.figure()
plt.subplot(2, 1, 1)
plt.hist(WHLdatCut['Angle ROC'][WHLAGNCond], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 1, 2)
plt.hist(RMdatCut['Angle ROC'][RMAGNCond], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between Radio-Optical-Cluster sources for AGN')
plt.tight_layout()


#%% Plot richness histograms
plt.figure()
plt.subplot(2, 1, 1)
plt.hist(WHLdat['Richness'], bins=18)
plt.xlabel('Richness')
plt.ylabel('Number of sources')
plt.title('WHL15')

plt.subplot(2, 1, 2)
plt.hist(RMdat['Richness'], bins=18)
plt.xlabel('Richness')
plt.ylabel('Number of sources')
plt.title('redMaPPer')
plt.suptitle('Richness of matched galaxy clusters')
plt.tight_layout()

# Plot angle against 2D distance from cluster centre
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(WHLdatCut['Angle ROC'], WHLdatCut['2D Distance'], 'o')
plt.xlabel('Angle (deg)')
plt.ylabel('2D Distance (Mpc)')
plt.title('WHL15')

plt.subplot(2, 1, 2)
plt.plot(RMdatCut['Angle ROC'], RMdatCut['2D Distance'],'o')
plt.xlabel('Angle (deg)')
plt.ylabel('2D Distance (Mpc)')
plt.title('redMaPPer')
plt.suptitle('2D Distance of radio galaxy from cluster centre against ROC angle')
plt.tight_layout()

# Plot normalised distance from cluster centre
normDist = np.log10(WHLdat['2D Distance']/WHLdat['r500'])

plt.figure()
plt.hist(normDist, bins=15)
plt.xlabel(r'$log(\frac{Dist_{2D}}{r_{500}})$')
plt.ylabel('Count')
plt.title('Normalised 2D Distance of radio galaxy from cluster centre (WHL15)')

# Plot 2D distance against radio-optical offset
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(WHLdat['2D Distance'], WHLdat['Radio-Optical Offset'], 'o')
plt.xlabel('2D Distance (Mpc)')
plt.ylabel('Offset (arcsec)')
plt.title('WHL15')

plt.subplot(2, 1, 2)
plt.plot(RMdat['2D Distance'], RMdat['Radio-Optical Offset'], 'o')
plt.xlabel('2D Distance (Mpc)')
plt.ylabel('Offset (arcsec)')
plt.title('redMaPPer')
plt.suptitle('2D Distance of radio galaxy from cluster centre against radio-optical offset')
plt.tight_layout()

# Plot angle against optical mass proxy
plt.figure()
plt.plot(WHLdat['Angle ROC'], WHLdat['M500'], 'o')
plt.xlabel('Angle (deg)')
plt.ylabel(r'$M_{500}$')
plt.title('ROC angle against Optical Mass Proxy (WHL15)')
