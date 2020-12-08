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


RMdatG = Table.read('Data\RM Angle data (Garon)')
RMdat3D = Table.read('Data\RM Angle data (3D)')
WHLdatG = Table.read('Data\WHL Angle data (Garon)')
WHLdat3D = Table.read('Data\WHL Angle data (3D)')

# Constrain radio-optical offset to within 1 and 100 arcseconds
# Constrain redshifts to be > 0.05
# Constrain 2D distance on sky to be < 5Mpc
def cut_cond(Table):
    cutCond = (Table['Radio-Optical Offset'].data >= 1) &\
              (Table['Radio-Optical Offset'].data < 100) &\
              (Table['z_best'].data > 0.05) &\
              (Table['2D Distance'].data <= 5)
    return Table[cutCond]
              
# Define BCG region - exclusive to WHL matches
def BCG_cut(Table):
    BCGcut = Table['2D Distance'].data <= 0.01*Table['r500']
    return Table[BCGcut]
            

# Apply constraints
RMdatGCut = cut_cond(RMdatG)
RMdat3DCut = cut_cond(RMdat3D)
WHLdatGCut = cut_cond(WHLdatG)
WHLdat3DCut = cut_cond(WHLdat3D)
WHLdatG_BCGs = BCG_cut(WHLdatG)
WHLdat3D_BCGs = BCG_cut(WHLdat3D)


# Import morphological data of sample of radio sources
MorphDat = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/Catalogue of morphological classifications.fits', format='fits')
MorphDat.remove_columns(['RA', 'DEC', 'LM_dec_size', 'LM_Flux'])
MorphDat.rename_column('Source_Name', 'Radio Source')

# Find radio sources in current tables that have morphological data and add those columns onto the end
RMMorphDatG = join(RMdatGCut, MorphDat, keys='Radio Source', join_type='left')
RMMorphDat3D = join(RMdat3DCut, MorphDat, keys='Radio Source', join_type='left')
WHLMorphDatG = join(WHLdatGCut, MorphDat, keys='Radio Source', join_type='left')
WHLMorphDat3D = join(WHLdat3DCut, MorphDat, keys='Radio Source', join_type='left')

# # Check that the matched sources in the tables look right (they do!)
# RMmatch = set(RMdatCut['Radio Source']) & set(MorphDat['Radio Source'])
# RMmatch = sorted(RMmatch)
# WHLmatch = set(WHLdatCut['Radio Source']) & set(MorphDat['Radio Source'])
# WHLmatch = sorted(WHLmatch)

# Since added columns are boolean, set any with -- to False, and write to fits tables
# (RMMorphDatG.filled(False)).write('RM cut morphological data (Garon)', format = 'fits')
# (RMMorphDat3D.filled(False)).write('RM cut morphological data (3D)', format = 'fits')
# (WHLMorphDatG.filled(False)).write('WHL cut morphological data (Garon)', format = 'fits')
# (WHLMorphDat3D.filled(False)).write('WHL cut morphological data (3D)', format = 'fits')

#%%
plt.close('all')

# Plot ROC angles with new constraints
plt.subplot(2, 2, 1)
plt.hist(WHLMorphDat3D['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15 (3D match)')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 2, 2)
plt.hist(RMMorphDat3D['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer (3D match)')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 2, 3)
plt.hist(WHLMorphDatG['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15 (Garon et al. match)')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 2, 4)
plt.hist(RMMorphDatG['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer (Garon et al. match)')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between cut Radio-Optical-Cluster sources')
plt.tight_layout()

def NAT_cut(Table):
    NATcut = Table['NAT'] == True
    return Table[NATcut]

def WAT_cut(Table):
    WATcut = Table['WAT'] == True
    return Table[WATcut]

# Plot ROC angles with new constraints for NATs
plt.figure()
plt.subplot(2, 2, 1)
plt.hist(NAT_cut(WHLMorphDat3D)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15 (3D match)')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 2, 2)
plt.hist(NAT_cut(RMMorphDat3D)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer (3D match)')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 2, 3)
plt.hist(NAT_cut(WHLMorphDatG)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15 (Garon et al. match)')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 2, 4)
plt.hist(NAT_cut(RMMorphDatG)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer (Garon et al. match)')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between cut Radio-Optical-Cluster sources for NATs')
plt.tight_layout()


# Plot ROC angles with new constraints for WATs
plt.figure()
plt.subplot(2, 2, 1)
plt.hist(WAT_cut(WHLMorphDat3D)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15 (3D match)')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 2, 2)
plt.hist(WAT_cut(RMMorphDat3D)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer (3D match)')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 2, 3)
plt.hist(WAT_cut(WHLMorphDatG)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15 (Garon et al. match)')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 2, 4)
plt.hist(WAT_cut(RMMorphDatG)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer (Garon et al. match)')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between cut Radio-Optical-Cluster sources for WATs')
plt.tight_layout()

def AGN_cut(Table):
    AGNcut = (Table['specAGN'] == 1.0) |\
             (Table['mqcAGN'] == True) |\
             (Table['XrayClass'] == 1.0)
    return Table[AGNcut]

# Plot ROC angles with new constraints for AGN
plt.figure()
plt.subplot(2, 2, 1)
plt.hist(AGN_cut(WHLMorphDat3D)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15 (3D match)')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 2, 2)
plt.hist(AGN_cut(RMMorphDat3D)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer (3D match)')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 2, 3)
plt.hist(AGN_cut(WHLMorphDatG)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15 (Garon et al. match)')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 2, 4)
plt.hist(AGN_cut(RMMorphDatG)['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer (Garon et al. match)')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between cut Radio-Optical-Cluster sources for LOFAR Classified AGN')
plt.tight_layout()


#%% Plot richness histograms
plt.figure()
plt.subplot(2, 2, 1)
plt.hist(WHLMorphDat3D['Richness'], bins=18)
plt.xlabel('Richness')
plt.ylabel('Number of sources')
plt.title('WHL15 (3D match)')

plt.subplot(2, 2, 2)
plt.hist(RMMorphDat3D['Richness'], bins=18)
plt.xlabel('Richness')
plt.ylabel('Number of sources')
plt.title('redMaPPer (3D match)')

plt.subplot(2, 2, 3)
plt.hist(WHLMorphDatG['Richness'], bins=18)
plt.xlabel('Richness')
plt.ylabel('Number of sources')
plt.title('WHL15 (Garon et al. match)')

plt.subplot(2, 2, 4)
plt.hist(RMMorphDatG['Richness'], bins=18)
plt.xlabel('Richness')
plt.ylabel('Number of sources')
plt.title('redMaPPer (Garon et al. match)')
plt.suptitle('Richness of cut matched galaxy clusters')
plt.tight_layout()

# Plot normalised distance from cluster centre
# normDist = np.log10(WHLdat['2D Distance']/WHLdat['r500'])

# plt.figure()
# plt.hist(normDist, bins=15)
# plt.xlabel(r'$log(\frac{Dist_{2D}}{r_{500}})$')
# plt.ylabel('Count')
# plt.title('Normalised 2D Distance of radio galaxy from cluster centre (WHL15)')







