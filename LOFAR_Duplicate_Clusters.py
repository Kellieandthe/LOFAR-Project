# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 13:33:17 2020

@author: -
"""

from astropy.table import Table
import numpy as np
import time

start_time = time.time()

# Import redMaPPer cluster catalogue data as table
RedMap = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/Redmapper DR8 Public v6.3 Cluster Catalog.fits.gz', format='fits')

# Put cluster names into array
RMIDOld = RedMap['NAME']
RMID = np.array([])

# Put cluster names into a form where they can be compared with the WHL15 cluster catalogue names
for i in range(0, len(RMIDOld)):
    oldName = RMIDOld[i]
    newName = oldName[2:] # remove 'RM' at beginning
    dec = str(int(round(float(newName[10:]),0))) # take 'dec' section of name which has one d.p. and round to nearest integer
    if len(dec) < 6: # make sure length of dec is correct
        l = 6 - len(dec) 
        dec = l*'0' + dec # add any preceding 0s at beginning of dec that may have been removed in rounding process
    np.append(RMID, (newName[:10] + dec)) # add new rounded dec back onto cluster name

RMID = np.array(RMID)

# Put other helpful cluster values into arrays
RMRA = RedMap['RA'].data
RMDec = RedMap['DEC'].data
RMZPhoto = RedMap['Z_LAMBDA'].data
RMZSpec = RedMap['Z_SPEC'].data
RMRich = RedMap['LAMBDA'].data

RM_z = np.array([])
RM_zSource = np.array([])

# -1 in the spectroscopic redshift column means that there is no value, so photometric redshift is taken instead in this case
for i in range(0, len(RMZPhoto)):
    if RMZSpec[i] == -1:
        np.append(RM_z, RMZPhoto[i])
        np.append(RM_zSource, 'Photometric') # keep track of whether the taken redshift is spectroscopic or photometric
    else:
        np.append(RM_z, RMZSpec[i])
        np.append(RM_zSource, 'Spectroscopic')

# Convert to arrays
RM_z = np.array(RM_z)
RM_zSource = np.array(RM_zSource)

# Import WHL15 cluster data. Updated data from the WHL12 catalogue is in WHL15, and newly identified clusters are in WHL15New
WHL15 = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/The WHL12 cluster catalog with updated parameters.txt', format='cds')
WHL15New = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/Newly identified rich clusters at high redshifts.txt', format='cds')

# Put cluster IDs into arrays to be able to compare with redMaPPer to find duplicates
WHL15IDOld = (WHL15['Name'].compressed()).data
WHL15NewIDOld = (WHL15New['Name'].compressed()).data
WHL15ID = np.array([])
WHL15NewID = np.array([])         

# Adjust names to be able to compare with redMaPPer        
for i in range(0, len(WHL15IDOld)):
    oldName = WHL15IDOld[i]
    newName = oldName[4:] # remove WHL and space at beginning of name
    if newName[-1] == '*': # remove asterisk at end of name which signifies the cluster data has been updated in the new catalogue
        np.append(WHL15ID, newName[:16])
    else:
        np.append(WHL15ID, newName)
        
for i in range(0, len(WHL15NewIDOld)):
    oldName = WHL15NewIDOld[i]
    WHL15NewID.append(oldName[3:]) # remove WHL at beginning of name

# Combine old and new catalogue data into single arrays
WHLID = np.array(WHL15ID + WHL15NewID)
WHLRA = np.concatenate((WHL15['RAdeg'].compressed(), WHL15New['RAdeg'].compressed()))
WHLDec = np.concatenate((WHL15['DEdeg'].compressed(), WHL15New['DEdeg'].compressed()))
WHLRich = np.concatenate((WHL15['RL*500'].compressed(), WHL15New['RL*500'].compressed()))

WHLZPhoto = (WHL15['zphot'].compressed()).data
WHLZSpec = (WHL15['zspec'].compressed()).data
WHLZOld = np.array([])
WHL_zSource = np.array([])

# Same as for redMaPPer, -1 signifies no spectroscopic redshift, so photometric is taken instead
for i in range(0, len(WHLZPhoto)):
    if WHLZSpec[i] == -1:
        np.append(WHLZOld, WHLZPhoto[i])
        np.append(WHL_zSource, 'Photometric') # keep track of which type of redshift
    else:
        np.append(WHLZOld, WHLZSpec[i])
        np.append(WHL_zSource, 'Spectroscopic')

# Combine old and new catalogue redshifts
WHL_z = np.concatenate((WHLZOld, WHL15New['zspec'].compressed())) # all new catalogue clusters have spectroscopic redshifts

for i in range(0, len(WHL15New['zspec'].compressed())):
    np.append(WHL_zSource, 'Spectroscopic')

# Convert to array    
WHL_zSource = np.array(WHL_zSource)
    
# Start putting all cluster data into single arrays to put into table
ClusID = WHLID
ClusRA = WHLRA
ClusDec = WHLDec
Clus_z = WHL_z
Clus_zSource = WHL_zSource
ClusCat = np.full(len(ClusID), 'WHL15')
ClusRich = WHLRich

# Look for duplicates in the WHL15 and redMaPPer data, and take the WHL15 if there is a duplicate
for i in range(0, len(RMID)):
    duplicate = False
    for j in range(0, len(ClusID)):
        if RMID[i] == ClusID[j]:
            duplicate = True
    if duplicate == False:
        np.append(ClusID, RMID[i])
        np.append(ClusRA, RMRA[i])
        np.append(ClusDec, RMDec[i])
        np.append(Clus_z, RM_z[i])
        np.append(Clus_zSource, RM_zSource[i])
        np.append(ClusCat, 'redMaPPer')
        np.append(ClusRich, RMRich[i])
        
# Put all relevant data into a table
AllClusData = Table([ClusID, ClusRA, ClusDec, Clus_z, Clus_zSource, ClusCat, ClusRich],
                    names = ('Cluster ID', 'RA', 'Dec', 'z', 'z Source',
                             'Catalogue', 'Richness'))
AllClusData['RA'].unit = 'deg'
AllClusData['Dec'].unit = 'deg'

# Write data to a FITS file
AllClusData.write('Cluster data', format='fits')

print("This program took", time.time() - start_time, "to run")