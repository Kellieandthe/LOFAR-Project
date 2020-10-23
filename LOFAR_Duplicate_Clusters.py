# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 13:33:17 2020

@author: -
"""

from astropy.table import Table
import numpy as np


RedMap = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/Redmapper DR8 Public v6.3 Cluster Catalog.fits.gz', format='fits')

RMIDOld = RedMap['NAME']
RMID = []

for i in range(0, len(RMIDOld)):
    oldName = RMIDOld[i]
    newName = oldName[2:]
    dec = str(int(round(float(newName[10:]),0)))
    if len(dec) < 6:
        l = 6 - len(dec)
        dec = l*'0' + dec
    RMID.append((newName[:10] + dec))

RMRA = RedMap['RA']
RMDec = RedMap['DEC']
RMZPhoto = RedMap['Z_LAMBDA']
RMZSpec = RedMap['Z_SPEC']
RM_z = []
RM_zSource = []

for i in range(0, len(RMZPhoto)):
    if RMZSpec[i] == -1:
        RM_z.append(RMZPhoto[i])
        RM_zSource.append('Photometric')
    else:
        RM_z.append(RMZSpec[i])
        RM_zSource.append('Spectroscopic')
        

WHL15 = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/The WHL12 cluster catalog with updated parameters.txt', format='cds')
WHL15New = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/Newly identified rich clusters at high redshifts.txt', format='cds')

WHL15IDOld = WHL15['Name'].compressed()
WHL15NewIDOld = WHL15New['Name'].compressed()
WHL15ID = []
WHL15NewID = []         
        
for i in range(0, len(WHL15IDOld)):
    oldName = WHL15IDOld[i]
    newName = oldName[4:]
    if newName[-1] == '*':
        WHL15ID.append(newName[:16])
    else:
        WHL15ID.append(newName)
        
for i in range(0, len(WHL15NewIDOld)):
    oldName = WHL15NewIDOld[i]
    WHL15NewID.append(oldName[3:])

WHLID = WHL15ID + WHL15NewID
WHLRA = np.concatenate((WHL15['RAdeg'].compressed(), WHL15New['RAdeg'].compressed()))
WHLDec = np.concatenate((WHL15['DEdeg'].compressed(), WHL15New['DEdeg'].compressed()))

WHLZPhoto = WHL15['zphot'].compressed()
WHLZSpec = WHL15['zspec'].compressed()
WHLZOld = []
WHL_zSource = []

for i in range(0, len(WHLZPhoto)):
    if WHLZSpec[i] == -1:
        WHLZOld.append(WHLZPhoto[i])
        WHL_zSource.append('Photometric')
    else:
        WHLZOld.append(WHLZSpec[i])
        WHL_zSource.append('Spectroscopic')
        
WHL_z = np.concatenate((WHLZOld, WHL15New['zspec'].compressed()))

for i in range(0, len(WHL15New['zspec'].compressed())):
    WHL_zSource.append('Spectroscopic')
    

ClusID = []
ClusRA = []
ClusDec = []
Clus_z = []
Clus_zSource = []
ClusCat = []

for i in range(0, len(WHLID)):
    ClusID.append(WHLID[i])
    ClusRA.append(WHLRA[i])
    ClusDec.append(WHLDec[i])
    Clus_z.append(WHL_z[i])
    Clus_zSource.append(WHL_zSource[i])
    ClusCat.append('WHL15')

for i in range(0, len(RMID)):
    duplicate = False
    for j in range(0, len(ClusID)):
        if RMID[i] == ClusID[j]:
            duplicate = True
    if duplicate == False:
        ClusID.append(RMID[i])
        ClusRA.append(RMRA[i])
        ClusDec.append(RMDec[i])
        Clus_z.append(RM_z[i])
        Clus_zSource.append(RM_zSource[i])
        ClusCat.append('redMaPPer')
        

AllClusData = Table([ClusID, ClusRA, ClusDec, Clus_z, Clus_zSource, ClusCat],
                    names = ('Cluster ID', 'RA', 'Dec', 'z', 'z Source',
                             'Catalogue'))
AllClusData['RA'].unit = 'deg'
AllClusData['Dec'].unit = 'deg'

AllClusData.write('Cluster data', format='fits')
