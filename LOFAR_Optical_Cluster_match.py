# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 09:32:48 2020

@author: kelli
"""

from astropy.table import Table
# import numpy as np
import LOFAR_Functions as lf
# import math

RadDat = Table.read('Radio source data', format='fits')
ClusDat = Table.read('Cluster data', format='fits')

gal_z = RadDat['z']
galRA = RadDat['Opt RA']
galDec = RadDat['Opt Dec']

clusID = ClusDat['Cluster ID']
clusRA = ClusDat['RA']
clusDec = ClusDat['Dec']
clus_z = ClusDat['z']
clus_zSource = ClusDat['z Source']
clusCat = ClusDat['Catalogue']

delta_z = []
dist = []

MatchID = []
MatchRA = []
MatchDec = []
Matchz = []
MatchzSource = []
MatchCat = []
   

for i in range(0, len(gal_z)):
    z_vals = []
    d_vals = []
    ClusIDs = []
    ClusRAs = []
    ClusDecs = []
    Cluszs = []
    CluszSources = []
    ClusCats = []
    for j in range(0, len(clus_z)):
        z_val = lf.z_diff_calc(gal_z[i], clus_z[j])
        if abs(z_val) < 0.04:
            offset = lf.offset(clusRA[j], galRA[i], clusDec[j], galDec[i])
            d_val = lf.distGC_calc(clus_z[j], offset)
            if d_val < 15:
                z_vals.append(z_val)
                d_vals.append(d_val)
                ClusIDs.append(clusID[j])
                ClusRAs.append(clusRA[j])
                ClusDecs.append(clusDec[j])
                Cluszs.append(clus_z[j])
                CluszSources.append(clus_zSource[j])
                ClusCats.append(clusCat[j])
    if len(z_vals) == 0:
        dist.append('None')
        delta_z.append('None')
        MatchID.append('None')
        MatchRA.append('None')
        MatchDec.append('None')
        Matchz.append('None')
        MatchzSource.append('None')
        MatchCat.append('None')
    else:
        dist.append(min(d_vals))
        delta_z.append(z_vals[d_vals.index(min(d_vals))])
        MatchID.append(ClusIDs[d_vals.index(min(d_vals))])
        MatchRA.append(ClusRAs[d_vals.index(min(d_vals))])
        MatchDec.append(ClusDecs[d_vals.index(min(d_vals))])
        Matchz.append(Cluszs[d_vals.index(min(d_vals))])
        MatchzSource.append(CluszSources[d_vals.index(min(d_vals))])
        MatchCat.append(ClusCats[d_vals.index(min(d_vals))])
        

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

ClusMatchData.write('Cluster match data', format = 'fits')















