# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 13:33:17 2020

@author: -
"""

from astropy.table import Table, join
import numpy as np
import matplotlib.pyplot as plt


# Import redMaPPer cluster catalogue data as table
RedMap = Table.read('Data\Full RM Cluster Catalogue', format='fits')

# Put cluster names into array
RMID = np.array([])

# Put cluster names into a form where they can be compared with the WHL15 cluster catalogue names
for i in range(0, len(RedMap)):
    if i%1000 == 0:
        print(i)
    oldName = RedMap['Cluster ID'][i]
    newName = oldName[2:] # remove 'RM' at beginning
    dec = str(int(round(float(newName[10:]),0))) # take 'dec' section of name which has one d.p. and round to nearest integer
    if len(dec) < 6: # make sure length of dec is correct
        l = 6 - len(dec) 
        dec = l*'0' + dec # add any preceding 0s at beginning of dec that may have been removed in rounding process
    RMID = np.append(RMID, (newName[:10] + dec)) # add new rounded dec back onto cluster name

# Import WHL cluster data
WHL = Table.read('Data\Full WHL Cluster Catalogue', format='fits')      

WHLID = np.array([])

# Adjust names to be able to compare with redMaPPer        
for i in range(0, len(WHL)):
    if i%1000 == 0:
        print(i)
    if i < 132684:
        oldName = WHL['Cluster ID'][i]
        newName = oldName[4:] # remove WHL and space at beginning of name
        if newName[-1] == '*': # remove asterisk at end of name which signifies the cluster data has been updated in the new catalogue
            WHLID = np.append(WHLID, newName[:16])
        else:
            WHLID = np.append(WHLID, newName)
    else:
        oldName = WHL['Cluster ID'][i]
        WHLID = np.append(WHLID, oldName[3:]) # remove WH at beginning of name
        
dupRM = np.array([])
dupWHL = np.array([])
dupWHL = dupWHL.astype(str)

# Look for duplicates in the WHL15 and redMaPPer data, and take the WHL15 if there is a duplicate
for i in range(0, len(RMID)):
    cond = WHLID == RMID[i]
    if len(WHLID[cond]) == 1:
        dupRM = np.append(dupRM, RedMap['Cluster ID'][i])
        dupWHL = np.append(dupWHL, WHL['Cluster ID'][np.where(cond == True)[0]])
        

dupTab = Table([dupRM, dupWHL], names=['RM Cluster ID', 'WHL Cluster ID'])

WHL.rename_columns(['Cluster ID', 'Richness'], ['WHL Cluster ID', 'WHL Richness'])
RedMap.rename_columns(['Cluster ID', 'Richness'], ['RM Cluster ID', 'RM Richness'])

newTab = join(dupTab, WHL, keys='WHL Cluster ID', join_type='left')
newTab2 = join(newTab, RedMap, keys = 'RM Cluster ID', join_type='left')

#%%
plt.scatter(newTab2['RM Richness'], newTab2['WHL Richness'], s=7)
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$R_{L*}$')
plt.xticks(range(25, 225, 25))
plt.yticks(range(0, 250, 50))
plt.xlim([13, 225])
plt.ylim([0, 250])
plt.title('Richness estimators for duplicate sources in WHL15 and RM14 cluster catalogues')






        
        