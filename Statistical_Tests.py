# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:16:14 2020

@author: ppykd1
"""

from astropy.table import Table, vstack
import LOFAR_Functions as lf
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kstest
from collections import Counter
from astropy.coordinates import SkyCoord
import random

RLAGN = Table.read('Data\RLAGN angle data old', format = 'fits')

#%% Cuts

cutDat = lf.delZ_cut(RLAGN, 0.01)
cutDat = lf.Mpc_cut(cutDat, 1.5*cutDat['r500'])
cutDat = lf.BCG_cut(cutDat, 0.01)


#%%

# Get list of all matched clusters in data set without any repeats
unique = np.array(list(Counter(RLAGN['Cluster ID']).keys()))
# Get list of how many times each cluster is repeated (how many galaxy matches it has)
count = np.array(list(Counter(RLAGN['Cluster ID']).values()))

# Only keep cluster with >= 3 galaxy matches
mult = unique[count >= 3]
num = count[count >= 3]

# Create a copy of the RLAGN table
cop = RLAGN.copy()

# Start by only keeping the first cluster in mult
cop = cop[cop['Cluster ID'] == mult[0]]
# for loop then runs through the rest of mult and appends the data for the rest of the clusters
for i in np.arange(1, len(mult)):
    temp = RLAGN[RLAGN['Cluster ID'] == mult[i]]
    cop = vstack([cop, temp])

# SkyCoords for all of the radio, optical and cluster co-ordinates in cop    
radCoord = SkyCoord(cop['Radio RA'], cop['Radio DEC'], unit='deg')
optCoord = SkyCoord(cop['Optical RA'], cop['Optical DEC'], unit='deg')
clusCoord = SkyCoord(cop['Cluster RA'], cop['Cluster DEC'], unit='deg')

plt.close('all')
fig = plt.figure()
# Plot for all 9 cluster systems
for i in np.arange(0, 9):
    j = 3*i
    if i >= 7: # The 7th cluster has 4 corresponding galaxies rather than 3
        j = 3*i + 1
    plt.subplot(3, 3, i+1)
    # plt.figure()
    plt.scatter(radCoord[j:j+num[i]].ra, radCoord[j:j+num[i]].dec)
    plt.scatter(optCoord[j:j+num[i]].ra, optCoord[j:j+num[i]].dec)
    plt.scatter(clusCoord[j].ra, clusCoord[j].dec)
    plt.xlabel('RA (deg)')
    plt.ylabel('DEC (deg)')
    # Finding max and min RA and DEC
    maxval = max(np.hstack((radCoord[j:j+num[i]].ra, optCoord[j:j+num[i]].ra, clusCoord[j].ra)).data)
    minval = min((np.hstack((radCoord[j:j+num[i]].ra, optCoord[j:j+num[i]].ra, clusCoord[j].ra))).data)
    plt.xlim([maxval+0.001, minval-0.001])
    plt.title('Cluster '+cop['Cluster ID'][j])
    # plt.legend(['Radio source', 'Optical source', 'Cluster centre'])
    # plt.tight_layout()
fig.legend(['Radio source', 'Optical source', 'Cluster centre'], loc='lower left')
plt.tight_layout()    

#%%

# Take just angle data
ang = RLAGN['Angle ROC'].data
# Cut to only keep data outside of BCG radius
angcut = lf.BCG_cut(RLAGN, 0.01)['Angle ROC'].data
# Sort data to put into cumulative distributions
x1 = np.sort(ang)
x2 = np.sort(angcut)
# Create y axis of K-S test graphs
y1 = np.arange(1, len(x1)+1)/len(x1)
y2 = np.arange(1, len(x2)+1)/len(x2)

# co-ordinates for straight line on K-S test graph
xu = np.array([0, 180])
yu = np.array([0, 1])

# Creat values to be put into K-S test
xtest1 = np.linspace(0, 180, len(x1), endpoint=True)
xtest2 = np.linspace(0, 180, len(x2), endpoint=True)

# Get random values within same range
r1 = np.array([])
for i in np.arange(0, len(x1)):
    r1 = np.append(r1, random.uniform(0, 180))

# Cut r1 down to create a random value array that's the same length as x2
r2 = r1[0:len(x2)]

# Sort random arrays
r1 = np.sort(r1)
r2 = np.sort(r2)

#%%
# Changing angle data for histograms
x2alt = x2.copy()
xtest2alt = xtest2.copy()
r2alt = r2.copy()
for i in np.arange(0, len(x2alt)):
    if x2alt[i] < 90:
        x2alt[i] = x2alt[i]+180
    if xtest2alt[i] < 90:
        xtest2alt[i] = xtest2alt[i]+180
    if r2alt[i] < 90:
        r2alt[i] = r2alt[i]+180
        
x2alt = np.sort(x2alt)
xtest2alt = np.sort(xtest2alt)
r2alt = np.sort(r2alt)
#%%

# Null hypothesis states that both cumulative distributions are similar
# Rejecting the null hypothesis means that the two cumulative distributions are different
# The test statistic is the largest vertical distance between the distributions
# The closer the TS is to 0, the more likely the two samples were drawn from the same distribution
# The p-value of 0.n means that there is an n% chance of obtaining a set of measurements at least this correlated
# if the underlying data is uncorrelated.
# It seems to be a common theme to reject the null hypothesis if p < 0.05.

# Perform ks tests for cut and uncut data vs straight line (uniform distribution)
ts1, p1 = kstest(x1, xtest1)
ts2, p2 = kstest(x2, xtest2)
# Perform ks tests for cut and uncut data vs random data
tsr1, pr1 = kstest(x1, r1)
tsr2, pr2 = kstest(x2, r2)

# Plot ks tests
plt.close('all')
plt.figure()
plt.subplot(1, 2, 1)
# Plot straight line
plt.plot(xu, yu)
# Plot uncut random data
plt.step(r1, y1)
# Plot uncut AGN data
plt.step(x1, y1)
plt.xticks(np.arange(0, 200, 20))
plt.xlabel('Angle (deg)')
plt.ylabel('Cumulative probability')
plt.title('Uncut AGN data')
# Include ks test values
textstr1 = '\n'.join((
    r'Uniform test statistic = %.2f' % (ts1, ),
    r'Uniform p-value = %.2f' % (p1, ),
    r'Random test statistic = %.2f' % (tsr1, ),
    r'Random p-value = %.2f' % (p1, )))
plt.text(0.05, 0.92, textstr1)
plt.legend(['Uniform distribution', 'Random data', 'AGN data'], loc='lower right')

plt.subplot(1, 2, 2)
# Plot straight line
plt.plot(xu, yu)
# Plot cut random data
plt.step(r2, y2)
# Plot cut AGN data
plt.step(x2, y2)
plt.xticks(np.arange(90, 290, 20))
plt.xlabel('Angle (deg)')
plt.ylabel('Cumulative probability')
plt.title('2D Distance > $0.01R_{500}$')
# Include ks test values
textstr2 = '\n'.join((
    r'Uniform test statistic = %.2f' % (ts2, ),
    r'Uniform p-value = %.2f' % (p2, ),
    r'Random test statistic = %.2f' % (tsr2, ),
    r'Random p-value = %.2f' % (p2, )))
plt.text(0.05, 0.92, textstr2)
plt.tight_layout()

#%%

# Find mean and sd of both cut and uncut data
mean1 = len(x1)/18
sd1 = np.sqrt(mean1)

mean2 = len(x2)/18
sd2 = np.sqrt(mean2)

# plt.close('all')
fig = plt.figure()
plt.subplot(2, 2, 1)
# Plot uncut angle data
plt.hist(ang, bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('Uncut AGN data')
plt.xticks(np.arange(0, 200, 20))
plt.axhline(mean1, color='k', linestyle='dashed', linewidth=1)
plt.axhline(mean1-sd1, color='r', linestyle='dashed', linewidth=1)
plt.axhline(mean1+sd1, color='r', linestyle='dashed', linewidth=1)

plt.subplot(2, 2, 2)
# Plot cut angle data
plt.hist(angcut, bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('2D Distance > $0.01R_{500}$')
plt.xticks(np.arange(0, 200, 20))
plt.axhline(mean2, color='k', linestyle='dashed', linewidth=1)
plt.axhline(mean2-sd2, color='r', linestyle='dashed', linewidth=1)
plt.axhline(mean2+sd2, color='r', linestyle='dashed', linewidth=1)

plt.subplot(2, 2, 3)
# Plot uncut random data
plt.hist(r1, bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('Random data (Same sample size as uncut AGN)')
plt.xticks(np.arange(0, 200, 20))
plt.axhline(mean1, color='k', linestyle='dashed', linewidth=1)
plt.axhline(mean1-sd1, color='r', linestyle='dashed', linewidth=1)
plt.axhline(mean1+sd1, color='r', linestyle='dashed', linewidth=1)

plt.subplot(2, 2, 4)
#Plot cut random data
plt.hist(r2, bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('Random data (Same sample size as cut AGN)')
plt.xticks(np.arange(0, 200, 20))
plt.axhline(mean2, color='k', linestyle='dashed', linewidth=1)
plt.axhline(mean2-sd2, color='r', linestyle='dashed', linewidth=1)
plt.axhline(mean2+sd2, color='r', linestyle='dashed', linewidth=1)
fig.legend(['mean', r'$1\sigma$'], loc='center')
plt.tight_layout()


