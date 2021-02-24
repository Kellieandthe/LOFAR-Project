# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 11:41:27 2021

@author: ppykd1
"""

from astropy.table import Table
import LOFAR_Functions as lf
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kstest
import random

#%% Pre-processing
AGN = Table.read('Data\RLAGN angle data', format = 'fits')
AGNold = Table.read('Data\RLAGN angle data old', format = 'fits')
WHL = Table.read('Data\WHL angle data', format = 'fits')

# notAGN = WHL.copy()
# for i in np.arange(0, len(AGN)):
#     if np.isin(AGN['Radio Source'][i], notAGN['Radio Source']) == True:
#         notAGN.remove_row(i)

WAT = lf.WAT_cut(WHL)
NAT = lf.NAT_cut(WHL)
FR1 = lf.FR1_cut(WHL)
FR2 = lf.FR2_cut(WHL)
# ext = lf.ext_cut(WHL)
#%%
def offset_cond(Table):
    cutCond = (Table['Radio-Optical Offset'].data >= 5) &\
              (Table['Radio-Optical Offset'].data < 15)
    return Table[cutCond]

#%% Decide parameters
dat = NAT.copy()

dat = lf.Mpc_cut(dat, 2)
# dat = lf.delZ_cut(dat, 0.01)
dat = lf.BCG_cut(dat, 0.05)
# dat = lf.spec_cut(dat)
# # dat = lf.rich_cut(dat, 20)
# # dat = lf.ltflux_cut(dat, 500)
# # dat = lf.mtflux_cut(dat, 0.5)
# dat = offset_cond(dat)

def coords(ind):
    print([dat['Optical RA'][ind], dat['Optical DEC'][ind]])

#%% Angle
plt.close('all')

plt.figure()
plt.hist(dat['Angle ROC'], bins=18)
plt.xticks(np.arange(0, 200, 20))
plt.xlabel('Theta')
plt.ylabel('Number of sources')
plt.title('Angle distribution of radio sources in clusters')

#%% 2D Distance
plt.close('all')

plt.figure()
plt.hist(dat['2D Distance'], bins=50)
plt.xlabel('2D Distance (Mpc)')
plt.ylabel('Number of sources')
plt.title('Distance between cluster centre and optical source')
#plt.xticks(np.arange(0, 16, 1))
plt.tight_layout()

#%% Radio-optical offset
plt.close('all')

plt.figure()
plt.hist(dat['Radio-Optical Offset'], bins=50)
plt.xlabel('Radio-Optical Offset (arcsec)')
plt.ylabel('Number of sources')
plt.title('Distance between galaxy and radio source')
plt.tight_layout()

#%% Angle vs no. of r500s
plt.close('all')

plt.figure()
plt.scatter(dat['2D Distance']/dat['r500'], dat['Angle ROC'])
plt.ylabel('Theta')
plt.xlabel(r'2D Distance/$r_{500}$')
plt.yticks(np.arange(0, 200, 20))
plt.xlim([0, 3])

#%% Richness
plt.close('all')

plt.figure()
plt.hist(dat['Richness'], bins=20)
plt.xlabel('Richness')
plt.ylabel('Number of sources')
plt.title('Richness of clusters')

#%% Flux
plt.close('all')

plt.figure()
plt.hist(dat['Total_flux'], bins=20)
plt.xlabel('Total flux')
plt.ylabel('Number of sources')
plt.title('Total flux of radio sources')

#%% For checking that angles look right
plt.close('all')

n = 5
plt.figure()
plt.scatter(dat['Optical RA'][n], dat['Optical DEC'][n])
plt.scatter(dat['Radio RA'][n], dat['Radio DEC'][n])
plt.scatter(dat['Cluster RA'][n], dat['Cluster DEC'][n])
plt.text(dat['Optical RA'][n], dat['Cluster DEC'][n], r'Angle = %.2f' % (dat['Angle ROC'][n], ))
plt.legend(['Optical', 'Radio', 'Cluster'])

#%% k-s test
plt.close('all')

x = np.sort(dat['Angle ROC'].data)
y = np.arange(1, len(x)+1)/len(x)

xu = np.array([0, 180])
yu = np.array([0, 1])
xtest = np.linspace(0, 180, len(x), endpoint=True)

r = np.array([])
for i in np.arange(0, len(x)):
    r = np.append(r, random.uniform(0, 180))
    
r = np.sort(r)

x_alt = lf.alter_angle(x)
xtest_alt = lf.alter_angle(xtest)
r_alt = lf.alter_angle(r)

ts, p = kstest(x, xtest)
tsr, pr = kstest(x, r)
ts_alt, p_alt = kstest(x_alt, xtest_alt)
tsr_alt, pr_alt = kstest(x_alt, r_alt)

# Plot ks tests
plt.figure()
# Plot straight line
plt.plot(xu, yu)
# Plot uncut random data
plt.step(r, y)
# Plot uncut AGN data
plt.step(x, y)
# Include ks test values
textstr1 = '\n'.join((
    r'Straight line vs data test statistic = %.2f' % (ts, ),
    r'Straight line vs data p-value = %.2f' % (p, ),
    r'Random vs data test statistic = %.2f' % (tsr, ),
    r'Random vs data p-value = %.2f' % (p, )))
plt.text(0.05, 0.92, textstr1)
plt.xticks(np.arange(0, 200, 20))
plt.xlabel('Angle (deg)')
plt.ylabel('Cumulative probability')
plt.title('K-S test (regular order)')
plt.legend(['Straight line', 'Random data', 'AGN data'], loc='lower right')
plt.tight_layout()

## Plot histograms to go with k-s test ##

# Find mean and sd of data
mean = len(x)/18
sd = np.sqrt(mean)

# plt.close('all')
fig = plt.figure()
plt.subplot(1, 2, 1)
# Plot uncut angle data
plt.hist(x, bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('Regular order')
plt.xticks(np.arange(0, 200, 20))
plt.axhline(mean, color='darkgreen', linestyle='dashed', linewidth=1)
plt.axhline(mean-sd, color='darkorange', linestyle='dashed', linewidth=1)
plt.axhline(mean+sd, color='darkorange', linestyle='dashed', linewidth=1)

plt.subplot(1, 2, 2)
#Plot cut random data
plt.hist(r, bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('Random data (Same sample size)')
plt.xticks(np.arange(0, 200, 20))
plt.axhline(mean, color='darkgreen', linestyle='dashed', linewidth=1)
plt.axhline(mean-sd, color='darkorange', linestyle='dashed', linewidth=1)
plt.axhline(mean+sd, color='darkorange', linestyle='dashed', linewidth=1)
fig.legend(['mean', r'$1\sigma$'], loc='lower left')
plt.tight_layout()






