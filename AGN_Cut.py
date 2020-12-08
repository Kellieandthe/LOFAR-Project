# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 12:52:50 2020

@author: ppykd1
"""

from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
# import LOFAR_Functions as lf
from astropy.cosmology import FlatLambdaCDM

# Set cosmology according to Garon et al. paper
cosmo = FlatLambdaCDM(H0=68, Om0=0.315)

LOFAR = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/HETDEX associations and optical IDs.fits', format='fits') # len = 318520

# Flux cut
FC = LOFAR[LOFAR['Total_flux'] > 0.5] # len = 239845
# Cut for those with Optical IDs
O = LOFAR[LOFAR['ID_name'] != ''] # len = 231765
# Intersection of both above
FCO = FC[FC['ID_name'] != ''] # len = 172947
# Cut for all with redshift value
Z = LOFAR[np.isnan(LOFAR['z_best']) == False] # len = 162252

# Good photometric redshift test described in Hardcastle et al. paper
def z_test(Table):
    Delta_z = (Table['z1_max'] - Table['z1_min'])/2
    zTest = Delta_z/(1 + Table['z1_median']) # z1_median is all the photometric redshift values
    return zTest

# Cut for all with good photometric redshifts or spectroscopic redshifts
ZG = Z[(z_test(Z) < 0.1) | (np.isnan(Z['z_spec']) == False)] # len = 89769
FCOZ = FCO[np.isnan(FCO['z_best']) == False]
FCOZG = FCOZ[(z_test(FCOZ) < 0.1) | (np.isnan(FCOZ['z_spec']) == False)]

#%% Right side of Figure 2

# Find all sources which have both a spectroscopic and photometric redshift
bothZ = LOFAR[(np.isnan(LOFAR['z_spec']) == False) & (np.isnan(LOFAR['z1_median']) == False)]
# Find 'good' photometric redshifts from z test
GPhot = bothZ[(z_test(bothZ) < 0.1)]
# Find 'bad' photometric redshifts from z test
BPhot = bothZ[(z_test(bothZ) > 0.1)]

plt.close('all')
plt.scatter(np.log(1+BPhot['z_spec']), np.log(1+BPhot['z1_median']), s=7, c='darkslateblue') # Plot bad redshifts
plt.scatter(np.log(1+GPhot['z_spec']), np.log(1+GPhot['z1_median']), s=7, c='teal') # Plot good redshifts
plt.xticks(np.log(1+np.arange(6)), np.arange(6))
plt.yticks(np.log(1+np.arange(6)), np.arange(6))
plt.xlim([0, np.log(6.5)])
plt.ylim([0, np.log(6.5)])
plt.xlabel('Spectroscopic redshift')
plt.ylabel('Photometric redshift')
plt.legend(['Bad photometric', 'Good photometric'])


#%% Left side of Figure 2

# Find all sources from O with reasonable magnitudes
OWISE = O[O['w1Mag'] < 22]

# Find all sources from OWISE with no redshift value
WNoZ = OWISE[np.isnan(OWISE['z_best']) == True]
# Find all sources from OWISE with a photometric redshift
WPhot = OWISE[np.isnan(OWISE['z1_median']) == False] # z1_median is all the photometric redshift values
# Use z test to find 'good' photometric sources from WPhot
WGPhot = WPhot[z_test(WPhot) < 0.1]
# Find all sources from OWISE with a spectroscopic redshift
WSpec = OWISE[np.isnan(OWISE['z_spec']) == False]

plt.close('all')
plt.hist([WSpec['w1Mag'], WGPhot['w1Mag'], WPhot['w1Mag'], WNoZ['w1Mag']],
         color=['lawngreen', 'teal', 'darkslateblue', 'mediumpurple'], stacked=True, bins=300)
plt.xlabel('WISE band 1 magnitude')
plt.ylabel('Number of sources')
plt.xlim([12, 23])
plt.legend(['Spectroscopic', 'Good Photometric', 'Photometric', 'No redshift'])

# This still looks slightly different to the reference figure and I don't know why!!


#%%
M = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/gal_info_dr7_v5_2.fit.gz', format='fits')
M.remove_row(170713)



Mco = SkyCoord(M['RA'], M['DEC'], unit='deg')
Lco = SkyCoord(FCOZG['ID_ra'], FCOZG['ID_dec'], unit='deg')

idx, d2d, d3d = Lco.match_to_catalog_sky(Mco)
d2d = d2d.arcsecond

# plt.hist(d2d[d2d<=1])
idxCut = idx[d2d<=1]
d2dCut = d2d[d2d<=1]
Mnew = M[idxCut]

FCOZGM = FCOZG[d2d<=1]

Cut1 = FCOZG[(FCOZG['Ks_rest'] < (-17)) & (FCOZG['Ks_rest'] > (-33))]
#%%
Mnew['SUBCLASS'] = np.char.strip(Mnew['SUBCLASS'])

FCOZGM_RLAGN = FCOZGM[(Mnew['SUBCLASS'] == 'AGN') | (Mnew['SUBCLASS'] == 'AGN BROADLINE')]
FCOZGM_SFG = FCOZGM[(Mnew['SUBCLASS'] == 'STARFORMING') |\
                    (Mnew['SUBCLASS'] == 'STARFORMING BROADLINE')]
    
Cut2 = Cut1[np.isin(Cut1['Source_Name'], FCOZGM_SFG['Source_Name'], invert=True)]

S_obs = Cut2['Total_flux']*10**(-3)*10**(-26)
D_L = cosmo.luminosity_distance(Cut2['z_best'])*3.0857*10**22/u.Mpc
D_L = D_L.astype(float)
z = Cut2['z_best']
nu = 150*10**6
#alpha = np.log10(S_obs)/np.log10(nu)
alpha = 0.7
L = (S_obs*4*np.pi*D_L**2)/((1+z)**(1+alpha))

keep1 = np.isin(Cut2['Source_Name'], FCOZGM_RLAGN['Source_Name'])
keep2 = (L > 10**25) & (Cut2['Ks_rest'] > (-25))
keep3 = (Cut2['Ks_rest'] < (-25)) & (np.log10(L) > (25.3 - 0.06*(25 + Cut2['Ks_rest'])))
#%%
keepAll = keep1 | keep2 | keep3

#Cut3 = Cut2[Cut2['AllWISE'] ]

#%%

W1vega = FCOZG['w1Mag'] - 2.699
W2vega = FCOZG['w2Mag'] - 3.339
W3vega = FCOZG['w3Mag'] - 5.174
W4vega = FCOZG['w4Mag'] - 6.620

W2W3 = W2vega - W3vega
W1W2 = W1vega - W2vega

plt.plot(W2W3, W1W2, 'o', color='palegreen')
plt.xlim([-0.5, 5])
plt.ylim([-0.5, 2])





