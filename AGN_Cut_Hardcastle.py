# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 12:52:50 2020

@author: ppykd1
"""

from astropy.table import Table, hstack, vstack, setdiff
import numpy as np
import matplotlib.pyplot as plt
# from astropy import units as u
from astropy.coordinates import SkyCoord
import LOFAR_Functions as lf
from astropy.cosmology import FlatLambdaCDM

# Set cosmology according to Garon et al. paper
cosmo = FlatLambdaCDM(H0=68, Om0=0.315)

LOFAR = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/HETDEX associations and optical IDs.fits', format='fits') # len = 318520

# Add columns with Vega WISE values
lf.Vega_conv(LOFAR)

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
plt.figure()
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

# plt.close('all')
plt.figure()
plt.hist([WSpec['w1Mag'], WGPhot['w1Mag'], WPhot['w1Mag'], WNoZ['w1Mag']],
         color=['lawngreen', 'teal', 'darkslateblue', 'mediumpurple'], stacked=True, bins=300)
plt.xlabel('WISE band 1 magnitude')
plt.ylabel('Number of sources')
plt.xlim([12, 23])
plt.legend(['Spectroscopic', 'Good Photometric', 'Photometric', 'No redshift'])


#%% Defining FCOZGM

# Import MPA-JHU data
M = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/gal_info_dr7_v5_2.fit.gz', format='fits')
MBS = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/gal_idxfix_dr7_v5_2.fit.gz', format='fits')
MSM = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/totlgm_dr7_v5_2.fit.gz', format='fits')
MEL = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/gal_line_dr7_v5_2.fit.gz', format='fits')
MSFR = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/gal_totsfr_dr7_v5_2.fits.gz', format='fits')
# This row seems to have dud data
M.remove_row(170713)
MBS.remove_row(170713)
MSM.remove_row(170713)
MEL.remove_row(170713)
MSFR.remove_row(170713)
#%%
# Add relevant data from various MPA-JHU tables
M.add_columns([MBS['D4000_N_SUB'], MSM['MEDIAN'], MEL['OIII_4363_FLUX'],
               MEL['OIII_4959_FLUX'], MEL['OIII_5007_FLUX'], MEL['H_BETA_FLUX'],
               MEL['NII_6548_FLUX'], MEL['NII_6584_FLUX'], MEL['H_ALPHA_FLUX'],
               MSFR['AVG']],
              names=['4000A Break Strength', 'Stellar Mass', 'OIII_4363',
                     'OIII_4959', 'OIII_5007', 'H_beta', 'NII_6548',
                     'NII_6584', 'H_alpha', 'SFR'])

# Produce SkyCoord objects for both MPA-JHU data and FCOZG data
Mco = SkyCoord(M['RA'], M['DEC'], unit='deg')
Lco = SkyCoord(FCOZG['ID_ra'], FCOZG['ID_dec'], unit='deg')

# 2D match of catalogues
idx, d2d, d3d = Lco.match_to_catalog_sky(Mco)
d2d = d2d.arcsecond

# Plot separations between catalogues
# plt.close('all')
plt.figure()
plt.hist(d2d, bins=50)
plt.xlabel('Offset (arcseconds)')
plt.ylabel('Number of Sources')
plt.title('Offset between matched MPA-JHU sources and FCOZG sources')

# Chose to make a cut of within 2 arcseconds to be matched
idxCut = idx[d2d<=2]
d2dCut = d2d[d2d<=2]
# Only keep matched MPA-JHU sources in matched order
Mnew = M[idxCut]
# Define FCOZGM as FCOZG sources with MPA-JHU match within 1arcsec
FCOZGM = FCOZG[d2d<=2]
# Attach relevant MPA-JHU data
FCOZGM = hstack([FCOZGM, Mnew])

# Write to FITS table
# FCOZGM.write('FCOZGM', format='fits')

#%% Figure 3
# Remove extra spaces
FCOZGM['SUBCLASS'] = np.char.strip(FCOZGM['SUBCLASS'])

# Dictate SFG from those classified by MPA-JHU
SFG = FCOZGM[(FCOZGM['SUBCLASS'] == 'STARFORMING') |\
                    (FCOZGM['SUBCLASS'] == 'STARFORMING BROADLINE') |\
                    (FCOZGM['SUBCLASS'] == 'STARBURST') |\
                    (FCOZGM['SUBCLASS'] == 'STARBURST BROADLINE')]
# Find luminosity of SFG
L150SFG = lf.Lum_calc(SFG['Total_flux']*10**(-29), SFG['z_best'])
# Take all those that aren't SFG as classified by MPA-JHU
notSFG = FCOZGM[(FCOZGM['SUBCLASS'] != 'STARFORMING') &\
                    (FCOZGM['SUBCLASS'] != 'STARFORMING BROADLINE') &\
                    (FCOZGM['SUBCLASS'] != 'STARBURST') &\
                    (FCOZGM['SUBCLASS']!= 'STARBURST BROADLINE')]
# Find luminosities of non-SFG
L150notSFG = lf.Lum_calc(notSFG['Total_flux']*10**(-29), notSFG['z_best'])

# plt.close('all')
plt.figure()
plt.scatter(SFG['SFR'][SFG['SFR'] > -99], np.log10(L150SFG)[SFG['SFR'] > -99], s=4, alpha=0.1)
plt.scatter(notSFG['SFR'][notSFG['SFR'] > -99], np.log10(L150notSFG)[notSFG['SFR'] > -99], s=4, alpha=0.1, c='r')
plt.plot([-1.8, 2.9], [20, 27], 'k', linewidth=1)
plt.xlim([-2, 3])
plt.ylim([20, 27])
plt.xticks(np.arange(-2, 4))
plt.yticks(np.arange(20, 28))
plt.xlabel(r'$log_{10}(SFR) (M_{sun}yr^{-1})$')
plt.ylabel(r'$log_{10}(L_{150}) (W Hz^{-1}$)')

#%% Figure 4

# Define function to calculate W2-W3 and W1-W2
def W_calc(Table):
    W2W3 = Table['w2Vega'] - Table['w3Vega']
    W1W2 = Table['w1Vega'] - Table['w2Vega']
    return W2W3, W1W2

# Plot for FCOZGM set and use to dictate FCOZGM AGN and SFG
plt.figure()
plt.scatter(W_calc(FCOZGM)[0], W_calc(FCOZGM)[1], s=4, alpha=0.1)
plt.plot(np.full(len(np.arange(-1, 1.1, 0.1)), 0.8), np.arange(-1, 1.1, 0.1), 'y')
plt.xlim([-2, 3])
plt.ylim([-1, 1])
plt.yticks(np.arange(-1, 1.1, 0.5))
plt.xlabel('W2-W3 (AB mags)')
plt.ylabel('W1-W2 (AB mags)')

# Dictate FCOZGM AGN and SFG based on being either side of the W2-W3 = 0.8 line
FCOZGM_SFG = FCOZGM[W_calc(FCOZGM)[0] > 0.8]
FCOZGM_AGN = FCOZGM[W_calc(FCOZGM)[0] < 0.8]

# Create subset that has WISE data
FCOZGW = FCOZG[FCOZG['AllWISE'] != 'N/A']

# Vertices of polygonal shape dictating SFG region from Hardcastle
xpoints = np.array([5, 3.82, 2.94, 2.73, 3.02, 3.23, 3.89, 5])
ypoints = np.array([0.78, 0.79, 0.65, 0.4, 0.15, 0.04, 0.04, 0.06])
# Find gradients of lines between points
m = (ypoints[np.arange(1, 8)] - ypoints[np.arange(0, 7)])/(xpoints[np.arange(1, 8)] - xpoints[np.arange(0, 7)])
# Find y intercepts of lines
c = ypoints[np.arange(0, 7)] - m*xpoints[np.arange(0, 7)]

# Plot Figure 4 of Hardcastle
plt.figure()
plt.scatter(W_calc(FCOZGW)[0], W_calc(FCOZGW)[1], s=4, alpha=0.1, c='palegreen')
plt.scatter(W_calc(FCOZGM_AGN)[0], W_calc(FCOZGM_AGN)[1], s=4, alpha=0.1, c='darkred')
plt.scatter(W_calc(FCOZGM_SFG)[0], W_calc(FCOZGM_SFG)[1], s=4, alpha=0.1, c='dodgerblue')
plt.plot(xpoints, ypoints, c='k', linewidth=1)
plt.xlim([-0.5, 5])
plt.ylim([-0.5, 2])
plt.xlabel('W2-W3 (Vega mags)')
plt.ylabel('W1-W2 (Vega mags)')
#%%
# Dictate x = W2-W3 and y = W1-W2 from FCOZGW set
x, y = W_calc(FCOZGW)

# Only keep data that lies within the polygonal region (the SFGs)
SFGcut = FCOZGW[(y < m[0]*x + c[0]) & (y < m[1]*x + c[1]) & (y < m[2]*x + c[2]) &\
            (y > m[3]*x + c[3]) & (y > m[4]*x + c[4]) & (y > m[5]*x + c[5]) & (y > m[6]*x + c[6])]

# Plot to show that only data within the polygonal region is kept
plt.figure()
plt.scatter(W_calc(SFGcut)[0], W_calc(SFGcut)[1], s=4, alpha=0.1)
plt.plot(xpoints, ypoints, c='k', linewidth=1)
plt.xlim([-0.5, 5])
plt.ylim([-0.5, 2])
plt.xlabel('W2-W3 (Vega mags)')
plt.ylabel('W1-W2 (Vega mags)')


#%%
# Step 1 in Hardcastle selection method - Only keep data within -33 < K_s < -17
# "Disposes of sources with outlier absolute magnitudes that presumably indicate aberrant redshifts"
Cut1 = FCOZG[(FCOZG['Ks_rest'] > (-33)) & (FCOZG['Ks_rest'] < (-17))]

# Step 2 in Hardcastle selection method - remove sources classed as SFG from FCOZGM   
Cut2 = setdiff(Cut1, FCOZGM_SFG, keys='Source_Name')

# Take data that has no attributed WISE data
noWISE = FCOZG[FCOZG['AllWISE'] == 'N/A']
# Combine with SFGcut data - This is all the data that I want to cut (Step 3) unless the exceptions are met
both = vstack([SFGcut, noWISE])

# Calculate luminosity of both. Observed flux density is in mJy, conversion to W/m^2Hz
L = lf.Lum_calc(both['Total_flux']*10**(-29), both['z_best'])

# Going through Step 3 exceptions
# Keep if sources in both are classified as AGN in FCOZGM
keep1 = np.isin(both['Source_Name'], FCOZGM_AGN['Source_Name'])
# Keep if L > 10^25 WHz^-1 and Ks-band magnitude is > -25 (non-quasars)
keep2 = (L > 10**25) & (both['Ks_rest'] > (-25))
# Keep if Ks-band magnitude is < -25 (quasars) and radio luminosity meets below condition
keep3 = (both['Ks_rest'] < (-25)) & (np.log10(L) > (25.3 - 0.06*(25 + both['Ks_rest'])))

# Combine all exceptions and apply to both - This is data I want to keep
SFGkeep = both[keep1 | keep2 | keep3]
# Remove the data I want to keep from both, so that when I then remove both
# from the total data set, the ones I want to keep are not removed
SFGcut2 = setdiff(both, SFGkeep, keys='Source_Name')
# Remove both (minus the data I want to keep) from the total data set
Cut3 = setdiff(Cut2, SFGcut2, keys='Source_Name') # length = 26,586
RLAGN = Cut3.copy()

# Plot to see that the cuts look about right
plt.figure()
plt.scatter(W_calc(Cut3)[0], W_calc(Cut3)[1], s=4, alpha=0.1)
plt.plot(xpoints, ypoints, c='k', linewidth=1)
plt.xlim([-0.5, 5])
plt.ylim([-0.5, 2])
plt.xlabel('W2-W3 (Vega mags)')
plt.ylabel('W1-W2 (Vega mags)')

# RLAGN.write('RLAGN', format='fits')








