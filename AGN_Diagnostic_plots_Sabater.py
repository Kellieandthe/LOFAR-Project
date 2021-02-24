# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 10:58:21 2020

@author: ppykd1
"""

from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
# from astropy import units as u
# from astropy.coordinates import SkyCoord
import LOFAR_Functions as lf
from astropy.cosmology import FlatLambdaCDM

# Set cosmology according to Garon et al. paper
cosmo = FlatLambdaCDM(H0=68, Om0=0.315)

FCOZGM = Table.read('Data\FCOZGM', format='fits')
#%%
# Observed flux density in mJy, conversion to W/m^2Hz
L150 = lf.Lum_calc(FCOZGM['Total_flux']*10**(-29), FCOZGM['z_best'])
L150SM = L150[FCOZGM['Stellar Mass'] > 0]
datSM = FCOZGM[FCOZGM['Stellar Mass'] > 0]
M = 10**datSM['Stellar Mass']

line5 = 1.45 - 0.55*(np.arange(10, 12.2, 0.1) - 12.2)

xvals = np.array([12.846153846153847, 13.298076923076923, 13.875,
                 14.326923076923077, 14.807692307692308, 15.086538461538462,
                 15.317307692307693, 15.625, 15.923076923076923])

yvals = np.array([1.4552845528455287, 1.3739837398373986, 1.3089430894308944,
                 1.280487804878049, 1.2601626016260166, 1.247967479674797,
                 1.2317073170731712, 1.2276422764227646, 1.2113821138211385])

z = np.polyfit(xvals, yvals, 3)
xdat = np.arange(12.2, 16, 0.1)
p = np.poly1d(z)
line6 = p(xdat)

plt.close('all')
plt.figure()
plt.scatter(np.log10(L150SM/M), datSM['4000A Break Strength'],  s=4)
plt.plot(np.arange(10, 12.2, 0.1), line5, 'y')
# plt.plot(xdat, line6, 'y')
plt.xlim([10, 16])
plt.ylim([0.5, 2.5])
plt.xlabel(r'$log_{10}(L_{150}/M^*)/W Hz^{-1} M_{sun}^{-1}$')
plt.ylabel(r'$D_n(4000)$')

#%%

x = np.log10(FCOZGM['NII_6584']/FCOZGM['H_alpha'])
y = np.log10(FCOZGM['OIII_5007']/FCOZGM['H_beta'])
line = 1.3 + 0.61/(np.arange(-1.5, 0.05, 0.02) - 0.05)

#plt.close('all')
plt.figure()
plt.scatter(x, y, s=4)
plt.plot(np.arange(-1.5, 0.05, 0.02), line, 'y')
plt.xlim([-1.5, 1])
plt.ylim([-1, 1.5])
plt.xlabel(r'$log_{10}$([NII]/H-alpha)')
plt.ylabel(r'$log_{10}$([OIII]/H-beta)')

#%%

# H_alpha flux given in 10^-17 erg/s/cm^2
# = 10^-17 ergs^-1cm^-2   erg = 10^-7 J
# = 10^-24 Js^-1cm^-2     J/s = W
# = 10^-24 Wcm^-2         cm^-2 = 10^4 m^-2
# = 10^-20 Wm^-2
# Lum_calc function takes flux in units of Wm^-2Hz^-1, so just carry
# through the Hz unit and it cancels to leave final L unit of W)
LHalpha = lf.Lum_calc(FCOZGM['H_alpha']*10**-20, FCOZGM['Z'])/(3.828*10**26)
LHalpha = np.array(LHalpha).astype(float)
cond = (LHalpha > 0) & (L150 > 0)

line2 = np.arange(20, 28, 0.1) - 16.9
line3 = np.arange(20, 28, 0.1) - 16.1

#plt.close('all')
plt.figure()
plt.scatter(np.log10(L150[cond]), np.log10(LHalpha[cond]), s=4)
plt.plot(np.arange(20, 28, 0.1), line2, 'y')
plt.plot(np.arange(20, 28, 0.1), line3, 'y')
plt.xlim([20, 26.5])
plt.ylim([4, 9.5])
plt.xlabel(r'$log_{10}(L_{150}/W Hz^{-1})$')
plt.ylabel(r'$log_{10}(L_{Ha}$/solar units)')

#%%

W2W3 = FCOZGM['w2Mag'] - FCOZGM['w3Mag']
W1W2 = FCOZGM['w1Mag'] - FCOZGM['w2Mag']

#plt.close('all')
plt.figure()
plt.scatter(W2W3, W1W2, s=4)
plt.plot(np.full(len(np.arange(-1, 1.1, 0.1)), 0.8), np.arange(-1, 1.1, 0.1), 'y')
plt.xlim([-2, 3])
plt.ylim([-1, 1])
plt.yticks(np.arange(-1, 1.1, 0.5))
plt.xlabel('W2-W3 (AB mags)')
plt.ylabel('W1-W2 (AB mags)')

FCOZGM_SFG = FCOZGM[W2W3 > 0.8]
FCOZGM_AGN = FCOZGM[W2W3 < 0.8]









