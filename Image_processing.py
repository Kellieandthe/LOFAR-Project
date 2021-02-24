# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 09:56:55 2021

@author: ppykd1
"""

from astropy.io import fits
from astropy.table import Table, join
from astropy.wcs import WCS
import LOFAR_Functions as lf
import matplotlib.pyplot as plt
import matplotlib.image as img
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

WHL = Table.read('Data\WHL angle data', format = 'fits')
NAT = lf.NAT_cut(WHL)

Mingo = Table.read('Data\Mingo NATs.fits', format='fits')
Mingo.rename_column('source', 'Radio Source')
NAT = join(NAT, Mingo, keys='Radio Source', join_type='left')

np.putmask(NAT['PA1'], NAT['PA1'] < 0, (360 + NAT['PA1']))
np.putmask(NAT['PA2'], NAT['PA2'] < 0, (360 + NAT['PA2']))

coords = np.transpose(np.array([NAT['Radio RA'], NAT['Radio DEC']]))
gal = SkyCoord(NAT['Optical RA'], NAT['Optical DEC'], unit='deg')
rad = SkyCoord(NAT['Radio RA'], NAT['Radio DEC'], unit='deg')

theta_gr = gal.position_angle(rad).degree

angDiff = abs(NAT['PA2'] - NAT['PA1'])
np.putmask(angDiff, angDiff > 180, (angDiff - 360))
perpBi = np.array([])
for i in np.arange(0, len(NAT)):
    perpBi = np.append(perpBi, min([NAT['PA1'][i], NAT['PA2'][i]]) + angDiff[i]/2)
np.putmask(perpBi, perpBi < 0, (360 + perpBi))

# Converting everything from 0<=theta<=360 to -180<=theta<=180
# np.putmask(theta_gr, theta_gr > 180, (theta_gr - 360))
# np.putmask(NAT['PA1'], NAT['PA1'] > 180, (NAT['PA1'] - 360))
# np.putmask(NAT['PA2'], NAT['PA2'] > 180, (NAT['PA2'] - 360))
# np.putmask(perpBi, perpBi > 180, (perpBi - 360))

plt.close('all')
plt.figure(figsize=(17,5))
plt.subplot(1, 3, 1)
plt.hist(theta_gr)
plt.xlabel(r'$\theta_r$')
plt.ylabel('Number')
# plt.xticks(np.arange(-180, 210, 30))
plt.xticks(np.arange(0, 390, 30))
plt.title('Angle to radio source from galaxy, \nmeasured anti-CW from North')
plt.subplot(1, 3, 2)
plt.hist(perpBi)
plt.xlabel(r'$\theta_{bi}$')
plt.ylabel('Number')
# plt.xticks(np.arange(-180, 210, 30))
plt.xticks(np.arange(0, 390, 30))
plt.title('Angle to bisector of brightest peaks from galaxy, \nmeasured anti-CW from North')
plt.subplot(1, 3, 3)
plt.hist(NAT['PA1'])
plt.xlabel(r'$\theta_{bp}$')
plt.ylabel('Number')
# plt.xticks(np.arange(-180, 210, 30))
plt.xticks(np.arange(0, 390, 30))
plt.title('Angle to brightest peak from galaxy, \nmeasured anti-CW from North')

Diff1 = abs(theta_gr - perpBi)
Diff2 = abs(theta_gr - NAT['PA1'])
np.putmask(Diff1, Diff1 > 180, (360 - Diff1))
np.putmask(Diff2, Diff2 > 180, (360 - Diff2))
plt.figure(figsize=(16,6))
plt.subplot(1, 2, 1)
plt.hist(Diff1)
plt.xlabel(r'$\Delta\theta$')
plt.ylabel('Number')
plt.xticks(np.arange(0, 210, 30))
plt.title(r'Difference between $\theta_r$ and $\theta_{bi}$, folded around $180^{\circ}$')
plt.subplot(1, 2, 2)
plt.hist(Diff2)
plt.xlabel(r'$\Delta\theta$')
plt.ylabel('Number')
plt.xticks(np.arange(0, 210, 30))
plt.title(r'Difference between $\theta_r$ and $\theta_{bp}$, folded around $180^{\circ}$')
    
# def brightCoord(theta, dist):
#     if theta == 0:
#         xbp = galPix_x
#         ybp = galPix_y + dist
#     elif theta > 0 and theta < 90:
#         xbp = galPix_x - dist*np.sin(theta)
#         ybp = galPix_y + dist*np.cos(theta)
#     elif theta == 90:
#         xbp = galPix_x - dist
#         ybp = galPix_y
#     elif theta > 90 and theta < 180:
#         theta = theta - 90
#         xbp = galPix_x - dist*np.cos(theta)
#         ybp = galPix_y - dist*np.sin(theta)
#     elif theta == 180:
#         xbp = galPix_x
#         ybp = galPix_y - dist
#     elif theta > 180 and theta < 270:
#         theta = theta - 180
#         xbp = galPix_x + dist*np.cos(theta)
#         ybp = galPix_y - dist*np.sin(theta)
#     elif theta == 270:
#         xbp = galPix_x + dist
#         ybp = galPix_y
#     elif theta > 270 and theta < 360:
#         theta = theta - 270
#         xbp = galPix_x + dist*np.cos(theta)
#         ybp = galPix_y + dist*np.sin(theta)
#     return xbp, ybp
#%%
plt.close('all')
n = 10
for i in np.arange(0, n):
    x = np.array([NAT['Radio RA'][i], NAT['Optical RA'][i], NAT['Cluster RA'][n]]) * u.deg
    y = np.array([NAT['Radio DEC'][i], NAT['Optical DEC'][i], NAT['Cluster DEC'][n]]) * u.deg

    # Image is of field size 0.03deg
    name = 'Images/NAT'+str(i+1)+'.fits'
    NATim = fits.getdata(name)
    HDU = fits.getheader(name)

    wmap = WCS(HDU)
    
    # galPix_x, galPix_y = wmap.world_to_pixel(gal[i])
    # BP1_x, BP1_y = brightCoord(NAT['PA1'][i], NAT['dist1'][i])
    # BP2_x, BP2_y = brightCoord(NAT['PA2'][i], NAT['dist2'][i])

    fig = plt.figure(figsize=(13,9))
    ax = plt.subplot(1, 2, 1, projection=wmap)
    ax.imshow(NATim, origin='lower')
    plt.xlabel('RA')
    plt.ylabel('DEC')

    ra = ax.coords[0]
    dec = ax.coords[1]
    ra.set_major_formatter('d.ddd')
    dec.set_major_formatter('d.ddd')

    ax.plot(x, y, transform=ax.get_transform('fk5'), color='lavender', linestyle='--')
    ax.scatter(x[0], y[0], transform=ax.get_transform('fk5'), color='coral', s=40, marker='s')
    ax.scatter(x[1], y[1], transform=ax.get_transform('fk5'), color='purple', s=40, marker='x')
    plt.legend(['Angle with cluster', 'Radio source', 'Galaxy'])
    plt.xlim([10, NATim.shape[1]-10])
    plt.ylim([10, NATim.shape[0]-10])
    plt.title(u'Radio source angle anti-cw from North = %.1f\N{DEGREE SIGN}' % theta_gr[i])
    
    Mingo_im = img.imread('Images/Mingo images/'+str(NAT['index'][i])+'.png')
    plt.subplot(1, 2, 2)
    plt.imshow(Mingo_im)
    plt.axis('off')
    plt.title(u'Perp. bisector angle anti-cw from North = %.1f\N{DEGREE SIGN}'  % perpBi[i] +
              u'\nBrightest source angle anti-cw from North = %.1f\N{DEGREE SIGN}' % NAT['PA1'][i])
    
    

