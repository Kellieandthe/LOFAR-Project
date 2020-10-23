# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 15:06:53 2020

@author: kelli
"""

from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import math
import LOFAR_Functions as lf


LOFAR = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/HETDEX associations and optical IDs.fits', format='fits')

RadSourceFull = LOFAR['Source_Name']
OptSourceFull = LOFAR['ID_name']
RadRAFull = LOFAR['RA']
RadDecFull = LOFAR['DEC']
OptRAFull = LOFAR['ID_ra']
OptDecFull = LOFAR['ID_dec']
RadRAErrFull = LOFAR['E_RA']
RadDecErrFull = LOFAR['E_DEC']
zFull = LOFAR['z_best']
zSourceFull = LOFAR['z_best_source']

RadSource = []
OptSource = []
RadRA = []
RadDec = []
OptRA = []
OptDec = []
RadRAErr = []
RadDecErr = []
Opt_z = []
zSource = []

for i in range(0, len(OptRAFull)):
    if np.isnan(zFull[i]) == False and np.isnan(OptRAFull[i]) == False:
        RadRA.append(RadRAFull[i])
        RadDec.append(RadDecFull[i])
        OptRA.append(OptRAFull[i])
        OptDec.append(OptDecFull[i])
        RadSource.append(RadSourceFull[i])
        OptSource.append(OptSourceFull[i])
        RadRAErr.append(RadRAErrFull[i])
        RadDecErr.append(RadDecErrFull[i])
        Opt_z.append(zFull[i])
        if zSourceFull[i] == 1.0:
            zSource.append('Spectroscopic')
        else:
            zSource.append('Photometric')
          
d_or = []
x_err = []
y_err = RadDecErr
d_err = []

for i in range(0, len(OptRA)):
    d_or.append(lf.offset(RadRA[i], OptRA[i], RadDec[i], OptDec[i]))
    x = lf.x_calc(RadRA[i], OptRA[i], OptDec[i])
    y = lf.y_calc(RadDec[i], OptDec[i])
    x_err.append(math.cos(math.radians(OptDec[i])) * RadRAErr[i])
    d_err.append(math.sqrt((x**2 * x_err[i]**2 + y**2 * y_err[i]**2) / 
                           (x**2 + y**2)))
    
    
AllRadData = Table([RadSource, RadRA, RadDec, RadRAErr, RadDecErr, OptSource, 
                    OptRA, OptDec, Opt_z, zSource, d_or, d_err],
                   names = ('Radio Source ID', 'Rad RA', 'Rad Dec', 'RA Error',
                            'Dec Error', 'Optical Source ID', 'Opt RA',
                            'Opt Dec', 'z', 'z Source', 'd_or', 'd_or Error'))

AllRadData['Rad RA'].unit = 'deg'
AllRadData['Rad Dec'].unit = 'deg'
AllRadData['RA Error'].unit = 'arcsec'
AllRadData['Dec Error'].unit = 'arcsec'
AllRadData['Opt RA'].unit = 'deg'
AllRadData['Opt Dec'].unit = 'deg'
AllRadData['d_or'].unit = 'arcsec'
AllRadData['d_or Error'].unit = 'arcsec'


AllRadData.write('Radio source data', format='fits')
    
plt.hist(d_or, bins=50, log=True)
plt.xlabel('Offset (arcseconds)')
plt.ylabel('Number of sources')
plt.title('Offset from radio source to optical source')

def plot_loghist(x, bins):
    hist, bins = np.histogram(x, bins=bins)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.hist(x, bins=logbins, log=True)
    plt.xscale('log')

plt.figure()    
plot_loghist(d_or, 50)
plt.xlabel('Offset (arcseconds)')
plt.ylabel('Number of sources')
plt.title('Offset from radio source to optical source')


RadIDColTen = []
OptIDColTen = []
offsetColTen = []

RadIDColHun = []
OptIDColHun = []
offsetColHun = []

for i in range(0, len(d_or)):
    if d_or[i] > 10:
        RadIDColTen.append(RadSource[d_or.index(d_or[i])])
        OptIDColTen.append(OptSource[d_or.index(d_or[i])])
        offsetColTen.append(d_or[i])
    if d_or[i] > 100:
        RadIDColHun.append(RadSource[d_or.index(d_or[i])])
        OptIDColHun.append(OptSource[d_or.index(d_or[i])])
        offsetColHun.append(d_or[i])
        
OverTen = Table([RadIDColTen, OptIDColTen, offsetColTen], names = ('Radio Source ID', 'Optical Source ID', 'offset (arcseconds)'))
OverHun = Table([RadIDColHun, OptIDColHun, offsetColHun], names = ('Radio Source ID', 'Optical Source ID', 'offset (arcseconds)'))





    

    