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

# Import FITS table of LoTTS radio & optical data
LOFAR = Table.read('C:/Users/ppykd1/Documents/PhD/LOFAR Project/Data/HETDEX associations and optical IDs.fits', format='fits')

# Remove sources that have no associated redshift value
condition = np.isnan((LOFAR['z_best'].compressed()).data) == False
LOFAR = LOFAR[condition]

# Lots of data in the imported table is unneeded, so separate the table into helpful arrays
# The arrays import as 'masked' so the compressed function fixes that
RadSource = (LOFAR['Source_Name'].compressed()).data
OptSource = (LOFAR['ID_name'].compressed()).data
RadRA = (LOFAR['RA'].compressed()).data
RadDec = (LOFAR['DEC'].compressed()).data
OptRA = (LOFAR['ID_ra'].compressed()).data
OptDec = (LOFAR['ID_dec'].compressed()).data
RadRAErr = (LOFAR['E_RA'].compressed()).data
RadDecErr = (LOFAR['E_DEC'].compressed()).data
Opt_z = (LOFAR['z_best'].compressed()).data
zSource = (LOFAR['z_best_source'].compressed()).data


d_or = [] # offset between optical and radio source
x_err = [] # error in x
y_err = RadDecErr # error in y
d_err = [] # error in offset

# Calculate offset, x, y, and errors
for i in range(0, len(OptRA)):
    d_or.append(lf.offset(RadRA[i], OptRA[i], RadDec[i], OptDec[i]))
    x = lf.x_calc(RadRA[i], OptRA[i], OptDec[i])
    y = lf.y_calc(RadDec[i], OptDec[i])
    x_err.append(math.cos(math.radians(OptDec[i])) * RadRAErr[i])
    d_err.append(math.sqrt((x**2 * x_err[i]**2 + y**2 * y_err[i]**2) / 
                            (x**2 + y**2)))

# Convert to arrays from lists    
d_or = np.array(d_or)
x_err = np.array(x_err)
d_err = np.array(d_err)

# Put all data in a table to be exported as FITS file
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


# AllRadData.write('Radio source data', format='fits')

plt.close('all')    

# Plot histogram with logarithmic y axis
plt.hist(d_or, bins=50, log=True)
plt.xlabel('Offset (arcseconds)')
plt.ylabel('Number of sources')
plt.title('Offset from radio source to optical source')

# Plot histogram with logarithmic x and y axes
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

#%%

RadIDColTen = []
OptIDColTen = []
offsetColTen = []

RadIDColHun = []
OptIDColHun = []
offsetColHun = []

# Determine potential outliers in offsets if >10 or >100 arcseconds
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





    

