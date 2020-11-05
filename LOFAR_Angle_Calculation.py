# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:02:52 2020

@author: -
"""

from astropy.table import Table
import LOFAR_Functions as lf
import matplotlib.pyplot as plt

# Import Cluster match data FITS file created in LOFAR_Optical_Cluster_match.py
dat = Table.read('Cluster match data')

# Remove sources that do not have a cluster match
# condition = np.isnan(dat['Cluster RA'].data) == False
condition = dat['Cluster RA'].data != b'None'
dat = dat[condition]

# Take relevant arrays for angle calc
RadRA = dat['Rad RA'].data
RadDec = dat['Rad Dec'].data
OptRA = dat['Opt RA'].data
OptDec = dat['Opt Dec'].data
ClusRA = (dat['Cluster RA'].data).astype(float)
ClusDec = (dat['Cluster Dec'].data).astype(float)


x_oc = [] # x offset between optical and cluster
y_oc = [] # y offset between optical and cluster
x_or = [] # x offset between optical and radio
y_or = [] # y offset between optical and radio

theta_oc = [] # angle between optical and cluster anti-clockwise from North
theta_or = [] # angle between optical and radio anti-clockwise from North
theta_diff = [] # angle between cluster, optical and radio

for i in range(0, len(ClusRA)):
    x_oc.append(lf.x_calc(ClusRA[i], OptRA[i], OptDec[i])) # store all x and y offsets
    y_oc.append(lf.y_calc(ClusDec[i], OptDec[i]))
    x_or.append(lf.x_calc(RadRA[i], OptRA[i], OptDec[i]))
    y_or.append(lf.y_calc(RadDec[i], OptDec[i]))
    
    # Calculate angles by splitting sign of offset values into four quadrants. See project book for detailed diagram.
    if x_oc[i] > 0 and y_oc[i] > 0:
        theta_oc.append(lf.vertAngle(x_oc[i], y_oc[i]))
    elif x_oc[i] > 0 and y_oc[i] < 0:
        theta_oc.append(180 - lf.vertAngle(x_oc[i], abs(y_oc[i])))
    elif x_oc[i] < 0 and y_oc[i] < 0:
        theta_oc.append(180 + lf.vertAngle(x_oc[i], y_oc[i]))
    elif x_oc[i] < 0 and y_oc[i] > 0:
        theta_oc.append(360 - lf.vertAngle(abs(x_oc[i]), y_oc[i]))
        
    if x_or[i] > 0 and y_or[i] > 0:
        theta_or.append(lf.vertAngle(x_or[i], y_or[i]))
    elif x_or[i] > 0 and y_or[i] < 0:
        theta_or.append(180 - lf.vertAngle(x_or[i], abs(y_or[i])))
    elif x_or[i] < 0 and y_or[i] < 0:
        theta_or.append(180 + lf.vertAngle(x_or[i], y_or[i]))
    elif x_or[i] < 0 and y_or[i] > 0:
        theta_or.append(360 - lf.vertAngle(abs(x_or[i]), y_or[i]))
    
    diff = abs(theta_or[i] - theta_oc[i]) # take difference between angles
    if diff > 180: # only angle < 180 is wanted
        theta_diff.append(360 - diff)
    else:
        theta_diff.append(diff)
      

plt.close('all')

plt.hist(theta_diff, bins=20)
plt.xlabel('Angle between Cluster-Optical-Radio (deg)')
plt.ylabel('Number of sources')

dist_oc = (dat['Dist_oc'].data).astype(float)
plt.figure()
plt.hist(dist_oc, bins=20)
plt.xlabel('Distance between cluster centre and optical source (Mpc)')
plt.ylabel('Number of sources')


# Adding to the FITS table to include the angle between ROC
dat.add_column(theta_diff, name='Angle ROC')
dat['Angle ROC'].unit = 'deg'

# dat.write('Angle data', format = 'fits')
    
    
    
    
    
    