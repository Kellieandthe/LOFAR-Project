# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:02:52 2020

@author: -
"""

from astropy.table import Table
import LOFAR_Functions as lf
import matplotlib.pyplot as plt


data = Table.read('Cluster match data')

RadRAFull = data['Rad RA']
RadDecFull = data['Rad Dec']
OptRAFull = data['Opt RA']
OptDecFull = data['Opt Dec']
ClusRAFull = data['Cluster RA']
ClusDecFull = data['Cluster Dec']

dist_ocFull = data['Dist_oc']
d_orFull = data['d_or']

RadRA = []
RadDec = []
OptRA = []
OptDec = []
ClusRA = []
ClusDec = []
dist_oc = []
d_or = []


for i in range(0, len(ClusRAFull)):
    if ClusRAFull[i] != 'None':
        RadRA.append(RadRAFull[i])
        RadDec.append(RadDecFull[i])
        OptRA.append(OptRAFull[i])
        OptDec.append(OptDecFull[i])
        ClusRA.append(float(ClusRAFull[i]))
        ClusDec.append(float(ClusDecFull[i]))
        d_or.append(d_orFull[i])
        dist_oc.append(float(dist_ocFull[i]))
        

d_oc = []
d_rc = []
x_oc = []
y_oc = []
x_or = []
y_or = []

theta_oc = []
theta_or = []
theta_diff = []
    
for i in range(0, len(ClusRA)):
    d_oc.append(lf.offset(ClusRA[i], OptRA[i], ClusDec[i], OptDec[i]))
    d_rc.append(lf.offset(RadRA[i], ClusRA[i], RadDec[i], ClusDec[i]))
    # theta.append(lf.angleCalc(d_oc[i], d_or[i], d_rc[i]))
    x_oc.append(lf.x_calc(ClusRA[i], OptRA[i], OptDec[i]))
    y_oc.append(lf.y_calc(ClusDec[i], OptDec[i]))
    x_or.append(lf.x_calc(RadRA[i], OptRA[i], OptDec[i]))
    y_or.append(lf.y_calc(RadDec[i], OptDec[i]))
    
    if x_oc[i] > 0 and y_oc[i] > 0:
        theta_oc.append(lf.vertAngle(x_oc[i], y_oc[i]))
    elif x_oc[i] > 0 and y_oc[i] < 0:
        theta_oc.append(90 + lf.horAngle(x_oc[i], abs(y_oc[i])))
    elif x_oc[i] < 0 and y_oc[i] < 0:
        theta_oc.append(180 + lf.vertAngle(x_oc[i], y_oc[i]))
    elif x_oc[i] < 0 and y_oc[i] > 0:
        theta_oc.append(270 + lf.horAngle(abs(x_oc[i]), y_oc[i]))
        
    if x_or[i] > 0 and y_or[i] > 0:
        theta_or.append(lf.vertAngle(x_or[i], y_or[i]))
    elif x_or[i] > 0 and y_or[i] < 0:
        theta_or.append(90 + lf.horAngle(x_or[i], abs(y_or[i])))
    elif x_or[i] < 0 and y_or[i] < 0:
        theta_or.append(180 + lf.vertAngle(x_or[i], y_or[i]))
    elif x_or[i] < 0 and y_or[i] > 0:
        theta_or.append(270 + lf.horAngle(abs(x_or[i]), y_or[i]))
    
    diff = abs(theta_or[i] - theta_oc[i])
    if diff > 180:
        theta_diff.append(360 - diff)
    else:
        theta_diff.append(diff)
    
plt.hist(theta_diff, bins=20)
plt.xlabel('Angle between Cluster-Optical-Radio (deg)')
plt.ylabel('Number')

plt.figure()
plt.hist(dist_oc, bins=20)
plt.xlabel('Distance between cluster centre and optical source (Mpc)')
plt.ylabel('Number')

# AngleData = Table([RadRA, RadDec, OptRA, OptDec, ClusRA, ClusDec, x_oc, y_oc,
#                     x_or, y_or, theta_oc, theta_or, theta_diff, d_or, d_oc, d_rc],
#                   names = ('Radio RA', 'Radio Dec', 'Optical RA', 'Optical Dec',
#                             'Cluster RA', 'Cluster Dec', 'x_oc', 'y_oc', 'x_or',
#                             'y_or', 'theta_oc', 'theta_or', 'theta_diff',
#                             'd_or', 'd_oc', 'd_rc'))

# AngleData.write('Angle data', format = 'fits')
    
    
    
    
    
    