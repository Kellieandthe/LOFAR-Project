# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:02:52 2020

@author: -
"""

from astropy.table import Table
import LOFAR_Functions as lf


data = Table.read('Cluster match data')

RadRA = ['Rad RA']
RadDec = ['Rad Dec']
OptRA = ['Opt RA']
OptDec = ['Opt Dec']
ClusRA = ['Cluster RA']
ClusDec = ['Cluster Dec']

d_or = data['d_or']
d_oc = []
d_rc = []
theta = []

for i in range(0, len(ClusRA)):
    d_oc.append(lf.offset(ClusRA[i], OptRA[i], ClusDec[i], OptDec[i]))
    d_rc.append(lf.offset(RadRA[i], ClusRA[i], RadDec[i], ClusDec[i]))
    theta.append(lf.angleCalc(d_oc[i], d_or[i], d_rc[i]))
    
    
    
    
    