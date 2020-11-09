# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 12:16:13 2020

@author: ppykd1
"""

from astropy.table import Table
# import LOFAR_Functions as lf
import matplotlib.pyplot as plt
import numpy as np


dat = Table.read('Angle data')

cutCond = (dat['d_or'].data >= 1) & (dat['d_or'].data <= 100) & \
    ((dat['Dist_oc'].data).astype(float) >= 0.5) & ((dat['Dist_oc'].data).astype(float) <= 5)

datCut = dat[cutCond]

plt.hist(datCut['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('Angle between Cluster-Optical-Radio sources')
plt.xticks(np.arange(0, 200, 20))
