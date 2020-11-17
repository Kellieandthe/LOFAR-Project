# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 12:16:13 2020

@author: ppykd1
"""

from astropy.table import Table
# import LOFAR_Functions as lf
import matplotlib.pyplot as plt
import numpy as np


RMdat = Table.read('RM Angle data')
WHLdat = Table.read('WHL Angle data')

RMcutCond = (RMdat['Radio-Optical Offset'].data >= 1) & (RMdat['Radio-Optical Offset'].data < 100) & (RMdat['Optical z'].data > 0.05)
WHLcutCond = (WHLdat['Radio-Optical Offset'].data >= 1) & (WHLdat['Radio-Optical Offset'].data < 100) & (WHLdat['Optical z'].data > 0.05)

RMdatCut = RMdat[RMcutCond]
WHLdatCut = WHLdat[WHLcutCond]

plt.close('all')

plt.subplot(2, 1, 1)
plt.hist(WHLdatCut['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('WHL15')
plt.xticks(np.arange(0, 200, 20))

plt.subplot(2, 1, 2)
plt.hist(RMdatCut['Angle ROC'], bins=18)
plt.xlabel('Angle (deg)')
plt.ylabel('Number of sources')
plt.title('redMaPPer')
plt.xticks(np.arange(0, 200, 20))
plt.suptitle('Angle between Cluster-Optical-Radio sources')
plt.tight_layout()

