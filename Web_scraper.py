# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 09:51:21 2021

@author: ppykd1
"""

from astropy.table import Table, join
import LOFAR_Functions as lf
import numpy as np
import requests

# Import data
WHL = Table.read('Data\WHL angle data', format = 'fits')
# Only keep NATs
NAT = lf.NAT_cut(WHL)

# Import Mingo data
Mingo = Table.read('Data\Mingo NATs.fits', format='fits')
Mingo.rename_column('source', 'Radio Source')
# Join Mingo data with NAT data
NAT = join(NAT, Mingo, keys='Radio Source', join_type='left')

# Rename mosaic IDs to replace + with _ for URL link
for j in np.arange(0, len(NAT)):
    if ('+' in NAT['Mosaic_ID'][j]) == True:
        NAT['Mosaic_ID'][j] = NAT['Mosaic_ID'][j].replace('+', '_')
        

for i in np.arange(0, len(NAT)):
    # Change Mosaic ID, RA and DEC for URL of each source, but keep size of image as sdec=0.03 and sra = 0.05
    image_url = 'https://vo.astron.nl/getproduct/hetdex/data/mosaics/'+str(NAT['Mosaic_ID'][i])+'-mosaic.fits?sdec=0.03&dec='+str(NAT['Radio DEC'][i])+'&ra='+str(NAT['Radio RA'][i])+'&sra=0.05'
    # Request URL
    r = requests.get(image_url, stream=True)
    # Dictate .fits file name and folder location
    image_name = 'Images/NAT'+str(NAT['index'][i])+'.fits'
    # Open and write to file name
    with open(image_name, 'wb') as f:
        f.write(r.content)
