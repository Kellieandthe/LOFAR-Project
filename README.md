# Kellie's LOFAR Project

Kellie de Vos

*University of Nottingham*

## Background

As radio galaxies orbit around their host galaxy cluster centre, the radio emissions they release are sometimes slowed down by the surrounding intracluster medium. This creates an offset in position between the observed radio source and its matching optical source. The angle between the cluster centre, optical source and radio source can then be calculated to determine the orbital direction of the galaxy around the galaxy cluster centre.

This has already been attempted by Garon et al. (2019) with radio sources with wide-angle tails (WATs), and looked at by Croston et al. (2018).

## Overview

I am using data from the LOFAR Two-Metre Sky Survey (LoTSS) Public Data Release 1, and cluster catalogues by Wen et al. (2015) and Rykoff et al. (2014), both of which utilise the Sloan Digital Sky Survey (SDSS).

I first analyse the LoTSS data by determining the offset between each radio source and its corresponding optical source. This is done taking the RA and Dec of the radio and optical sources and using an offset formula to determine the separation between them in arcseconds. This has been done in the module LOFAR_Radio_Optical_match.py.

I found that the most efficient way to be able to call the functions I have made is to put them in their own module, so they have been put into a module called LOFAR_Functions.py.

Then, since I am using two cluster catalogues, I cross-reference them to check for any duplicates. This has been done in module LOFAR_Duplicate_Clusters.py.

I then hope to match each optical source to its potential galaxy cluster by using parameters from Garon et al. I determine if the optical source is within |Δz| < 0.04 (where Δz is defined in Garon et al.) of the cluster centre. If it is, I determine whether it is within a 15Mpc radius of the cluster centre. The matching cluster for each optical source is then taken as the closest cluster centre within that 15Mpc radius. This is being attempted in module LOFAR_Optical_Cluster_match.py.

Once each optical source has been assigned a host galaxy cluster, I then take the angle between the cluster centre, the the optical source, and the radio source using the cosine rule. The lengths in the cosine rule are given by the offsets between pairs of each of the three points. This is going to be attempted in module LOFAR_Angle_Calculation.py.

Throughout my project so far I have used modules such as astropy.table, numpy, math, and matplotlib.pyplot.

## Example code

```python

```
