#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 2019 09:40:50 by Hannah Holt

This programs makes a cheat sheet for switching between pressure lvls and geopotential altitudes

"""

import numpy as np
from scipy import stats

species = 'HE'
alpha = ['-0.38']

infolder = '/Users/hannahholt/Documents/MATLAB/TIEGCM/Contour_textfiles/pressure/' + species + '_alpha_' + alpha[0] + '/'
plvls = np.loadtxt(infolder + "p_lvls.txt", delimiter=',')         # the alpha = 0 has less pressure points.
mode_z0 = np.zeros(len(plvls))

# now read in the global geopotential altitudes for each p and take the average
i=0
for p in plvls:
    p0_infolder = infolder + 'plvl_' + str(p) + '/'
    z0_alts = np.loadtxt(p0_infolder + "zp_altitudes_km_pdrag.txt", delimiter=',')  # the actual altitude values [km] for the p0 lvls
    mode = stats.mode(np.round(z0_alts), axis=None)
    mode_z0[i] = mode[0]
    i=i+1

f = open('./press2alt.txt', 'w')
f.write("%s\t%s\n" % ('P Lev', 'Zpot Global Mode [km]'))
for i in range(len(plvls)):
    f.write("%.3f \t %i\n" % (plvls[i], mode_z0[i]))
f.close()