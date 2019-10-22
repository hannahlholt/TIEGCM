#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
"""
Created by Hannah L Holt 1/8/2019


This program uses output of the TIEGCM_Diff_Contour_plots_V4.m to find the altitudes that thermal diffusion alone causes a 10%
change in a species density at 400 km. THis is done by subtracting the percent differences of model vs. pure diffusive eq.
files ((alpha = -0.38) - (alpha = 0)) for every lower boundary integration altitude z0 and seeing where this new lat, long contour map causes a 10% change
Subtracting the files with  no thermal diffusion + winds from the ones w/ thermal diffusion + winds isolates the effects of thermal
diffusion.
"""


import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage

species = 'HE'                  # which species
alpha = ['-0.38', '0']          # thermal diffusion factors


bot_lim = 10             # percent difference limits (we want to be around 10%)
#top_lim = 12

if species == 'HE' and alpha == '0':
    name = '3.125'
elif species == 'HE' and alpha == '-0.38':
    name = '3.375'
else:
    name = '3.125'


# --------------------------------------------------------------------------------------------------------------------

# Load magnetic equator data points to be overlayed on plot
mag_equator = np.loadtxt("Magnetic_equator_lat_lon.txt", delimiter=',')
marklon = 0
magnetic_x = np.linspace(0, 143 / 6, 360)

LT_ticks = np.array([12., 15., 18., 21., 0., 3., 6., 9., 12.])
alt_10_perc = np.zeros((72, 144))       # altitudes at which lower boundary z0 stuff caused a 10% change away from diffusive equilibrium


x = np.linspace(0, 143 / 6, 144)    # longitude values [deg]
y = np.linspace(-88.75, 88.75, 72)  # latitude values [deg]
X, Y = np.meshgrid(x, y)

infolder1 = '/Users/hannahholt/Documents/MATLAB/TIEGCM/Contour_textfiles/pressure/' + species + '_alpha_' + alpha[0] + '/'
infolder2 = '/Users/hannahholt/Documents/MATLAB/TIEGCM/Contour_textfiles/pressure/' + species + '_alpha_' + alpha[1] + '/'

# ---------------------------------------------------------------------------------------------------------------------

# first load the number of lower boundary altitude points
print('Loading pressure altitude points.\n')
# first load the of lower boundary pressure level points
plvls = np.loadtxt(infolder1 + "p_lvls.txt", delimiter=',')

# reverse ordered going from 400km down to z0
plvls = plvls[::-1]

for p0 in plvls:
    p0_infolder1 = infolder1 + 'plvl_' + str(p0) + '/'  # input folders for text files from every z0 point for thermal diff.
    p0_infolder2 = infolder2 + 'plvl_' + str(p0) + '/'  # input folders for text files from every z0 point for NO therm diff

    # first, load the all the data
    data_therm = np.loadtxt((p0_infolder1 + species + "_dens_perc_diff_400_pdrag.txt"), delimiter=',')
    data_notherm = np.loadtxt((p0_infolder2 + species + "_dens_perc_diff_400_pdrag.txt"), delimiter=',')

    data_sub = np.subtract(data_therm, data_notherm)
    z0_alts = np.loadtxt(p0_infolder1 + "zp_altitudes_km_pdrag.txt", delimiter=',')  # the actual altitude values [km] for the p0 lvls

    # find the altitudes for this specific z0 map that caused an ~10% change after subtracting the no therm data from therm data
    for i in range(0, 72):
        for j in range(0, 144):
            if (abs(data_sub[i, j]) > bot_lim) and (alt_10_perc[i, j] == 0):
                alt_10_perc[i,j] = z0_alts[i,j]


# smooth function with 1 sigma gaussian filter
alt_smooth = ndimage.gaussian_filter(alt_10_perc, sigma=1, order=0)

# -----  Percent differences -----------------------
plt.figure()
cont = plt.contourf(X, Y, alt_smooth, 10, colors='k')
levels = np.linspace(110, np.amax(alt_smooth))
myplot = plt.contourf(X, Y, alt_smooth, levels, cmap='jet')
cbar = plt.colorbar(myplot, format='%.0f')
cbar.ax.set_ylabel('Lower Boundary Altitude [km]')
plt.title(str(bot_lim) + '% ' + species + ' Density Change at ~400km due to Thermal Diffusion')

plt.xticks(np.arange(0., 24., 3.), LT_ticks)
plt.plot(magnetic_x, mag_equator[:, 1], 'r')  # Add magnetic equator line
plt.xlabel('Local Solar Time [hr]')
plt.ylabel('Latitude [deg]')
plt.savefig('./10percChange_Figs/' + species + '_alt10perc_thermDiff_ONLY_V2.pdf')
plt.show()

