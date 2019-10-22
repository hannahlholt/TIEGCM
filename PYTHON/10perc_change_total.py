#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
"""
Created by Hannah L Holt 1/4/2019


This program uses output of the TIEGCM_Diff_Contour_plots.m to find a global lower integration altitude value where model
vs. pure diffusive equilibrium values vary by 10% at 400 km

"""




import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage


species = 'N2'          # which species
alpha = '0'             # thermal diffusion factor

if species == 'HE' and alpha == '0':
    name = '3.375'
elif species == 'HE' and alpha == '-0.38':
    name = '3.375'
else:
    name = '3.375'



bot_lim = 10             # percent difference limits (we want to be around 10%)
#top_lim = 12


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

infolder = '/Users/hannahholt/Documents/MATLAB/TIEGCM/Contour_textfiles/pressure/' + species + '_alpha_' + alpha + '/'

# ---------------------------------------------------------------------------------------------------------------------

print('Loading ' + infolder + '\n')

# first load the of lower boundary pressure level points
plvls = np.loadtxt(infolder + "p_lvls.txt", delimiter=',')


# reverse ordered going from 400km down to p0
plvls = plvls[::-1]

for p0 in plvls:
    p0_infolder = infolder + 'plvl_' + str(p0) + '/'

    # first, load the all the data
    data = np.loadtxt(p0_infolder + species + "_dens_perc_diff_400_pdrag.txt", delimiter=',')
    z0_alts = np.loadtxt(p0_infolder + "zp_altitudes_km_pdrag.txt", delimiter=',')              # the actual altitude values [km] for the p0 lvls

    # find the altitudes for this specific p0 map that caused an ~10% change
    for i in range(0, 72):
        for j in range(0, 144):
            if (abs(data[i, j]) > bot_lim) and (alt_10_perc[i, j] == 0):
                alt_10_perc[i,j] = z0_alts[i,j]


# smooth function with 1 sigma gaussian filter
alt_smooth = ndimage.gaussian_filter(alt_10_perc, sigma=1, order=0)

# -----  Percent differences -----------------------
plt.figure()
cont = plt.contourf(X, Y, alt_smooth, 10, colors='k')
if species == 'HE':
    levels = np.linspace(110, np.max(alt_smooth))
elif species == 'N2':
    levels = np.linspace(110, np.max(alt_smooth))
myplot = plt.contourf(X, Y, alt_smooth, levels, cmap='jet')
cbar = plt.colorbar(myplot, format='%.0f')
cbar.ax.set_ylabel('Lower Boundary Integration Altitude [km]')
plt.title('Altitudes for ' + str(bot_lim) + '% Change of ' + species + ' Diffusive Equilibrium at ~400km')

plt.xticks(np.arange(0., 24., 3.), LT_ticks)
plt.plot(magnetic_x, mag_equator[:, 1], 'r')  # Add magnetic equator line
plt.xlabel('Local Solar Time [hr]')
plt.ylabel('Latitude [deg]')
plt.savefig('./altitude_cont/10percChange_Figs/' + species + '_alt10perc_alpha_' + alpha + '.png')
plt.show()
