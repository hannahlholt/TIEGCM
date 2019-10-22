#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 17:19:23 2017

Updated By HLH - 1/4/2019 to account for multiple starting altitudes for integrating the diffusive species density values up to 400. These are compared with the model helium ouput at 400 km to see how they differ depending on the starting integration altitude

USE with make_video.py to show how the

"""


import matplotlib.pyplot as plt
import numpy as np
import os, sys
import scipy.ndimage as ndimage

species = 'O'              # which species
alpha = '0'             # thermal diffusion factor

if species == 'HE' and alpha == '0':
    mind = -170
    maxd = 170
elif species == 'HE' and alpha == '-0.38':
    mind = -75
    maxd = 75
elif species == 'N2':
    mind = -50
    maxd = 50
elif species == 'O':
    mind = -70
    maxd = 70

name = '3.375'



infolder = '/Users/hannahholt/Documents/MATLAB/TIEGCM/Contour_textfiles/pressure/' + species + '_alpha_' + alpha + '/'


# first load the of lower boundary pressure level points
plvls = np.loadtxt(infolder + "p_lvls.txt", delimiter=',')

# Load magnetic equator data points to be overlayed on plot
mag_equator = np.loadtxt("Magnetic_equator_lat_lon.txt", delimiter=',')
marklon = 0
magnetic_x = np.linspace(0, 143 / 6, 360)

LT_ticks = np.array([12., 15., 18., 21., 0., 3., 6., 9., 12.])
x = np.linspace(0, 143 / 6, 144)
y = np.linspace(-88.75, 88.75, 72)
X, Y = np.meshgrid(x, y)

i = 0
for p in plvls:
    p0_infolder = infolder + 'plvl_' + str(p) + '/'          # input folders for text files from every pressure level point
    p0_outfolder = '/Users/hannahholt/Documents/PYTHON/TIEGCM/Contour_Movie_Figs/pressureV2/' + species + '/alpha_' + alpha + '/'     # where to save the figures

    if os.path.exists(p0_outfolder) == 0:
        os.mkdir(p0_outfolder)

    # first, load the all the data
    data_differenced = np.loadtxt(p0_infolder + species + "_dens_perc_diff_400_pdrag.txt", delimiter=',')
    data_integrated = np.loadtxt(p0_infolder + species + "_dens_Diffusive_400_pdrag.txt", delimiter=',')
    data_model_top = np.loadtxt(infolder + species + "_dens_model_400_pdrag.txt", delimiter=',')
    data_model_bottom = np.loadtxt(p0_infolder + species + "_dens_model_z0_pdrag.txt", delimiter=',')

    fig = plt.figure(figsize=(14, 7.75)) # makes the figure bigger
    plt.suptitle('TIEGCM ' + species + ' DENSITIES, P-Coord, [UT = 0, pdrag, alpha = ' + alpha + ']')

    wspace = 0.4  # the amount of width reserved for blank space between subplots
    hspace = 0.4  # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=wspace, hspace=hspace)

    # -----  Percent differences -----------------------
    plt.subplot(221)
    # set a resonable amount of contour lines
    if int( np.round((np.amax(data_differenced) - np.amin(data_differenced)) / 7) ) < 5:
        num_cont = 6
    else:
        num_cont = 12
    cont = plt.contour(X, Y, data_differenced, num_cont, colors='k')

    levels = np.linspace(mind, maxd, 300)                       # set your beginging limits for colorbar
    myplot = plt.contourf(X, Y, data_differenced, levels, cmap='seismic')
    cbar = plt.colorbar(myplot, format='%.0f')
    cbar.ax.set_ylabel('Percent Difference')
    plt.title('Deviation from Diffusive Equilibrium at ~400km, Zp = ' + name)

    plt.xticks(np.arange(0., 24., 3.), LT_ticks)
    plt.plot(magnetic_x, mag_equator[:, 1], 'r')  # Add magnetic equator line
    plt.xlabel('Local Solar Time [hr]')
    plt.ylabel('Latitude [deg]')

    # ------ Integrated Species Density Assuming Pure Diffusive Profiles -----
    plt.subplot(222)
    levels = np.linspace(np.amin(data_integrated), np.amax(data_integrated), 100)
    cont = plt.contour(X, Y, data_integrated, 10, colors='k')
    myplot = plt.contourf(X, Y, data_integrated, levels, cmap='jet')
    cbar = plt.colorbar(myplot, format='%.1e')
    cbar.ax.set_ylabel('Number Density [#/cm^3]')
    plt.title('Integrated Diffusive Profile at ~400km, Zp = ' + name)

    plt.xticks(np.arange(0., 24., 3.), LT_ticks)
    plt.plot(magnetic_x, mag_equator[:, 1], 'r')  # Add magnetic equator line
    plt.xlabel('Local Solar Time [hr]')
    plt.ylabel('Latitude [deg]')

    # ---- Model Species Density at ~400 km ------------------------------
    plt.subplot(223)
    levels = np.linspace(np.amin(data_model_top), np.amax(data_model_top), 100)
    cont = plt.contour(X, Y, data_model_top, 10, colors='k')
    myplot = plt.contourf(X, Y, data_model_top, levels, cmap='jet')
    cbar = plt.colorbar(myplot, format='%.1e')
    cbar.ax.set_ylabel('Number Density [#/cm^3]')
    plt.title('Model Output at ~400km, Zp = ' + name)

    plt.xticks(np.arange(0.,24.,3.), LT_ticks)
    plt.plot(magnetic_x,mag_equator[:,1],'r') #Add magnetic equator line
    plt.xlabel('Local Solar Time [hr]')
    plt.ylabel('Latitude [deg]')

    # ---- Model Species Density at z0 km ------------------------------
    plt.subplot(224)
    levels = np.linspace(np.amin(data_model_bottom), np.amax(data_model_bottom), 100)
    cont = plt.contour(X, Y, data_model_bottom, 10, colors='k')
    myplot = plt.contourf(X, Y, data_model_bottom, levels, cmap='jet')
    cbar = plt.colorbar(myplot, format='%.1e')
    cbar.ax.set_ylabel('Number Density [#/cm^3]')
    plt.title(' Model Output at Lower Boundary Zp = ' + str(p))

    plt.xticks(np.arange(0., 24., 3.), LT_ticks)
    plt.plot(magnetic_x, mag_equator[:, 1], 'r')  # Add magnetic equator line
    plt.xlabel('Local Solar Time [hr]')
    plt.ylabel('Latitude [deg]')

    #plt.show()
    if i < 10:
        plt.savefig(p0_outfolder + '0' + str(i) + '_p0_' + str(p) + '.png')
    else:
        plt.savefig(p0_outfolder + str(i) + '_p0_' + str(p) + '.png')

    plt.close()

    i = i+1

# end for p in plvls




