#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 2019 09:40:50 by Hannah Holt

This program runs very similar to TIEGCM_contour_percDiff.py except it compares the integrated profiles w/ and w/o thermal
diffusion and looks at the percent difference between the two.

USE with make_video.py
"""


import matplotlib.pyplot as plt
import numpy as np
import os, sys
import scipy.ndimage as ndimage

species = 'HE'                  # which species
alpha = ['-0.38', '0']          # thermal diffusion factor
name = '3.375'
mind = -100
maxd = 100

infolder1 = '/Users/hannahholt/Documents/MATLAB/TIEGCM/Contour_textfiles/pressure/' + species + '_alpha_' + alpha[0] + '/'
infolder2 = '/Users/hannahholt/Documents/MATLAB/TIEGCM/Contour_textfiles/pressure/' + species + '_alpha_' + alpha[1] + '/'

# first load the of lower boundary pressure level points
plvls = np.loadtxt(infolder1 + "p_lvls.txt", delimiter=',')

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
    p0_infolder1 = infolder1 + 'plvl_' + str(p) + '/'  # input folders for text files from every z0 point for thermal diff.
    p0_infolder2 = infolder2 + 'plvl_' + str(p) + '/'  # input folders for text files from every z0 point for NO therm diff

    # first, load the all the data
    data_integrated_therm = np.loadtxt(p0_infolder1 + species + "_dens_Diffusive_400_pdrag.txt", delimiter=',')
    data_integrated_notherm = np.loadtxt(p0_infolder2 + species + "_dens_Diffusive_400_pdrag.txt", delimiter=',')
    data_model_bottom = np.loadtxt(p0_infolder1 + species + "_dens_model_z0_pdrag.txt", delimiter=',')

    # find the percDiff between the integrated plots
    data_differenced = ( np.divide(data_integrated_therm, data_integrated_notherm) - 1 ) * 100

    if i > 30:
       data_differenced = ndimage.gaussian_filter(data_differenced, sigma=1, order=0)

    p0_outfolder = '/Users/hannahholt/Documents/PYTHON/TIEGCM/Contour_Movie_Figs/pressureV2/' + species + '/thermal_diff_percDiff/'     # where to save the figures

    if os.path.exists(p0_outfolder) == 0:
        os.mkdir(p0_outfolder)

    fig = plt.figure(figsize=(14, 7.75)) # makes the figure bigger
    plt.suptitle('TIEGCM ' + species + ' DENSITIES, P-Coord, [UT = 0, pdrag]')

    wspace = 0.4  # the amount of width reserved for blank space between subplots
    hspace = 0.4  # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=wspace, hspace=hspace)

    # ------------  Percent differences -----------------------
    plt.subplot(221)
    levels = np.linspace(mind, maxd, 300)                       # set your beginging limits for colorbar
    cont = plt.contour(X, Y, data_differenced, 10, colors='k')
    myplot = plt.contourf(X, Y, data_differenced, levels, cmap='seismic')
    cbar = plt.colorbar(myplot, format='%.0f')
    cbar.ax.set_ylabel('Percent Difference')
    plt.title('Percent Difference - Therm vs No Therm')

    plt.xticks(np.arange(0., 24., 3.), LT_ticks)
    plt.plot(magnetic_x, mag_equator[:, 1], 'r')  # Add magnetic equator line
    plt.xlabel('Local Solar Time [hr]')
    plt.ylabel('Latitude [deg]')

    # ------ Integrated Species Density Assuming Pure Diffusive Profiles - WITH THERM DIFF -----
    plt.subplot(222)
    levels = np.linspace(np.amin(data_integrated_therm), np.amax(data_integrated_therm), 100)
    cont = plt.contour(X, Y, data_integrated_therm, 10, colors='k')
    myplot = plt.contourf(X, Y, data_integrated_therm, levels, cmap='jet')
    cbar = plt.colorbar(myplot, format='%.1e')
    cbar.ax.set_ylabel('Number Density [#/cm^3]')
    plt.title('Diffusive Profile WITH Thermal Diff. at ~400km, Zp = ' + name)

    plt.xticks(np.arange(0., 24., 3.), LT_ticks)
    plt.plot(magnetic_x, mag_equator[:, 1], 'r')  # Add magnetic equator line
    plt.xlabel('Local Solar Time [hr]')
    plt.ylabel('Latitude [deg]')

    # ---- Model Species Density at z0 km -------------------------------------
    plt.subplot(223)
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

    # ----Integrated Species Density Assuming Pure Diffusive Profiles - NO THERM DIFF ------------------------------
    plt.subplot(224)
    levels = np.linspace(np.amin(data_integrated_notherm), np.amax(data_integrated_notherm), 100)
    cont = plt.contour(X, Y, data_integrated_notherm, 10, colors='k')
    myplot = plt.contourf(X, Y, data_integrated_notherm, levels, cmap='jet')
    cbar = plt.colorbar(myplot, format='%.1e')
    cbar.ax.set_ylabel('Number Density [#/cm^3]')
    plt.title('Diffusive Profile NO Thermal Diff. at ~400km, Zp = ' + name)

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


