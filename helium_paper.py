#!/usr/local/bin/python3

import pyglow
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np



# # #==============================================================================
# # # Creating the TIEGCM Vs. MSIS Helium Number Density at 400 km geometric
# # #==============================================================================
# #
# # # LOAD TIEGCM GEOMETRIC DATA
# # #==============================================================================
#
# HE number density [1/cm^3] for 2.5 degree size with ion drag at 400 km geometric altitude
TIEGCMfile = '/Users/hannahholt/Documents/MATLAB/TIEGCM/Contour_textfiles/geom_coord/HE_alpha_-0.38/HE_dens_model_400km_pdrag.txt'

# Load TIEGCM data file
TIEGCM_data = np.loadtxt(TIEGCMfile,delimiter=',')


# CALL MSIS
#==============================================================================
lat = 0.
lon = 0.
dn = datetime(1992, 3, 20, 0, 2)#year,month,day,hour,minute

MSIS_data = np.zeros((72,144))
marklon = 0
gridspace = 2.5         # 2.5 degree grid spacing

#Loop through lat and lon. select what MSIS values are desired
for Lon in range(-72,72,1) :
    marklat = 0
    for Lat in range(-36,36,1):
        pt = pyglow.Point(dn, Lat*2.5, Lon*2.5, 400) #2.5 degree grid size
        result = pt.run_msis()
        MSIS_data[marklat,marklon] = result.nn['HE']   # helium number density
        marklat = marklat+1
    marklon = marklon+1

F107 = pt.f107
Ap = pt.ap


# LOAD MAGNETIC EQUATOR PTS
#==============================================================================
mag_equator = np.loadtxt("Magnetic_equator_lat_lon.txt", delimiter=',')
magnetic_x = np.linspace(0, 143 / 6, 360)

# Create grid
x = np.linspace(0, 143 / 6, 144)
y = np.linspace(-88.75, 88.75, 72)
X, Y = np.meshgrid(x, y)

LT_ticks = np.array([12., 15., 18., 21., 0., 3., 6., 9., 12.])

# PLOT FIGURE 1 - TIEGCM vs MSIS at 400 km geometric
# ==============================================================================

DATmin = np.amin(TIEGCM_data)
DATmax = np.amax(TIEGCM_data)

fig = plt.figure(figsize=(14, 5))  # makes the figure bigger

plt.subplot(121)
levels = np.linspace(DATmin, DATmax, 100)
myplot = plt.contourf(X, Y, TIEGCM_data, levels, cmap='jet')
cont = plt.contour(X, Y, TIEGCM_data, 10, colors='k')
cbar = plt.colorbar(myplot, format='%.2e')
cbar.ax.set_ylabel(r'Number Density [cm$^{-3}$]')

plt.plot(magnetic_x,mag_equator[:,1],'r') #Add magnetic equator line
plt.xticks(np.arange(0., 24., 3.), LT_ticks)
plt.title('TIEGCM Helium Density 400km UT=0 with Ion Drag')
plt.xlabel('Local Solar Time [hr]')
plt.ylabel('Latitude [deg]')

plt.subplot(122)
levels = np.linspace(DATmin, DATmax, 100)
myplot = plt.contourf(X, Y, MSIS_data, levels, cmap='jet')
cont = plt.contour(X, Y, MSIS_data, 10, colors='k')
cbar = plt.colorbar(myplot, format='%.2e')
cbar.ax.set_ylabel(r'Number Density [cm$^{-3}$]')

plt.plot(magnetic_x, mag_equator[:,1],'r') #Add magnetic equator line
plt.xticks(np.arange(0., 24., 3.), LT_ticks)
plt.title('MSIS Helium Density 400km UT=0, F10.7=167.2, Ap=3')
plt.xlabel('Local Solar Time [hr]')
plt.savefig('./helium_paper_figs/MSISvsTIEGCM_He_400km.pdf')
plt.show()



#==============================================================================
# Comparing the TIEGCM number density for O, N2 and HE at ~400 km (Zp = 3.375)
#==============================================================================

#
# # Number densities [1/cm^3] for 2.5 degree size with ion drag at ~400 km geopotential altitude *** THIS IS IN PRESSURE LVLS!!!
# TIEGCMfolder = '/Users/hannahholt/Documents/MATLAB/TIEGCM/Contour_textfiles/pressure/'
# HE_file = 'HE_alpha_-0.38/HE_dens_model_400_pdrag.txt'
# N2_file = 'N2_alpha_0/N2_dens_model_400_pdrag.txt'
# O_file = 'O_alpha_0/O_dens_model_400_pdrag.txt'
#
# # LOAD TIEGCM DATA FILES
# #==============================================================================
# HE_data = np.loadtxt(TIEGCMfolder + HE_file, delimiter=',')
# N2_data = np.loadtxt(TIEGCMfolder + N2_file, delimiter=',')
# O_data = np.loadtxt(TIEGCMfolder + O_file, delimiter=',')
#
# # LOAD MAGNETIC EQUATOR PTS
# #==============================================================================
# mag_equator = np.loadtxt("Magnetic_equator_lat_lon.txt", delimiter=',')
# magnetic_x = np.linspace(0, 143 / 6, 360)
#
# # Create grid
# x = np.linspace(0, 143 / 6, 144)
# y = np.linspace(-88.75, 88.75, 72)
# X, Y = np.meshgrid(x, y)
#
# LT_ticks = np.array([12., 15., 18., 21., 0., 3., 6., 9., 12.])
#
# # PLOT FIGURES
# #==============================================================================
#
#
# fig = plt.figure(figsize=(15, 4.5))  # makes the figure bigger
#
# plt.subplot(131)
# DATmin = np.amin(HE_data)
# DATmax = np.amax(HE_data)
# levels = np.linspace(DATmin, DATmax, 100)
# myplot = plt.contourf(X, Y, HE_data, levels, cmap='jet')
# cont = plt.contour(X, Y, HE_data, 10, colors='k')
# cbar = plt.colorbar(myplot, format='%.1e')
#
#
# plt.plot(magnetic_x,mag_equator[:,1],'r') #Add magnetic equator line
# plt.xticks(np.arange(0., 25., 3.), LT_ticks)
# plt.title('He Density')
# plt.xlabel('Local Solar Time [hr]')
# plt.ylabel('Latitude [deg]')
#
#
# plt.subplot(132)
# DATmin = np.amin(N2_data)
# DATmax = np.amax(N2_data)
# levels = np.linspace(DATmin, DATmax, 100)
# myplot = plt.contourf(X, Y, N2_data, levels, cmap='jet')
# cont = plt.contour(X, Y, N2_data, 10, colors='k')
# cbar = plt.colorbar(myplot, format='%.1e')
#
# plt.plot(magnetic_x,mag_equator[:,1],'r') #Add magnetic equator line
# plt.xticks(np.arange(0., 25., 3.), LT_ticks)
# plt.title('N2 Density')
# plt.xlabel('Local Solar Time [hr]')
#
# plt.subplot(133)
# DATmin = np.amin(O_data)
# DATmax = np.amax(O_data)
# levels = np.linspace(DATmin, DATmax, 100)
# myplot = plt.contourf(X, Y, O_data, levels, cmap='jet')
# cont = plt.contour(X, Y, O_data, 10, colors='k')
# cbar = plt.colorbar(myplot, format='%.1e')
# cbar.ax.set_ylabel('Number Density [1/cm^3]')
#
# plt.plot(magnetic_x,mag_equator[:,1],'r') #Add magnetic equator line
# plt.xticks(np.arange(0., 25., 3.), LT_ticks)
# plt.title('O Density')
# plt.xlabel('Local Solar Time [hr]')
#
# plt.suptitle('TIEGCM Number Densities at ~400km UT=0 with Ion Drag', fontsize=12, fontweight='bold')
# plt.tight_layout(pad=3,w_pad=.7)
# plt.savefig('./helium_paper_figs/He_N2_O_~400km.pdf')
# plt.show()


#==============================================================================
# TIGCM Neutral Temerature and Vertical Winds at ~400 km (Zp = 3.375)
#==============================================================================

# Adapted from Tiegcm_contour.py
# Use Contour_Plots_V2.m to produce the text files

# foldername = '/Users/hannahholt/Documents/MATLAB/TIEGCM/Contour_textfiles/pressure/'
#
# Zp_name = '3,375'
# name1 = 'Neutral_Temp_'
# name2 = 'Vert_Winds[ms]_'
#
# LT_ticks = np.array([12., 15., 18., 21., 0., 3., 6., 9., 12.])
# Zp_real = Zp_name.replace(',', '.')
# filetype = 'pdrag'
# description = 'with Ion Drag'
#
# # LOAD TIEGCM DATA FILES
# #==============================================================================
# Temp_data = np.loadtxt(foldername + name1 + "Zp_" + filetype + '_' + Zp_name + ".txt", delimiter=',')
# Wind_data = np.loadtxt(foldername + name2 + "Zp_" + filetype + '_' + Zp_name + ".txt", delimiter=',')
# label = '$Z_{p}=$'+ Zp_real  # for graph labeling
#
#
# # LOAD MAGNETIC EQUATOR PTS
# #==============================================================================
# mag_equator = np.loadtxt("Magnetic_equator_lat_lon.txt", delimiter=',')
# magnetic_x = np.linspace(0, 143 / 6, 360)
#
#
# # Create contour plot mesh of correct size
# x = np.linspace(0,143/6, 144)               # longitudes
# y = np.linspace(-88.75,88.75,72)            # latitudes
# X, Y = np.meshgrid(x, y)
#
# # PLOT FIGURES
# #==============================================================================
# fig = plt.figure(figsize=(13, 4.5))  # makes the figure bigger
#
# plt.subplot(121)
# DATmin = np.amin(Temp_data)
# DATmax = np.amax(Temp_data)
#
# levels = np.linspace(DATmin, DATmax, 100)
# # levels = np.linspace(890,1325,100)
# myplot = plt.contourf(X, Y, Temp_data, levels, cmap='jet')
# cont = plt.contour(X, Y, Temp_data, 10, colors='k')
# cbar = plt.colorbar(myplot, format='%d')
# cbar.ax.set_ylabel('Neutral Temp [K]')
# plt.title('TIEGCM Neutral Temp at ~400km UT=0 ' + description)
# plt.plot(magnetic_x, mag_equator[:,1],'r') #Add magnetic equator line
# plt.xticks(np.arange(0.,25, 3.), LT_ticks)
# plt.xlabel('Local Solar Time [hr]')
# plt.ylabel('Latitude [deg]')
#
#
# plt.subplot(122)
# DATmin = np.amin(Wind_data)
# DATmax = np.amax(Wind_data)
#
# levels = np.linspace(DATmin, DATmax, 100)
# # levels = np.linspace(-26,26,100)
# myplot = plt.contourf(X, Y, Wind_data, levels, cmap='seismic')
# cont = plt.contour(X, Y, Wind_data, 10, colors='k')
# cbar = plt.colorbar(myplot, format='%.1f')
# cbar.ax.set_ylabel('Vertical Winds [m/s]')
# plt.title('TIEGCM Vertical Winds at ~400km UT=0 ' + description)
# plt.plot(magnetic_x,mag_equator[:,1],'r') # Add magnetic equator line
# plt.xticks(np.arange(0.,25, 3.), LT_ticks)
# plt.xlabel('Local Solar Time [hr]')
#
# plt.savefig('./helium_paper_figs/Neut_Temp+Winds_~400km.pdf')
# plt.show()


#==============================================================================
# TIEGCM He density percent difference w/ and w/o thermal diffusion
#==============================================================================


species = 'HE'                  # which species
alpha = ['-0.38', '0']          # thermal diffusion factor

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

i = 4       # use plvl = -3.625


p0_infolder1 = infolder1 + 'plvl_' + str(plvls[i]) + '/'  # input folders for text files from every z0 point for thermal diff.
p0_infolder2 = infolder2 + 'plvl_' + str(plvls[i]) + '/'  # input folders for text files from every z0 point for NO therm diff

# first, load the all the data
data_integrated_therm = np.loadtxt(p0_infolder1 + species + "_dens_Diffusive_400_pdrag.txt", delimiter=',')
data_integrated_notherm = np.loadtxt(p0_infolder2 + species + "_dens_Diffusive_400_pdrag.txt", delimiter=',')

# find the percDiff between the integrated plots
data_differenced = ( np.divide(data_integrated_therm, data_integrated_notherm) - 1 ) * 100
mind = 30
maxd = 60

fig = plt.figure() # makes the figure bigger

# ------------  Percent differences -----------------------
levels = np.linspace(mind, maxd, 300)                       # set your beginging limits for colorbar
cont = plt.contour(X, Y, data_differenced, 12, colors='k')
myplot = plt.contourf(X, Y, data_differenced, levels, cmap='jet')
cbar = plt.colorbar(myplot, format='%.0f')
cbar.ax.set_ylabel('Perc. Diff. With and Without Thermal Diffusion')
plt.title('Helium Density Percent Difference at ~400 km\n' + r'Z$_0$ = ' + str(plvls[i]))

plt.xticks(np.arange(0., 24., 3.), LT_ticks)
plt.plot(magnetic_x, mag_equator[:, 1], 'r')  # Add magnetic equator line
plt.xlabel('Local Solar Time [hr]')
plt.ylabel('Latitude [deg]')

plt.savefig('./helium_paper_figs/He_dens_ThermVSnoTherm.pdf')

plt.show()
plt.close()


