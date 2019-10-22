#!/usr/local/bin/python2

'''
*** REQUIRES PYTHON 2 TO RUN CORRECTLY ***

make_video.py created by Hannah Holt 1/4/2019

This program combines multiple (consequetive) .png figures from a directory given and combines them into a .mov saved in
working directory of the program. This is currently being used with the Python diffusive equilibrium contour figs output

'''

import cv2
import os

species = 'O'              # which species
alpha = '0'                 # thermal diffusion factor

# for exporting the thermal diffusion percent difference
# image_folder = '/Users/hannahholt/Documents/PYTHON/TIEGCM/Contour_Movie_Figs/pressureV2/' + species + '/thermal_diff_percDiff/'     # where the individual figures are saves
# video_name = './Movies/pressure/Thermal_Diff_percDiff_' + species + '.mov'

# for exporting normal percent differences
image_folder = '/Users/hannahholt/Documents/PYTHON/TIEGCM/Contour_Movie_Figs/pressureV2/' + species + '/alpha_' + alpha + '/'
video_name = './Movies/pressure/Diffusive_Eq_Changes_V2_' + species + '_alpha_' + alpha + '.mov'

# for comparing the Scale heights
# image_folder = '/Users/hannahholt/Desktop/' + species + '/'
# video_name = './Movies/pressure/ProfileComp_' + species + '_alpha_' + alpha + '.mov'


# grabs the images
images = [img for img in os.listdir(image_folder) if img.endswith(".png")]

frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape                                 # assuming all figures are same size. If not, this rescales them.

fps = 2     # frames per seconds
video = cv2.VideoWriter(video_name, -1, fps, (width, height))       # Videowriter(filename, fourcc, fps, frameSize[, isColor]])

for image in images:
    video.write(cv2.imread(os.path.join(image_folder, image)))

cv2.destroyAllWindows()
video.release()
