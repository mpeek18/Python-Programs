# -*- coding: utf-8 -*-
"""
Created on Tue May 29 21:05:05 2018

@author: Matthew Peek
Last Modified: 8 June 2018
"""
import numpy as np
from astropy.io import ascii
import astropy.io.fits as fits
from matplotlib import pyplot as plt
"""
Algorithm:
    read in Absorption_Data.dat file
    get galID column
    
    for i in galID:
        fileName = 'CROP-SDSS-J120639.85+025308.3-G141_' + str(ID).zfill(5) + '.2d.fits'
        fitsImage = fits file data
        finalImage = numpy 2d image sum
        fits.writeto(image name)   
"""

#Open Absorber_Data.dat file and get galaxy id's, get image file name, open fits data,
#use numpy to stack images by summing them together. Write stacked image to new fits file.
count = 0
absorberFile = ascii.read('Absorption_Data.dat', delimiter='|')
ID = absorberFile['col2']

for i in range(1, len(ID)):
    #if (ID[i] == '341' or ID[i] == '410' or ID[i] == '352' or ID[i] == '333'):
    try:
        fileName = 'CROP-SDSS-J120639.85+025308.3-G141_00' + ID[i] + '.2d.fits'
        f = fits.open(fileName)
        image = f[0].data
        #print (image)
        fileList = [fits.getdata(fileName) for image in fileName[i]]
        #fileList = [fits.getdata(fileName)]
        #print (fileList)
        print (ID[i])
        print (image.shape)
        
        count += 1
    except IOError:
        print ("Image ID " + ID[i] + " not found!")
"""
for x in range(0, len(image)):
    for y in range(0, len(image[i])):
        image[x][y] += image[x][y]
        finalImage = image
print ("Final", finalImage)
"""
print ("fileList:", fileList)
stack = np.stack(fileList, axis=0)
finalImage = np.mean(stack, axis=0)

#finalImage = np.sum(fileList, axis=0)
fits.writeto('Field9_Stacked_Image.fits', finalImage, f[0].header, overwrite=True)
f.close()
print (count)

#Show stacked image.
plt.imshow(finalImage)
plt.colorbar()