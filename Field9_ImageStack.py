# -*- coding: utf-8 -*-
"""
Created on Tue May 29 21:05:05 2018

@author: Matthew Peek
Last Modified: 10 June 2018
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

fileList = []
for i in range(1, len(ID)):
    try:
        fileName = 'CROP-SDSS-J120639.85+025308.3-G141_00' + ID[i] + '.2d.fits'
        file = fits.open(fileName)
        image = file[0].data
        
        fileList = ([np.array(image)])
        
        print (ID[i])
        #print (fileList)
        #print (image.shape)
        
        count += 1
    except IOError:
        print ("Image ID " + ID[i] + " not found!")
file.close()

print (count)  
print ("Total fileList:", fileList)

finalImage = np.sum(fileList, axis=0)
fits.writeto('Field9_Stacked_Image.fits', finalImage, overwrite=True)
plt.imshow(finalImage)
plt.colorbar()
"""
for x in range(0, len(data)):
    for y in range(0, len(data[i])):
        data[x][y] += data[x][y]
        finalImage = data
print ("Final", finalImage)
"""
