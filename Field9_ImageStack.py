# -*- coding: utf-8 -*-
"""
Created on Tue May 29 21:05:05 2018

@author: Matthew Peek
Last Modified: 21 June 2018
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

"""
Stack function, takes numpy array as argument, loops through array argument
getting fits data and combining all image data into same array.

Stack image data standard, mean averaging, and median averaging, compare results.
Write results to new fits files.
"""
def stack(fileList):
    imageData = [fits.getdata(file) for file in fileList]
    
    print ("Total image data:", imageData,'\n')
    meanImage = np.mean(imageData, axis=0)
    medianImage = np.median(imageData, axis=0)
    imageStack = np.sum(imageData, axis=0)
    
    fits.writeto('Field9_Stacked_Image_Mean.fits', meanImage, overwrite=True)
    fits.writeto('Field9_Stacked_Image_Median.fits', medianImage, overwrite=True)
    fits.writeto('Field9_Stacked_Image.fits', imageStack, overwrite=True)
    
    print ("Image Mean:", meanImage,'\n')
    print ("Image Median:", medianImage,'\n')
    print ("Image Sum:", imageStack,'\n')
    
    plt.clf()
    plt.imshow(meanImage)
    plt.colorbar()
    
    plt.clf()
    plt.imshow(medianImage)
    plt.colorbar()
    
    plt.clf()
    plt.imshow(imageStack)
    plt.colorbar()
    print ("Stacking complete!")


#Open Absorber_Data.dat file and get galaxy id's, get image file name, open fits data.
#Start program by reading in id's and appending them to new list.
count = 0
absorberFile = ascii.read('Absorption_Data.dat', delimiter='|')
ID = absorberFile['col2']

fileList = []
for i in range(1, len(ID)):
    try:
        fileName = 'CROP-SDSS-J120639.85+025308.3-G141_00' + ID[i] + '.2d.fits'
        file = fits.open(fileName)
        fileList.append(fileName)
        image = file[0].data
        file.close()
        
        print (ID[i])
        print (image.shape)
        
        count += 1
    except IOError:
        print ("Image ID " + ID[i] + " not found!")
print ("Number of images processed:", count,'\n')

#Call Stack function
stack(fileList) 
