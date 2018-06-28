# -*- coding: utf-8 -*-
"""
Created on Tue May 29 21:05:05 2018

@author: Matthew Peek
Last Modified: 28 June 2018
All Fields Image Stack
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
Normalize image function, takes image as an argument, gets data and stores in
numpy array. Sum all data in image then divides each pixel by image sum.
returns normalized image as numpy array.
"""
def imageNorm(fileName):
    dataList = []
    data = [fits.getdata(fileName)]
    dataList.append(data)
    sumData = np.sum(data)
    print ("Summed Image Data:", sumData,'\n')
    print ("Data:", data,'\n')
        
    for i in range(0, len(dataList)):
        for j in range(0, len(dataList[0])):
            print ("Before:", dataList[i][j])
            normed = (dataList[i][j] / sumData)
    
    print ("Normed:", normed)    
    print ("Normalization complete!")     
    return normed
#End imageNorm function
    

"""
Stack function, takes numpy array as argument, loops through array argument
getting fits data and combining all image data into same array.

Stack image data standard, mean averaging, and median averaging, compare results.
Write results to new fits files.
"""
def stack(fileList):
    imageData = [fits.getdata(file) for file in fileList]
    #imageData = [file for file in fileList]
    
    print ("Total image data:", imageData,'\n')
    meanImage = np.mean(imageData, axis=0)
    medianImage = np.median(imageData, axis=0)
    imageStack = np.sum(imageData, axis=0)
    
    fits.writeto('Stacked_Image_Mean.fits', meanImage, overwrite=True)
    fits.writeto('Stacked_Image_Median.fits', medianImage, overwrite=True)
    fits.writeto('Stacked_Image.fits', imageStack, overwrite=True)
    
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
#End Stack function
    
"""
Define list containing field numbers, go through list and read in fields. Call
stack function to stack all fields. 
"""
count = 0
fileList = []
galList = [1, 2, 3, 4, 5, 7, 8, 9]
for i in range(0, len(galList)):
    try:
        fileName = 'Field' + str(galList[i]) + '_Stacked_Image.fits'
        
        #Call imageNorm function.
        #normed = imageNorm(fileName)
        
        #Append normalized image to fileList to pass as argument to stack function.
        fileList.append(fileName)
        file = fits.open(fileName)
        image = file[0].data
        file.close()
        
        print (galList[i])
        print (image.shape,'\n')
        
        count += 1
    except IOError:
        print ("Image ID " + (galList[i]) + " not found!")
print ("Number of images processed:", count,'\n')

#Call Stack function
stack(fileList) 
