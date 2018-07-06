# -*- coding: utf-8 -*-
"""
Created on Tue May 29 21:05:05 2018

@author: Matthew Peek
Last Modified: 6 July 2018
Field 9 Image Stack
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
def imageNormAbsorber(fileName):
    dataList = []
    data = [fits.getdata(fileName)]
    dataList.append(data)
    sumData = np.sum(data)
    print ("Summed Image Data:", sumData,'\n')
    print ("Data:", data,'\n')
        
    for i in range(0, len(dataList)):
        for j in range(0, len(dataList[0])):
            print ("Before:", dataList[i][j],'\n')
            normed = (dataList[i][j] / sumData)
    
    print ("Normed:", normed)   
    print ("Normed Sum:", np.sum(normed))
    print ("Normalization complete!")     
    return normed
#End Absorber imageNorm function

def imageNormNonAbsorber(fileName):
    dataList = []
    data = [fits.getdata(fileName)]
    dataList.append(data)
    sumData = np.sum(data)
    print ("Summed Image Data:", sumData,'\n')
    print ("Data:", data,'\n')
        
    for i in range(0, len(dataList)):
        for j in range(0, len(dataList[0])):
            print ("Before:", dataList[i][j],'\n')
            normed = (dataList[i][j] / sumData)
    
    print ("Normed:", normed)   
    print ("Normed Sum:", np.sum(normed))
    print ("Normalization complete!")     
    return normed
#End Non-Absorber imageNorm function
    

"""
Stack function, takes numpy array as argument, loops through array argument
getting fits data and combining all image data into same array.

Stack image data standard, mean averaging, and median averaging, compare results.
Write results to new fits files.
"""
def stackAbsorber(fileList):
    #imageData = [fits.getdata(file) for file in fileList]
    imageData = [file for file in fileList]
    
    print ("Total Absorber image data:", imageData,'\n')
    meanImage = np.mean(imageData, axis=0)
    medianImage = np.median(imageData, axis=0)
    imageStack = np.sum(imageData, axis=0)
    
    fits.writeto('Field9_Stacked_Image_Mean_Normed_Absorber.fits', meanImage, overwrite=True)
    fits.writeto('Field9_Stacked_Image_Median_Normed_Absorber.fits', medianImage, overwrite=True)
    fits.writeto('Field9_Stacked_Image_Normed_Absorber.fits', imageStack, overwrite=True)
    
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
    
    print ("Stacking Absorbers complete!")
#End Absorber Stack function

def stackNonAbsorber(fileList):
    #imageData = [fits.getdata(file) for file in fileList]
    imageData = [file for file in fileList]
    
    print ("Total Non-Absorber image data:", imageData,'\n')
    meanImage = np.mean(imageData, axis=0)
    medianImage = np.median(imageData, axis=0)
    imageStack = np.sum(imageData, axis=0)
    
    fits.writeto('Field9_Stacked_Image_Mean_Normed_NonAbsorber.fits', meanImage, overwrite=True)
    fits.writeto('Field9_Stacked_Image_Median_Normed_NonAbsorber.fits', medianImage, overwrite=True)
    fits.writeto('Field9_Stacked_Image_Normed_NonAbsorber.fits', imageStack, overwrite=True)
    
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
    
    print ("Stacking Non-Absorbers complete!")
#End Non-Absorber Stack function

    
"""
Open Absorber_Data.dat file and get galaxy id's, get image file name, open fits data.
Start program by reading in id's and appending them to new list.
"""
absorberFile = ascii.read('Absorption_Data.dat', delimiter='|')
ID = absorberFile['col2']
absorber = absorberFile['col7']
countAbsorber = 0
countNonAbsorber = 0
fileListAbsorb = []
fileListNonAbsorb = []
for i in range(1, len(ID)):
    if (absorber[i] == 'Yes'):
        try:
            fileName = 'CROP-SDSS-J120639.85+025308.3-G141_00' + ID[i] + '.2d.fits'
            #Call imageNorm function.
            normed = imageNormAbsorber(fileName)
        
            #Append normalized image to fileList to pass as argument to stack function.
            fileListAbsorb.append(normed)
            file = fits.open(fileName)
            image = file[0].data
            file.close()
        
            print (ID[i])
            print (image.shape,'\n')
        
            countAbsorber += 1
        except IOError:
                print ("Image ID " + ID[i] + " not found!")

    else:
        try:
            fileName = 'CROP-SDSS-J120639.85+025308.3-G141_00' + ID[i] + '.2d.fits'
            #Call imageNorm function.
            normed = imageNormNonAbsorber(fileName)
        
            #Append normalized image to fileList to pass as argument to stack function.
            fileListNonAbsorb.append(normed)
            file = fits.open(fileName)
            image = file[0].data
            file.close()
        
            print (ID[i])
            print (image.shape,'\n')
        
            countNonAbsorber += 1
        except IOError:
                print ("Image ID " + ID[i] + " not found!")
        
#Call stack functions
stackAbsorber(fileListAbsorb)
stackNonAbsorber(fileListNonAbsorb)
print ("Number of Absorbers Processed:", countAbsorber)
print ("Number of Non-Absorbers Processed:", countNonAbsorber)
