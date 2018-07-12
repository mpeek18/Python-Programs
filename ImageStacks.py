# -*- coding: utf-8 -*-
"""
Created on Tue May 29 21:05:05 2018

@author: Matthew Peek
Last Modified: 12 July 2018
All Fields Image Stack
"""
import numpy as np
import astropy.io.fits as fits
import scipy.ndimage as ndimage
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
StackAll function, takes numpy array as argument, loops through array argument
getting fits data and combining all image data into same array.

Stack image data standard, write results to new fits files.
"""
def stackAll(fileListAll):
    imageData = [fits.getdata(file) for file in fileListAll]
    #imageData = [file for file in fileList]
    
    print ("Total image data:", imageData,'\n')
    imageStack = np.sum(imageData, axis=0)
    
    fits.writeto('Stacked_Image_All.fits', imageStack, overwrite=True)
    
    print ("Image Sum:", imageStack,'\n')
    
    plt.clf()
    plt.imshow(imageStack)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    #Gaussian Blur
    plt.clf()
    betterImage = ndimage.gaussian_filter(imageStack, sigma=(2,2), order=0)
    plt.imshow(betterImage)
    plt.savefig('Stacked_Image_All_Blur', dpi=100)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    print ("StackAll Function Complete!", '\n')
#End Stack function
    

"""
StackMeanAll function, takes numpy array as argument, loops through array argument
getting fits data and combining all image data into same array.

Stack image data mean, write results to new fits files.

""" 
def stackMeanAll(fileListMeanAll):
    imageData = [fits.getdata(file) for file in fileListMeanAll]
    print ("Total image data:", imageData,'\n')
    meanImage = np.mean(imageData, axis=0)
    
    fits.writeto('Stacked_Image_Mean_All.fits', meanImage, overwrite=True)
    
    print ("Image Mean:", meanImage,'\n')
    
    plt.clf()
    plt.imshow(meanImage)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    #Gaussian Blur
    plt.clf()
    betterImage = ndimage.gaussian_filter(meanImage, sigma=(2,2), order=0)
    plt.imshow(betterImage)
    plt.savefig('Stacked_Image_Mean_All_Blur', dpi=100)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    print ("stackMeanAll Function Complete!", '\n')
#End stackMean function
    

"""
StackMedianAll function, takes numpy array as argument, loops through array argument
getting fits data and combining all image data into same array.

Stack image data median, write results to new fits files.
"""    
def stackMedianAll(fileListMedianAll):
    imageData = [fits.getdata(file) for file in fileListMedianAll]
    print ("Total image data:", imageData,'\n')
    medianImage = np.median(imageData, axis=0)
    
    fits.writeto('Stacked_Image_Median_All.fits', medianImage, overwrite=True)
    
    print ("Image Median:", medianImage,'\n')
    
    plt.clf()
    plt.imshow(medianImage)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    #Gaussian Blur
    plt.clf()
    betterImage = ndimage.gaussian_filter(medianImage, sigma=(2,2), order=0)
    plt.imshow(betterImage)
    plt.savefig('Stacked_Image_Median_All_Blur', dpi=100)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    print ("stackMedianAll Function Complete!")
#End stackMedian function

"""
StackMeanAbsorb function, takes numpy array as argument, loops through array argument
getting fits data and combining all image data into same array.

Stack image data median, write results to new fits files.
"""    
def stackMeanAbsorb(fileListMeanAbsorb):
    imageData = [fits.getdata(file) for file in fileListMeanAbsorb]
    print ("Total image data:", imageData,'\n')
    meanAbsorbImage = np.mean(imageData, axis=0)
    
    fits.writeto('Stacked_Image_Mean_Absorb.fits', meanAbsorbImage, overwrite=True)
    
    print ("Image Mean Absorb:", meanAbsorbImage,'\n')
    
    plt.clf()
    plt.imshow(meanAbsorbImage)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    #Gaussian Blur
    plt.clf()
    betterImage = ndimage.gaussian_filter(meanAbsorbImage, sigma=(2,2), order=0)
    plt.imshow(betterImage)
    plt.savefig('Stacked_Image_Mean_Absorb_Blur', dpi=100)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    print ("stackMeanAbsorb Function Complete!")
#End stackMedian function


"""
Define list containing field numbers, go through list and read in fields. Call
stack function to stack all fields. 
"""
fileListAll = []
fileListMeanAll = []
fileListMedianAll = []
fileListMeanAbsorb = []
fileListMedianAbsorb = []
fileListMeanNonAbsorb = []
fileListMedianNonAbsorb = []
fileListStandardAbsorb = []
fileListStandardNonAbsorb = []
galList = [1, 2, 3, 4, 5, 7, 8, 9]

for i in range(0, len(galList)):
    try:
        fileNameAll = 'Field' + str(galList[i]) + '_Stacked_Image_Normed_All.fits'
        fileMeanAll = 'Field' + str(galList[i]) + '_Stacked_Image_Mean_Normed_All.fits'
        fileMedianAll = 'Field' + str(galList[i]) + '_Stacked_Image_Median_Normed_All.fits'
        fileMeanAbsorb = 'Field' + str(galList[i]) + '_Stacked_Image_Mean_Normed_Absorber.fits'
        fileNameStandAbsorb = 'Field' + str(galList[i]) + '_Stacked_Image_Normed_Absorber.fits'
        fileNameStandNonAbsorb = 'Field' + str(galList[i]) + '_Stacked_Image_Normed_NonAbsorber.fits'
        fileNameMedianAbsorb = 'Field' + str(galList[i]) + '_Stacked_Image_Median_Normed_Absorber.fits'
        fileNameMeanNonAbsorb = 'Field' + str(galList[i]) + '_Stacked_Image_Mean_Normed_NonAbsorber.fits'
        fileNameMedianNonAbsorb = 'Field' + str(galList[i]) + '_Stacked_Image_Median_Normed_NonAbsorber.fits'
        
        
        #Append normalized image to fileList to pass as argument to stack function.
        fileListAll.append(fileNameAll)
        fileListMeanAll.append(fileMeanAll)
        fileListMedianAll.append(fileMedianAll)
        fileListMeanAbsorb.append(fileMeanAbsorb)
        fileListMedianAbsorb.append(fileNameMedianAbsorb)
        fileListMeanNonAbsorb.append(fileNameMeanNonAbsorb)
        fileListMedianNonAbsorb.append(fileNameMedianNonAbsorb)
        fileListStandardAbsorb.append(fileNameStandAbsorb)
        fileListStandardNonAbsorb.append(fileNameStandNonAbsorb)
        
    except IOError:
        print ("Image ID " + str(galList[i]) + " not found!")

print (len(fileListAll))
print (len(fileListMeanAll))
print (len(fileListMedianAll))
print (len(fileListMeanAbsorb))
print (len(fileListMedianAbsorb))
print (len(fileListMeanNonAbsorb))
print (len(fileListMedianNonAbsorb))
print (len(fileListStandardAbsorb))
print (len(fileListStandardNonAbsorb))
    
#Call Stack functions
stackAll(fileListAll)
stackMeanAll(fileListMeanAll)
stackMedianAll(fileListMedianAll)
stackMeanAbsorb(fileListMeanAbsorb)
#stack(fileList)
#stackMean(fileListMean)
#stackMedian(fileListMedian) 
