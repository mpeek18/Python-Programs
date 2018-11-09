# -*- coding: utf-8 -*-
"""
Created on Tue May 29 21:05:05 2018

@author: Matthew Peek
Last Modified: 9 November 2018
All Fields Image Stack
"""
import numpy as np
import astropy.io.fits as fits
import scipy.ndimage as ndimage
from matplotlib import pyplot as plt
    
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
    plt.savefig('Stacked_Image_All', dpi=100)
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
#End StackAll function
    

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
    plt.savefig('Stacked_Image_Mean_All', dpi=100)
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
#End stackMeanAll function
    

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
    plt.savefig('Stacked_Image_Median_All', dpi=100)
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
#End stackMedianAll function


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
    plt.savefig('Stacked_Image_Mean_Absorb', dpi=100)
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
#End stackMeanAbsorb function
    
    
"""
StackMeanNonAbsorb function, takes numpy array as argument, loops through array argument
getting fits data and combining all image data into same array.

Stack image data median, write results to new fits files.
"""    
def stackMeanNonAbsorb(fileListMeanNonAbsorb):
    imageData = [fits.getdata(file) for file in fileListMeanNonAbsorb]
    print ("Total image data:", imageData,'\n')
    meanNonAbsorbImage = np.mean(imageData, axis=0)
    
    fits.writeto('Stacked_Image_Mean_NonAbsorb.fits', meanNonAbsorbImage, overwrite=True)
    
    print ("Image Mean Non-Absorb:", meanNonAbsorbImage,'\n')
    
    plt.clf()
    plt.imshow(meanNonAbsorbImage)
    plt.savefig('Stacked_Image_Mean_NonAbsorb', dpi=100)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    #Gaussian Blur
    plt.clf()
    betterImage = ndimage.gaussian_filter(meanNonAbsorbImage, sigma=(2,2), order=0)
    plt.imshow(betterImage)
    plt.savefig('Stacked_Image_Mean_NonAbsorb_Blur', dpi=100)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    print ("stackMeanNonAbsorb Function Complete!")
#End stackMeanNonAbsorb function

    
"""
StackMedianAbsorb function, takes numpy array as argument, loops through array argument
getting fits data and combining all image data into same array.

Stack image data median, write results to new fits files.
"""    
def stackMedianAbsorb(fileListMedianAbsorb):
    imageData = [fits.getdata(file) for file in fileListMedianAbsorb]
    print ("Total image data:", imageData,'\n')
    medianAbsorbImage = np.median(imageData, axis=0)
    
    fits.writeto('Stacked_Image_Median_Absorb.fits', medianAbsorbImage, overwrite=True)
    
    print ("Image Median Absorb:", medianAbsorbImage,'\n')
    
    plt.clf()
    plt.imshow(medianAbsorbImage)
    plt.savefig('Stacked_Image_Median_Absorb', dpi=100)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    #Gaussian Blur
    plt.clf()
    betterImage = ndimage.gaussian_filter(medianAbsorbImage, sigma=(2,2), order=0)
    plt.imshow(betterImage)
    plt.savefig('Stacked_Image_Median_Absorb_Blur', dpi=100)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    print ("stackMedianAbsorb Function Complete!")
#End stackMedianAbsorb function


"""
StackMedianNonAbsorb function, takes numpy array as argument, loops through array argument
getting fits data and combining all image data into same array.

Stack image data median, write results to new fits files.
"""    
def stackMedianNonAbsorb(fileListMedianNonAbsorb):
    imageData = [fits.getdata(file) for file in fileListMedianNonAbsorb]
    print ("Total image data:", imageData,'\n')
    medianNonAbsorbImage = np.median(imageData, axis=0)
    
    fits.writeto('Stacked_Image_Median_NonAbsorb.fits', medianNonAbsorbImage, overwrite=True)
    
    print ("Image Median Non-Absorb:", medianNonAbsorbImage,'\n')
    
    plt.clf()
    plt.imshow(medianNonAbsorbImage)
    plt.savefig('Stacked_Image_Median_NonAbsorb', dpi=100)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    #Gaussian Blur
    plt.clf()
    betterImage = ndimage.gaussian_filter(medianNonAbsorbImage, sigma=(2,2), order=0)
    plt.imshow(betterImage)
    plt.savefig('Stacked_Image_Median_NonAbsorb_Blur', dpi=100)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    print ("stackMedianNonAbsorb Function Complete!")
#End stackMedianNonAbsorb function


"""
StackStandardAbsorb function, takes numpy array as argument, loops through array argument
getting fits data and combining all image data into same array.

Stack image data median, write results to new fits files.
"""    
def stackStandardAbsorb(fileListStandardAbsorb):
    imageData = [fits.getdata(file) for file in fileListStandardAbsorb]
    print ("Total image data:", imageData,'\n')
    standardAbsorbImage = np.sum(imageData, axis=0)
    
    fits.writeto('Stacked_Image_Standard_Absorb.fits', standardAbsorbImage, overwrite=True)
    
    print ("Image Standard Absorb:", standardAbsorbImage,'\n')
    
    plt.clf()
    plt.imshow(standardAbsorbImage)
    plt.savefig('Stacked_Image_Standard_Absorb', dpi=100)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    #Gaussian Blur
    plt.clf()
    betterImage = ndimage.gaussian_filter(standardAbsorbImage, sigma=(2,2), order=0)
    plt.imshow(betterImage)
    plt.savefig('Stacked_Image_Standard_Absorb_Blur', dpi=100)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    print ("stackStandardAbsorb Function Complete!")
#End stackStandardAbsorb function


"""
StackStandardNonAbsorb function, takes numpy array as argument, loops through array argument
getting fits data and combining all image data into same array.

Stack image data median, write results to new fits files.
"""    
def stackStandardNonAbsorb(fileListStandardNonAbsorb):
    imageData = [fits.getdata(file) for file in fileListStandardNonAbsorb]
    print ("Total image data:", imageData,'\n')
    standardNonAbsorbImage = np.sum(imageData, axis=0)
    
    fits.writeto('Stacked_Image_Standard_NonAbsorb.fits', standardNonAbsorbImage, overwrite=True)
    
    print ("Image Standard Non-Absorb:", standardNonAbsorbImage,'\n')
    
    plt.clf()
    plt.imshow(standardNonAbsorbImage)
    plt.savefig('Stacked_Image_Standard_NonAbsorb', dpi=100)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    #Gaussian Blur
    plt.clf()
    betterImage = ndimage.gaussian_filter(standardNonAbsorbImage, sigma=(2,2), order=0)
    plt.imshow(betterImage)
    plt.savefig('Stacked_Image_Standard_NonAbsorb_Blur', dpi=100)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.colorbar()
    plt.show()
    
    print ("stackStandardNonAbsorb Function Complete!")
#End stackStandardNonAbsorb function
    
# =============================================================================
# Program's 'main'. Begin reading in all images for stacking.
# =============================================================================
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
stackMedianAbsorb(fileListMedianAbsorb)
stackMeanNonAbsorb(fileListMeanNonAbsorb)
stackMedianNonAbsorb(fileListMedianNonAbsorb)
stackStandardAbsorb(fileListStandardAbsorb)
stackStandardNonAbsorb(fileListStandardNonAbsorb)
