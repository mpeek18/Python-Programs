# -*- coding: utf-8 -*-
"""
Created on Tue May 29 21:05:05 2018

@author: Matthew Peek
Last Modified: 27 July 2018
Galaxy Simulator Image Stack
"""
import timeit
import numpy as np
from astropy.io import ascii
import scipy.ndimage as image
import astropy.io.fits as fits
from matplotlib import pyplot as plt
from skimage.transform import rotate, resize

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

Imported timeit module to compare efficiency of one liner normalizing data and 
nested loops as normalizing data.
"""
def imageNorm(fileName):
    data = fits.getdata(fileName)
    sumData = np.sum(data)
    print ("Summed Image Data:", sumData,'\n')
    print ("Data:", data,'\n')
    
    startNormed = timeit.default_timer()
    normed = (data / sumData)
    print ("Normed:", normed)
    print ("Normed Sum:", np.sum(normed))
    stopNormed = timeit.default_timer()
    
    startLoop = timeit.default_timer()
    for i in range(0, len(data)):
        for j in range(0, len(data[0])):
            data[i][j] = (data[i][j] / sumData)
    print ("Loop Data:", data)
    print ("Loop Summed Data:", np.sum(data))
    stopLoop = timeit.default_timer()
    
    totalNormed = stopNormed - startNormed
    totalLoop = stopLoop - startLoop
    print ("Normed time elapsed:", totalNormed)
    print ("Loop time elapsed:", totalLoop)
    
    if (totalNormed < totalLoop):
        print ("One liner normed is faster!")
    else:
        print ("Nested loops are faster!")
    return normed
#End imageNorm function


"""
RotateImage function takes normed image and rotates by user defined degrees.
Returns rotated fits image.
"""
def rotateImage(dataList):
    rotated = image.rotate(dataList, 25.0)
    rotated = rotate(dataList, 25.0, True)
    print ("Image Rotated!")
    #derotate = image.rotate(rotated, -25.0)
    #rescaled = resize(derotate, (34,34))
    #print ("Image Derotated!")
    return rotated
#End rotateImage function
    

"""
Stack function, takes numpy array as argument, loops through array argument
getting fits data and combining all image data into same array.

Stack image data standard, mean averaging, and median averaging, compare results.
Write results to new fits files.
"""
def stack(fileList):
    #imageData = [fits.getdata(file) for file in fileList]
    imageData = [file for file in fileList]
    
    #print ("Total image data:", imageData,'\n')
    meanImage = np.mean(imageData, axis=0)
    medianImage = np.median(imageData, axis=0)
    imageStack = np.sum(imageData, axis=0)
    
    #fits.writeto('GalSim_Stacked_Image_Mean.fits', meanImage, overwrite=True)
    #fits.writeto('GalSim_Stacked_Image_Median.fits', medianImage, overwrite=True)
    #fits.writeto('GalSim_Stacked_Image.fits', imageStack, overwrite=True)
    
    fits.writeto('GalSim_Stacked_Image_Mean_Normed.fits', meanImage, overwrite=True)
    fits.writeto('GalSim_Stacked_Image_Median_Normed.fits', medianImage, overwrite=True)
    fits.writeto('GalSim_Stacked_Image_Normed.fits', imageStack, overwrite=True)
    
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
Open Absorber_Data.dat file and get galaxy id's, get image file name, open fits data.
Start program by reading in id's and appending them to new list.
"""
count = 0
absorberFile = ascii.read('SimGalaxyID.dat', delimiter='|')
ID = absorberFile['col2']

fileList = []
for i in range(1, len(ID)):
    try:
        fileName = 'Galaxy_noise_' + ID[i] + '.fits'
        
        #Call imageNorm function.
        dataList = imageNorm(fileName)
        #newImage = np.concatenate(dataList)
        #print ("New Image:", newImage)
        
        #Call rotateImage function.
        rotatedImage = rotateImage(dataList)
        
        #Append normalized image to fileList to pass as argument to stack function.
        fileList.append(rotatedImage)
        file = fits.open(fileName)
        fitsImage = file[0].data
        file.close()
        
        print (ID[i])
        print (fitsImage.shape)
        
        count += 1
    except IOError:
        print ("Image ID " + ID[i] + " not found!")
print ("Number of images processed:", count,'\n')

#Call Stack function
stack(fileList) 
