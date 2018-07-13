# -*- coding: utf-8 -*-
"""
Created on Tue May 29 21:05:05 2018

@author: Matthew Peek
Last Modified: 4 July 2018
Galaxy Simulator Image Stack
"""
import numpy as np
from astropy.io import ascii
import astropy.io.fits as fits
from matplotlib import pyplot as plt
from sklearn.preprocessing import normalize

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
#Try making new copy of image and normalizing it. 
#Possible problem with changing pixel values but never writing new data
#back into normalized image.
def imageNorm(fileName):
    dataList = []
    #sumCheck = 0.0
    data = fits.getdata(fileName)
    dataList.append(data)
    sumData = np.sum(data)
    print ("Summed Image Data:", sumData,'\n')
    print ("Data:", data,'\n')
    
    data2 = (data / sumData)
    """
    for i in range(0, len(data)):
        for j in range(0, len(data[0])):
            #("Before:", data[i][j])
            sumCheck += data[i][j]
            data[i][j] = (data[i][j] / sumData)
            #data = data[i][j]
    print ("Sum Check:", sumCheck)       
    #print ("Loop Summed Data:", np.sum(data[i][j]))
            
    #print ("Normed:", data)  
    normData = np.sum(data[i][j])
    print ("Normed Sum:", normData)
    #print ("Normalization complete!")     
    #return data
    """
    print ("Data 2", data2)
    print ("Summed data 2:", np.sum(data2))
    data3 = data / data2
    print ("Data3:", data3)
    #normed = normalize(data, axis=1, norm='l1')
    #print ("Normed Sum:", np.sum(normed),'\n')
    #print ("Normed Data:", normed)
    return data2
#End imageNorm function


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
        
        #Append normalized image to fileList to pass as argument to stack function.
        fileList.append(dataList)
        file = fits.open(fileName)
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
