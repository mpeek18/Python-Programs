# -*- coding: utf-8 -*-
"""
Created on Tue May 29 21:05:05 2018

@author: Matthew Peek
Last Modified: 5 August 2018
Field 4 Image Stack
"""
import numpy as np
from astropy.io import ascii
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
"""
def imageNormAbsorber(fileName):
    data = fits.getdata(fileName)
    sumData = np.sum(data)
    print ("Summed Image Data:", sumData,'\n')
    print ("Data:", data,'\n')
        
    normed = (data / sumData)
    """
    for i in range(0, len(data)):
        for j in range(0, len(data[0])):
            data[i][j] = (data[i][j] / sumData)
    print ("Loop Data:", data)
    print ("Loop Summed Data:", np.sum(data))
    """
    print ("Normed:", normed)   
    print ("Normed Sum:", np.sum(normed))
    print ("Normalization complete!")     
    return normed
#End imageNormAbsorb function

def imageNormNonAbsorber(fileName):
    data = fits.getdata(fileName)
    sumData = np.sum(data)
    print ("Summed Image Data:", sumData,'\n')
    print ("Data:", data,'\n')
        
    normed = (data / sumData)
    """
    for i in range(0, len(data)):
        for j in range(0, len(data[0])):
            data[i][j] = (data[i][j] / sumData)
    print ("Loop Data:", data)
    print ("Loop Summed Data:", np.sum(data))
    """
    print ("Normed:", normed)   
    print ("Normed Sum:", np.sum(normed))
    print ("Normalization complete!")     
    return normed
#End imageNormNonAbsorber function


"""
alignImages function takes normed numpy array and image ID number and rotates
numpy array using rotate function. Returns rotated numpy image. Else, returns
error statement.
"""
shape = []
def alignImages(normed, ID):
    if (ID == '427'):
        rotImage = rotate(normed, -45.7301, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '413'):
        rotImage = rotate(normed, -18.434, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '452'):
        rotImage = rotate(normed, -52.5584, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '428'):
        rotImage = rotate(normed, -6.5608, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '404'):
        rotImage = rotate(normed, 16.3298, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '417'):
        rotImage = rotate(normed, -59.7914, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '343'):
        rotImage = rotate(normed, -86.3324, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '461'):
        rotImage = rotate(normed, -82.6544, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '428'):
        rotImage = rotate(normed, 26.9037, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '380'):
        rotImage = rotate(normed, 37.7994, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '427'):
        rotImage = rotate(normed, -45.7301, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '413'):
        rotImage = rotate(normed, -18.434, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '452'):
        rotImage = rotate(normed, -52.5584, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '404'):
        rotImage = rotate(normed, 16.3298, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '417'):
        rotImage = rotate(normed, -59.7914, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '343'):
        rotImage = rotate(normed, -86.3324, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '461'):
        rotImage = rotate(normed, -82.6544, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '428'):
        rotImage = rotate(normed, 26.9037, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '380'):
        rotImage = rotate(normed, 37.7994, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    
    else:
        print ("Align Images Function Error!",'\n',"Image " + ID + " Not Found!")
#End alignImages Function


"""
Stack function, takes numpy array as argument, loops through array argument
getting fits data and combining all image data into same array.

Stack image data standard, mean averaging, and median averaging, compare results.
Write results to new fits files.
"""
def stackAbsorber(fileListAbsorb):
    imageData = [file for file in fileListAbsorb]
    
    print ("Total image data:", imageData,'\n')
    meanImage = np.mean(imageData, axis=0)
    medianImage = np.median(imageData, axis=0)
    imageStack = np.sum(imageData, axis=0)
    
    fits.writeto('Field4_Stacked_Image_Mean_Normed_Absorber.fits', meanImage, overwrite=True)
    fits.writeto('Field4_Stacked_Image_Median_Normed_Absorber.fits', medianImage, overwrite=True)
    fits.writeto('Field4_Stacked_Image_Normed_Absorber.fits', imageStack, overwrite=True)
    
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
    print ("Stacking Absorbers Complete!")
#End StackAbsorb function

def stackNonAbsorber(fileListNonAbsorb):
    imageData = [file for file in fileListNonAbsorb]
    
    print ("Total image data:", imageData,'\n')
    meanImage = np.mean(imageData, axis=0)
    medianImage = np.median(imageData, axis=0)
    imageStack = np.sum(imageData, axis=0)
    
    fits.writeto('Field4_Stacked_Image_Mean_Normed_NonAbsorber.fits', meanImage, overwrite=True)
    fits.writeto('Field4_Stacked_Image_Median_Normed_NonAbsorber.fits', medianImage, overwrite=True)
    fits.writeto('Field4_Stacked_Image_Normed_NonAbsorber.fits', imageStack, overwrite=True)
    
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
    print ("Stacking Non-Absorbers Complete!")
#End StackNonAbsorb function

def stackAll(fileListAll):
    imageData = [file for file in fileListAll]
    
    print ("Total image data:", imageData,'\n')
    meanImage = np.mean(imageData, axis=0)
    medianImage = np.median(imageData, axis=0)
    imageStack = np.sum(imageData, axis=0)
    
    fits.writeto('Field4_Stacked_Image_Mean_Normed_All.fits', meanImage, overwrite=True)
    fits.writeto('Field4_Stacked_Image_Median_Normed_All.fits', medianImage, overwrite=True)
    fits.writeto('Field4_Stacked_Image_Normed_All.fits', imageStack, overwrite=True)
    
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
    print ("Stacking All Complete!")
#End StackAll function
   
    
"""
Open Absorber_Data.dat file and get galaxy id's, get image file name, open fits data.
Start program by reading in id's and appending them to new list.
"""
absorberFile = ascii.read('Absorption_Data.dat', delimiter='|')
ID = absorberFile['col2']
absorber = absorberFile['col7']
totalCount = 0
countAbsorber = 0
countNonAbsorber = 0
fileListAll = []
fileListAbsorb = []
fileListNonAbsorb = []
for i in range(1, len(ID)):
    if (ID[i] != '445' and ID[i] != '548' and ID[i] != '360'):
        if (absorber[i] == 'Yes'):
            try:
                fileName = 'CROP-SDSS-J091730.18+324105.5-G141_00' + ID[i] + '.2d.fits'
        
                #Call imageNorm function.
                normed = imageNormAbsorber(fileName)
            
                #Call alignImages function and resize to stack.
                rotImage = alignImages(normed, ID[i])
                resized = resize(rotImage, (48,48))
        
                #Append normalized image to fileList to pass as argument to stack function.
                fileListAbsorb.append(resized)
                file = fits.open(fileName)
                image = file[0].data
                file.close()
        
                print (ID[i])
                print (image.shape)
        
                countAbsorber += 1
            except IOError:
                print ("Image ID " + ID[i] + " not found!")
            
        else:
            try:
                fileName = 'CROP-SDSS-J091730.18+324105.5-G141_00' + ID[i] + '.2d.fits'
                
                #Call imageNormNonAbsorb function.
                normed = imageNormNonAbsorber(fileName)
                
                #Call alignImages function and resize to stack.
                rotImage = alignImages(normed, ID[i])
                resized = resize(rotImage, (48,48))
        
                #Append normalized image to fileList to pass as argument to stack function.
                fileListNonAbsorb.append(resized)
                file = fits.open(fileName)
                image = file[0].data
                file.close()
        
                print (ID[i])
                print (image.shape)
            
                countNonAbsorber += 1
            except IOError:
                print ("Image ID " + ID[i] + " not found!")
        totalCount += 1

#Combine both lists to stack all images
fileListAll = fileListAbsorb + fileListNonAbsorb

#Call stack functions
stackAbsorber(fileListAbsorb)
stackNonAbsorber(fileListNonAbsorb)
stackAll(fileListAll)
print ("Number of Absorbers Processed:", countAbsorber)
print ("Number of Non-Absorbers Processed:", countNonAbsorber)
print ("Total Number of Galaxies Processed:", totalCount)
