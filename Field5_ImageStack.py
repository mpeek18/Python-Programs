# -*- coding: utf-8 -*-
"""
Created on Tue May 29 21:05:05 2018

@author: Matthew Peek
Last Modified: 5 August 2018
Field 5 Image Stack
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
    if (ID == '373'):
        rotImage = rotate(normed, 55.5931, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '295'):
        rotImage = rotate(normed, 4.3917, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '387'):
        rotImage = rotate(normed, -41.3679, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '208'):
        rotImage = rotate(normed, 32.8373, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '338'):
        rotImage = rotate(normed, -30.4187, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '361'):
        rotImage = rotate(normed, 88.24, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '329'):
        rotImage = rotate(normed, -64.8267, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '416'):
        rotImage = rotate(normed, -58.6607, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '317'):
        rotImage = rotate(normed, 57.9865, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '340'):
        rotImage = rotate(normed, 75.9299, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '264'):
        rotImage = rotate(normed, 24.8136, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '376'):
        rotImage = rotate(normed, 83.2816, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '198'):
        rotImage = rotate(normed, -7.2072, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '263'):
        rotImage = rotate(normed, -29.0901, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '354'):
        rotImage = rotate(normed, 21.2342, True)
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
    
    fits.writeto('Field5_Stacked_Image_Mean_Normed_Absorber.fits', meanImage, overwrite=True)
    fits.writeto('Field5_Stacked_Image_Median_Normed_Absorber.fits', medianImage, overwrite=True)
    fits.writeto('Field5_Stacked_Image_Normed_Absorber.fits', imageStack, overwrite=True)
    
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
    
    fits.writeto('Field5_Stacked_Image_Mean_Normed_NonAbsorber.fits', meanImage, overwrite=True)
    fits.writeto('Field5_Stacked_Image_Median_Normed_NonAbsorber.fits', medianImage, overwrite=True)
    fits.writeto('Field5_Stacked_Image_Normed_NonAbsorber.fits', imageStack, overwrite=True)
    
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
    
    fits.writeto('Field5_Stacked_Image_Mean_Normed_All.fits', meanImage, overwrite=True)
    fits.writeto('Field5_Stacked_Image_Median_Normed_All.fits', medianImage, overwrite=True)
    fits.writeto('Field5_Stacked_Image_Normed_All.fits', imageStack, overwrite=True)
    
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
    if (ID[i] != '214' and ID[i] != '393'):
        if (absorber[i] == 'Yes'):
            try:
                fileName = 'CROP-SDSS-J095432.63+354027.7-G141_00' + ID[i] + '.2d.fits'
                
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
                fileName = 'CROP-SDSS-J095432.63+354027.7-G141_00' + ID[i] + '.2d.fits'
                
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
