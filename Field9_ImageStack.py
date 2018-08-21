# -*- coding: utf-8 -*-
"""
Created on Tue May 29 21:05:05 2018

@author: Matthew Peek
Last Modified: 21 August 2018
Field 9 Image Stack
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
    
    print ("Normed:", normed)   
    print ("Normed Sum:", np.sum(normed))
    print ("Normalization complete!")     
    return normed
#End Absorber imageNorm function

def imageNormNonAbsorber(fileName):
    data = fits.getdata(fileName)
    sumData = np.sum(data)
    print ("Summed Image Data:", sumData,'\n')
    print ("Data:", data,'\n')
        
    normed = (data / sumData)
    
    print ("Normed:", normed)   
    print ("Normed Sum:", np.sum(normed))
    print ("Normalization complete!")     
    return normed
#End Non-Absorber imageNorm function


"""
alignImages function takes normed numpy array and image ID number and rotates
numpy array using rotate function. Returns rotated numpy image. Else, returns
error statement.
"""
shape = []
def alignImages(normed, ID):
    if (ID == '336'):
        rotImage = rotate(normed, 79.9759, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '331'):
        rotImage = rotate(normed, -68.5231, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '352'):
        rotImage = rotate(normed, 53.5262, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '333'):
        rotImage = rotate(normed, -81.5679, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '377'):
        rotImage = rotate(normed, 56.8418, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '360'):
        rotImage = rotate(normed, 70.2982, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '341'):
        rotImage = rotate(normed, -80.5149, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '410'):
        rotImage = rotate(normed, 23.2531, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '371'):
        rotImage = rotate(normed, -30.3901, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '454'):
        rotImage = rotate(normed, -7.6972, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '316'):
        rotImage = rotate(normed, -66.9743, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '396'):
        rotImage = rotate(normed, 45.1241, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '391'):
        rotImage = rotate(normed, -35.4654, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '207'):
        rotImage = rotate(normed, -30.6527, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '370'):
        rotImage = rotate(normed, -60.7307, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '416'):
        rotImage = rotate(normed, 55.9022, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '421'):
        rotImage = rotate(normed, -29.1317, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '363'):
        rotImage = rotate(normed, 61.1025, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '409'):
        rotImage = rotate(normed, -35.5607, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '439'):
        rotImage = rotate(normed, 63.9652, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '342'):
        rotImage = rotate(normed, 78.592, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '295'):
        rotImage = rotate(normed, 19.1239, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '247'):
        rotImage = rotate(normed, 24.3813, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '373'):
        rotImage = rotate(normed, -7.7152, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '322'):
        rotImage = rotate(normed, -70.9729, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '484'):
        rotImage = rotate(normed, 37.7005, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '383'):
        rotImage = rotate(normed, 52.4581, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '435'):
        rotImage = rotate(normed, -72.0903, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '355'):
        rotImage = rotate(normed, -58.5689, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '314'):
        rotImage = rotate(normed, -16.5677, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '434'):
        rotImage = rotate(normed, 51.9941, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '384'):
        rotImage = rotate(normed, -85.0748, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    
    else:    
        print ("Align Images Function Error!",'\n',"Image " + ID + " Not Found!")
#End alignAbsorber function
  
      
"""
Stack functions, takes numpy array as argument, loops through array argument
getting fits data and combining all image data into same array.

Stack image data standard, mean averaging, and median averaging.
Write results to new fits files.
"""
def stackAbsorber(fileListAbsorb):
    #imageData = [fits.getdata(file) for file in fileList]
    imageData = [file for file in fileListAbsorb]
    
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
    
    print ("Stacking Absorbers Complete!")
#End Absorber Stack function


def stackNonAbsorber(fileListNonAbsorb):
    #imageData = [fits.getdata(file) for file in fileList]
    imageData = [file for file in fileListNonAbsorb]
    
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
    
    print ("Stacking Non-Absorbers Complete!")
#End Non-Absorber Stack function


def stackAll(fileListAll):
    imageData = [file for file in fileListAll]
    
    print ("Total Non-Absorber image data:", imageData,'\n')
    meanImage = np.mean(imageData, axis=0)
    medianImage = np.median(imageData, axis=0)
    imageStack = np.sum(imageData, axis=0)
    
    fits.writeto('Field9_Stacked_Image_Mean_Normed_All.fits', meanImage, overwrite=True)
    fits.writeto('Field9_Stacked_Image_Median_Normed_All.fits', medianImage, overwrite=True)
    fits.writeto('Field9_Stacked_Image_Normed_All.fits', imageStack, overwrite=True)
    
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
#End stackAll function
    
#################################################################################
"""
Program's Main Begins Here.
"""    
#################################################################################
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
    #Exclude galaxy ID's 325, 371, & 377 due to grism over subtraction.
    if (ID[i] != '377'): 
        if (absorber[i] == 'Yes'):
            try:
                fileName = 'CROP-SDSS-J120639.85+025308.3-G141_00' + ID[i] + '.2d.fits'
                
                #Call imageNormAbsorber function.
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
                print (image.shape,'\n')
        
                countAbsorber += 1
            except IOError:
                print ("Image ID " + ID[i] + " not found!")

        else:
            try:
                fileName = 'CROP-SDSS-J120639.85+025308.3-G141_00' + ID[i] + '.2d.fits'
            
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
                print (image.shape,'\n')
        
                countNonAbsorber += 1
            except IOError:
                print ("Image ID " + ID[i] + " not found!")
        totalCount += 1    

#Combine absorber and non-absorber lists to stack all images
fileListAll = fileListAbsorb + fileListNonAbsorb

#Call stack functions
stackAbsorber(fileListAbsorb)
stackNonAbsorber(fileListNonAbsorb)
stackAll(fileListAll)
print ("Number of Absorbers Processed:", countAbsorber)
print ("Number of Non-Absorbers Processed:", countNonAbsorber)
print ("Total Number of Galaxies Processed:", totalCount)
