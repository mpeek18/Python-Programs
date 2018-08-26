# -*- coding: utf-8 -*-
"""
Created on Tue May 29 21:05:05 2018

@author: Matthew Peek
Last Modified: 26 August 2018
Field 8 Image Stack
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
#End imageNormAbsorb function

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
#End imageNormNonAbsorber function


"""
alignImages function takes normed numpy array and image ID number and rotates
numpy array using rotate function. Returns rotated numpy image. Else, returns
error statement.
"""
shape = []
def alignImages(normed, ID):
    if (ID == '418'):
        rotImage = rotate(normed, 69.9865, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '395'):
        rotImage = rotate(normed, -43.6707, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '322'):
        rotImage = rotate(normed, 33.5314, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '419'):
        rotImage = rotate(normed, -11.1677, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '428'):
        rotImage = rotate(normed, -59.3438, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '363'):
        rotImage = rotate(normed, 5.9664, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '426'):
        rotImage = rotate(normed, -61.0051, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '332'):
        rotImage = rotate(normed, -87.7122, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '470'):
        rotImage = rotate(normed, -37.26, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '402'):
        rotImage = rotate(normed, -20.484, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '312'):
        rotImage = rotate(normed, 72.83, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '264'):
        rotImage = rotate(normed, 8.4771, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '382'):
        rotImage = rotate(normed, -15.0421, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '369'):
        rotImage = rotate(normed, -20.2335, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '425'):
        rotImage = rotate(normed, 10.8038, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '390'):
        rotImage = rotate(normed, -36.5522, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '314'):
        rotImage = rotate(normed, 78.3745, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '481'):
        rotImage = rotate(normed, 55.7708, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '443'):
        rotImage = rotate(normed, -85.6534, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '478'):
        rotImage = rotate(normed, -5.0887, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '488'):
        rotImage = rotate(normed, 74.8462, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '304'):
        rotImage = rotate(normed, 82.8102, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '277'):
        rotImage = rotate(normed, 36.9347, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '396'):
        rotImage = rotate(normed, 35.7361, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '461'):
        rotImage = rotate(normed, 61.88, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '288'):
        rotImage = rotate(normed, 14.4325, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '417'):
        rotImage = rotate(normed, -47.7795, True)
        print ("Image " + ID + " Rotated!")
        print ("Image shape:", rotImage.shape)
        return rotImage
    elif (ID == '528'):
        rotImage = rotate(normed, 32.6782, True)
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
    
    fits.writeto('Field8_Stacked_Image_Mean_Normed_Absorber.fits', meanImage, overwrite=True)
    fits.writeto('Field8_Stacked_Image_Median_Normed_Absorber.fits', medianImage, overwrite=True)
    fits.writeto('Field8_Stacked_Image_Normed_Absorber.fits', imageStack, overwrite=True)
    
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
    
    fits.writeto('Field8_Stacked_Image_Mean_Normed_NonAbsorber.fits', meanImage, overwrite=True)
    fits.writeto('Field8_Stacked_Image_Median_Normed_NonAbsorber.fits', medianImage, overwrite=True)
    fits.writeto('Field8_Stacked_Image_Normed_NonAbsorber.fits', imageStack, overwrite=True)
    
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
    
    fits.writeto('Field8_Stacked_Image_Mean_Normed_All.fits', meanImage, overwrite=True)
    fits.writeto('Field8_Stacked_Image_Median_Normed_All.fits', medianImage, overwrite=True)
    fits.writeto('Field8_Stacked_Image_Normed_All.fits', imageStack, overwrite=True)
    
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
    
#################################################################################
"""
Program's Main Begins Here.
"""    
#################################################################################
"""
Open Absorber_Data.dat file and get galaxy id's, get image file name, open fits data.
Start program by reading in id's and appending them to new list.
"""
absorberFile = ascii.read('Absorption_Data_Field7.dat', delimiter='|')
ID = absorberFile['col2']
absorber = absorberFile['col7']
totalCount = 0
countAbsorber = 0
countNonAbsorber = 0
fileListAll = []
fileListAbsorb = []
fileListNonAbsorb = []
for i in range(1, len(ID)):
    if (ID[i] != '308' and ID[i] != '480'):
        if (absorber[i] == 'Yes'):
            try:
                fileName = 'CROP-SDSS-J120342.24+102831.8-G141_00' + ID[i] + '.2d.fits'
        
                #Call imageNorm function.
                normed = imageNormAbsorber(fileName)
                
                """
                Call alignImages function and resize to stack.
                Note, images are not all the same size after rotating them, must resize
                in order to stack images.
                """
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
                fileName = 'CROP-SDSS-J120342.24+102831.8-G141_00' + ID[i] + '.2d.fits'
                
                #Call imageNormNonAbsorb function.
                normed = imageNormNonAbsorber(fileName)
                
                """
                Call alignImages function and resize to stack.
                Note, images are not all the same size after rotating them, must resize
                in order to stack images.
                """
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
