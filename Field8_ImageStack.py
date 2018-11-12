# -*- coding: utf-8 -*-
"""
Created on Tue May 29 21:05:05 2018

@author: Matthew Peek
Last Modified: 11 November 2018
Field 8 Image Stack
"""
import numpy as np
from astropy.io import ascii
import astropy.io.fits as fits
from astropy.table import Table
from matplotlib import pyplot as plt
from skimage.transform import rotate, resize

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
getGalAngle function takes galaxy ID as argument, goes through field7IDs list 
and matches passed argument. Finds associated galaxy angle and returns angle.
"""
def getGalAngle(ID):
    for i in range(0, len(field7IDs)):
        if (ID == field7IDs[i]):
            galAngle = field7Angles[i]
    print ("ID " + ID + " " + "angle " + galAngle)
    return galAngle
#End getGalAngle Function
    
"""
alignImages function takes normed numpy array and image ID number and rotates
numpy array using rotate function. Returns rotated numpy image. Else, returns
error statement.
"""
shape = []
def alignImages(normed, ID, galAngle):
    rotImage = rotate(normed, float(galAngle), True)
    print ("Image " + ID + " Rotated!")
    print ("Image shape:", rotImage.shape)
    return rotImage  
#End alignImages function


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
Read in All_Galaxy_Angles file and get fields, galaxy ID's, and Angles.
Find only galaxies for field 7 and append ID's and Angles to new list
for processing in alignImages function.
"""
fields = []
galaxyIDs = []
angles = []
try:
    angleFile = open('All_Galaxy_Angles_M.txt', 'r')
    for line in angleFile:
        fields.append(line.split()[0])
        galaxyIDs.append(line.split()[2])
        angles.append(line.split()[4])
    angleFile.close()

except IOError:
    print ("File could not be found in current directory!")
    
print ("Fields", fields,'\n')
print ("Galaxy ID's", galaxyIDs,'\n')
print ("Angles", angles,'\n')   
print ("Fields", fields,'\n')
print ("Galaxy ID's", galaxyIDs,'\n')
print ("Angles", angles,'\n')   

field7IDs = []
field7Angles = []
for i in range(0, len(fields)):
    if (fields[i] == '7'):
        field7IDs.append(galaxyIDs[i])
        field7Angles.append(angles[i])
        
print ("Field 7 ID's", field7IDs,'\n')
print ("Field 7 Angles", field7Angles,'\n')

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
redshift = absorberFile['col3']
absorber = absorberFile['col7']
wavelength = absorberFile['col10']

totalCount = 0
countAbsorber = 0
countNonAbsorber = 0
haAbsorb = []
haNonAbsorb = []
objIDAbsorb = []
fileListAll = []
objIDNonAbsorb = []
fileListAbsorb = []
objRedshiftAbsorb = []
fileListNonAbsorb = []
objRedshiftNonAbsorb = []
for i in range(1, len(ID)):
    if (ID[i] != '308' and ID[i] != '480'):
        if (absorber[i] == 'Yes'):
            
            #Append absorber ID's, redshifts, and wavelength to lists for ascii
            #table output.
            haAbsorb.append(wavelength[i])
            objIDAbsorb.append(ID[i])
            objRedshiftAbsorb.append(redshift[i])
        
            try:
                fileName = 'CROP-SDSS-J120342.24+102831.8-G141_00' + ID[i] + '.2d.fits'
        
                #Call imageNorm function.
                normed = imageNormAbsorber(fileName)
                
                """
                Try and match ID's from Absorbtion data file with ID's from
                All Galaxy Angles file. If match found, call getGalAngles function
                and pass current matching ID as argument.
                """
                galAngle = float(0.0)
                if (ID[i] in field7IDs):
                    galAngle = getGalAngle(ID[i])
                #print ("ID " + ID[i] + " " + "Angle " + galAngle) 
                
                """
                Call alignImages function and resize to stack.
                Note, images are not all the same size after rotating them, must resize
                in order to stack images.
                """
                rotImage = alignImages(normed, ID[i], galAngle)
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
            
            #Append non-absorber ID's, redshifts, and wavelength to list for ascii
            #table output.
            haNonAbsorb.append(wavelength[i])
            objIDNonAbsorb.append(ID[i])
            objRedshiftNonAbsorb.append(redshift[i])
        
            try:
                fileName = 'CROP-SDSS-J120342.24+102831.8-G141_00' + ID[i] + '.2d.fits'
                
                #Call imageNormNonAbsorb function.
                normed = imageNormNonAbsorber(fileName)
                
                """
                Try and match ID's from Absorbtion data file with ID's from
                All Galaxy Angles file. If match found, call getGalAngles function
                and pass current matching ID as argument.
                """
                galAngle = float(0.0)
                if (ID[i] in field7IDs):
                    galAngle = getGalAngle(ID[i])
                #print ("ID " + ID[i] + " " + "Angle " + galAngle) 
                
                """
                Call alignImages function and resize to stack.
                Note, images are not all the same size after rotating them, must resize
                in order to stack images.
                """
                rotImage = alignImages(normed, ID[i], galAngle)
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

# =============================================================================
# Write data to ascii table
# =============================================================================
stackDataAbsorbers = (Table([objIDAbsorb, objRedshiftAbsorb, haAbsorb], 
                            names=['ID Absorber', 'Redshift Absorber', 'Wavelength']))
ascii.write(stackDataAbsorbers, 'Field8_Stack_Data_Absorb.dat', format='fixed_width', overwrite=True)

stackDataNonAbsorbers = (Table([objIDNonAbsorb, objRedshiftNonAbsorb, haNonAbsorb], 
                               names=['ID Non-Absorber', 'Redshift Non-Absorber', 'Wavelength']))
ascii.write(stackDataNonAbsorbers, 'Field8_Stack_Data_NonAbsorb.dat', format='fixed_width', overwrite=True)
print ("Field8_Stack_Data file has been written")
