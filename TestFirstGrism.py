# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 17:52:58 2017

@author: Matthew Peek
Purpose: Test first grism file
Last Modified: 20 Feb. 2017
"""

import astropy.io.fits as fits
import numpy as np
from matplotlib import pylab as plt
from matplotlib.colors import LogNorm
###############################################################################
# Read in origional fits file
grism1 = 'goodsn-14-G141_23662.2D.fits'
hdulist = fits.open(grism1)
hdulist.info()
image, header = hdulist[5].data, hdulist[5].header
#print (image,header)
print ('NAXIS1:', header['NAXIS1'])
print ('NAXIS2:', header['NAXIS2'])
print ('image shape:', image.shape)
###############################################################################
# Create new fits file
"""
hdulist = fits.PrimaryHDU()
hdulist.data = np.zeros((37,282))
hdulist.writeto('new.fits', clobber=True)
print (hdulist.header)
print ("I'm done!")
"""
###############################################################################
# Clean up grism final
plt.clf()
plt.imshow(image, cmap='gray')
#plt.subplots_adjust(right=2.0)
plt.colorbar()

weightMap = hdulist[6].data
plt.clf()
plt.imshow(weightMap, cmap='gray', norm=LogNorm())
#plt.subplots_adjust(right=2.0)
plt.colorbar()

contam = hdulist[8].data
plt.clf()
plt.imshow(contam, cmap='gray')
#plt.subplots_adjust(right=2.0)
plt.colorbar()

model, header = hdulist[7].data, hdulist[7].header

final = image-contam-model
plt.clf()
plt.imshow(final, cmap='gray')
#plt.subplots_adjust(right=2.0)
plt.colorbar()
plt.savefig('final')
print ()
###############################################################################
#Algorithm
#take pixels in row, for every pixel in row, take diffArray of pixel and
#the mean in row. flux - mean or flux / mean.
###############################################################################

# This is final clf ONLY. No mean calculations or difference
values = []
f = fits.open('newFits.fits')
newData = f[0].data

for i in range(0, header['NAXIS2']):
    for j in range(0, header['NAXIS1']):
        newValues = final[i,j]
        values.append(newValues)
        newData[i][j] = newValues
fits.writeto('new2.fits', newData, f[0].header, clobber=True)
f.close()

###############################################################################
"""
f = fits.open('newFits.fits')
newData = f[0].data
    
meanValues = [] #Declare empty array for median values
#pixelValues = []  #Declare empty array for pixel values
diffValues = []
#count = 0
for i in range(0, header['NAXIS2']):
    row = final[i, :]
    #print (len(row))
    rowMean = np.mean(final[i, :])
    #rowSTDEV = std(row) #Find standard deviation for each row
    for j in range(0, header['NAXIS1']):
        diff = (final[i ,j] - rowMean) 
        #if diff > np.std(final[i,:]):
        #if (diff > 0.02):
        diffValues.append(diff)
        newData[i][j] = diff
        
        #print (diffValues)
        #Write values to new fits file here
        #writeFits('new.fits', hdulist[0].data, diffValues)
        
    #print (row)
    #print ("Row Mean:", rowMean)
    meanValues.append(rowMean) #Write median values into array
    #pixelValues.append(i) #Write pixel values into array
    #count+=1
fits.writeto('new1.fits', newData, f[0].header, clobber=True)
f.close()
#print ("I'm done!")
#print ()
#print ('Mean Array', meanValues)
#print ("Count: ", count) 
#print ()
#print('Pixel Values:', pixelValues)
#print ()
#print ('diffArray Array:', diffValues)
#print ("Array size: ", len(diffValues))
#print ()
"""
###############################################################################
#Print out the wavelength
f = fits.open('new2.fits')
newData = f[0].data

minPix = 0
maxPix = 0
lamdaArray = []
wavelength, header = hdulist[9].data, hdulist[9].header
zGal = 1.35
hAlphaWavelength = 6563 #Angstroms
hAlphaGal = (zGal+1) * hAlphaWavelength
print (hAlphaGal)
print ()

for m in range(0, len(wavelength)):
    if (abs(wavelength[m] - hAlphaGal) < 400):
        #lamda = wavelength[m]
        #lamdaArray.append(lamda)
        if (minPix == 0):
            minPix = m
        maxPix = m
        #print (minPix, maxPix)
        print ("Wavelength: ", wavelength[m], " Index: ", m)
print (minPix)
print (maxPix)
newImage = newData[:,minPix:maxPix]
fits.writeto('Crop.fits', newImage, f[0].header, clobber=True)
f.close()
"""
plt.errorbar(pixelValues, meanValues, yerr=rowSTDEV) #Create standard deviation graph
plt.subplots_adjust(right=1.5)  #Adust size for graph/ makes easier to see
plt.show()
"""

###############################################################################
"""
# Create histogram for Mean Values
count = 0
meanArray = []
for i in range(0, len(meanValues)):
    meanArray.append(meanValues[i])
    count+=1
#print ("Mean Values: ", meanArray)
#print ('Count: ', count)
#print ()
plt.clf()
plt.hist(meanArray, bins=37, color='orange', label='Mean Values')
plt.subplots_adjust(right=2.0)
plt.legend()
plt.show()
###############################################################################
# Create histogram for difference values
print ()
count = 0
diffArray = []
for i in range(0, len(diffValues)):
    diffArray.append(diffValues[i])
    count+=1
#print ("diffArray Values: ", diffArray)
#print ("Count: ", count)
#print ()
plt.clf()
plt.hist(diffArray, bins=37, color='orange', label='diffArray Values')
plt.subplots_adjust(right=2.0)
plt.legend()
plt.show()
"""