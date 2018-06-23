# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 17:45:25 2018

@author: Matthew Peek
Last Modified: 22 June 2018
Galaxy Simulator
"""
import numpy as np
from astropy.io import ascii
import astropy.io.fits as fits
from astropy.table import Table

"""
makeCopy function makes copies of given image. Takes numpy array containing
image data as a paramter and number of copies to make as a parameter
"""
ID = []
def makeCopy(imageData, num):
    count = 0
    while (count < num):
        fits.writeto('Galaxy_' + str(count) + '.fits', imageData, overwrite=True)
        ID.append(count)
        count += 1
    print ("Image copies complete!",'\n')
#End makeCopy function

"""
addNoise function creates random distribution of noise and adds it to an image.
Writes new, noisy image, to fits file.

Parameters take a fits image and the current count which image is being
processed.
"""
def addNoise(imageCopy, count):
    noise = np.random.randn(34, 34)
    imageNoise = imageCopy + noise
    fits.writeto('Galaxy_noise_' + str(count) + '.fits', imageNoise, overwrite=True)
    print ("Adding noise to image complete!", '\n')
    print ("Noise added to " + str(count + 1) + " images.", '\n')
#End addNoise function
    
################################################################################# 
"""
Read in single fits image, get image data and call appropriate function.
Throws error if fits image cannot be found.
"""       
try:
    fileName = 'CROP-SDSS-J120639.85+025308.3-G141_00360.2d.fits'
    image = fits.open(fileName)
    imageData = image[0].data
    image.close()
    print (imageData,'\n')
    
    #Call makeCopy function
    makeCopy(imageData, 500)
    
except IOError:
    print ("Image not found!")
    
#################################################################################    
"""
Read in copied galaxy images, pass images to addNoise function.
Keep count of which image is currently being processed and pass 
to addNoise function.
"""
galID = []
count = 0
for i in range(0, len(ID)):
    try:
        copyName = 'Galaxy_' + str(ID[i]) + '.fits'
        file = fits.open(copyName)
        imageCopy = file[0].data
        file.close()
        print (imageCopy,'\n')
        
        #Call addNoise function
        addNoise(imageCopy, count)
        galID.append(ID[i])
        count += 1
        
    except IOError:
        print ("Galaxy_" + str(ID[i]) + " could not be found!")


#Write galaxy ID's to ascii file
data = (Table([galID], names=['Galaxy ID']))
ascii.write(data, 'SimGalaxyID.dat', format='fixed_width', overwrite=True)