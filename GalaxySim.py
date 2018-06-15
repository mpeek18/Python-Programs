# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 17:45:25 2018

@author: Matthew Peek
Last Modified: 15 June 2018
"""
import numpy as np
import astropy.io.fits as fits

"""
makeCopy function makes copies of given image. Takes numpy array containing
image data as a paramter and number of copies to make as a parameter
"""
def makeCopy(imageData, num):
    count = 0
    while (count < num):
        fits.writeto('Galaxy ' + str(count) + ' .fits', imageData, overwrite=True)
        count += 1
    print ("Image copies complete!")

 
"""
Read in single fits image, get image data and call appropriate function.
Throws error if fits image cannot be found.
"""       
try:
    fileName = 'CROP-SDSS-J120639.85+025308.3-G141_00360.2d.fits'
    image = fits.open(fileName)
    image = image[0].data
    print (image,'\n')
    
    #Call makeCopy function
    makeCopy(image, 5)
except IOError:
    print ("Image not found!")