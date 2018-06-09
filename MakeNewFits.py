# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 21:16:54 2017

@author: Matthew Peek
Purpose: Create new fits file
Last Modified: 16 Feb. 2017
"""

import astropy.io.fits as fits
import numpy as np

hdu = fits.PrimaryHDU()
hdu.data = np.zeros((39, 312))
print (hdu.header)
hdu.writeto('blankFits.fits', clobber=True)
print ("I'm done!")