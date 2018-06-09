# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 00:23:59 2017

@author: Matthew Peek
Last modified: 10 March 2017
Purpose: Use this to get fits file info prior to using RunGrismFiles.py
"""
import astropy.io.fits as fits
###############################################################################
def fitsInfo(fileName):
    hdulist = fits.open(fileName)
    return hdulist.info()
###############################################################################
galaxyID = [23662]
for i in range(0, len(galaxyID)):
    galID = galaxyID[i]
    fileName = 'goodsn-14-G141-big_23662.new_zfit'
    fitsInfo(fileName)
    print ("End of Header:", fileName)
    print ()
