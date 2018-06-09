# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 22:07:40 2017

@author: Matthew Peek
Last Modified 14 April 2017

Algorithm:
    open grism files, get hdu info
    
def processFits()    
    open empty fits file
    clean up grism using hdu dimensions, 
    final = sci-contam-model
    iterate through final
    write finalData to newImage
    close empty fits file

def analyzeFits()    
    open newImage file
    enter equations to calculate H-Alpha
    get wavelength from fits dimension
    iterate through wavelength calculating H-Alpha data
    crop newImage using found H-Alpha data
    write cropImage to new fits image
    close newImage file
"""
import numpy as np
from astropy.io import ascii
import astropy.io.fits as fits
from astropy.table import Table
import scipy.ndimage as ndimage
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from photutils import CircularAnnulus
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from photutils import CircularAperture, aperture_photometry
###############################################################################
# Open Grism, read and process sci, contam, model. Write to new final fits file
def processFits(fitsName):
    hdulist = fits.open(fitsName)
    hdulist.info()
    print()
    sci, header = hdulist[5].data, hdulist[5].header
    sciHeight = header['NAXIS2']
    """
    plt.clf()
    plt.imshow(sci, cmap='gray', norm=LogNorm())
    plt.colorbar()
    plt.savefig('Sci')
    """
    contam, header = hdulist[8].data, hdulist[8].header
    """
    plt.clf()
    plt.imshow(contam, cmap='gray', norm=LogNorm())
    plt.colorbar()
    plt.savefig('contam')
    """
    model, header = hdulist[7].data, hdulist[7].header
    """
    plt.clf()
    plt.imshow(model, cmap='gray', norm=LogNorm())
    plt.colorbar()
    plt.savefig('model')
    """
    final = sci-contam-model
    """
    plt.clf()
    plt.imshow(final, cmap='gray', norm=LogNorm())
    plt.colorbar()
    plt.savefig('Final')
    """
    """
    #Gaussian Blur  
    imgBlur = ndimage.gaussian_filter(final, sigma=(2,2), order=0)
    plt.imshow(imgBlur, interpolation='nearest')
    plt.savefig('ImageBlur' + str(galID), dpi=300)
    plt.subplots_adjust(right=2.0)
    plt.show()
    """
###############################################################################
    # This section creates a blank fits file
    xAxis = header['NAXIS2']
    yAxis = header['NAXIS1']
    
    hdu = fits.PrimaryHDU()
    hdu.data = np.zeros((xAxis, yAxis))
    hdu.writeto('BlankFits.fits', clobber=True)
    print ("I finished making the new fits!")
    print ()
###############################################################################    
# Iterate through cleaned fits and write data into new fits file
    values = []
    f = fits.open('BlankFits.fits')
    newData = f[0].data
        
    for i in range(0, header['NAXIS2']):
        for j in range(0, header['NAXIS1']):
            finalImage = final[i,j]
            newData[i][j] = finalImage  #Use newData[i][j] for flux calculation
            values.append(newData[i][j])
    fits.writeto('Final-goodsn-14-G141_' + str(galID) + '.2D.fits', newData, f[0].header, clobber=True)
    f.close() 
    return sciHeight
#End processFits Function
###############################################################################
###############################################################################    
def analyzeFits(finalName, fitsName, zGal, sciHeight):
    hdulist = fits.open(fitsName)
    f = fits.open(finalName)
    newData = f[0].data
    
    xAxisMin = 0
    xAxisMax = 0
    wavelength, header = hdulist[9].data, hdulist[9].header
    hAlphaWavelength = 6563 #Angstroms
    OIIIWavelength = 5008 #Angstroms
    
    if (zGal < 1.6):    #If redshift is less then z=1.6, check H-Alpha lines
        hAlphaGal = (zGal+1) * hAlphaWavelength
        print ("H-Alpha:", hAlphaGal)
        print ()
    
        for m in range(0, len(wavelength)):
            if (abs(wavelength[m] - hAlphaGal) < 400):
                if (xAxisMin == 0):
                    xAxisMin = m
                xAxisMax = m
                print ("Wavelength: ", wavelength[m], " Index: ", m)
        print (xAxisMin)
        print (xAxisMax)
        print()
        newImage = newData[:,xAxisMin:xAxisMax]
        
    elif (zGal > 1.6):  #If redshift is greater than z=1.6, check O-III lines
        OIIIGal = (zGal+1) * OIIIWavelength
        print ("O-III:", OIIIGal)
        print ()
    
        for m in range(0, len(wavelength)):
            if (abs(wavelength[m] - OIIIGal) < 400):
                if (xAxisMin == 0):
                    xAxisMin = m
                xAxisMax = m
                print ("Wavelength: ", wavelength[m], " Index: ", m)
        print (xAxisMin)
        print (xAxisMax)
        print ()
        newImage = newData[:,xAxisMin:xAxisMax] #Wavelength data
        
    #Begin placing aperture and annulus
    xAxis = (xAxisMax - xAxisMin) / 2 #Center point X-axis
    yAxis = sciHeight / 2             #Center point Y-axis
    positions = [(xAxis, yAxis)]      #Center point plotted
    aperture = CircularAperture(positions, r=1)
    phot_table_ap = aperture_photometry(newImage, aperture)
    r_in = np.linspace(1, 11, 11)
    r_out = np.linspace(2, 12, 11)
    
    print ( "r_in =", r_in)
    print ("r_out =", r_out)
    print ()
    
    fluxArray = [phot_table_ap['aperture_sum'] / aperture.area()]
    radArray = [1]
    for i in range(0, len(r_in)):
        rIn = r_in[i]
        rOut = r_out[i]
        annulus_apertures = CircularAnnulus(positions, rIn, rOut)
        phot_table = aperture_photometry(newImage, annulus_apertures)
        fluxArray.append(phot_table['aperture_sum'] / annulus_apertures.area())
        rMean = (rOut + rIn) / 2
        radArray.append(rMean)
        
    print ("Radius:", radArray)
    print ("Flux:", fluxArray)
    print ()
    
    sumFlux = sum(fluxArray)
    print (sumFlux)
    
    #Plot flux as radius increases
    plt.clf()
    plt.plot(radArray, fluxArray)
    plt.show()
    
    fits.writeto('Crop-goodsn-14-G141_' + str(galID) + '.2D.fits', newImage, f[0].header, clobber=True)
    f.close()
    
    print ()
    print ("End of Fits Processing!")
    print ()
    """
    R90 Algorithm:
        fluxFraction = []
        percentSum = 0
        radiusCount = -1
        for f in fluxArray:
            percentSum += f
            fluxFraction = percentSum / SumFlux
            
            optional,
            if percentSum / sumFlux > .90:
                print radiusCount
                print rMean[radiusCount]
    """
    fluxFrac = []
    percentSum = 0
    for f in fluxArray:
        percentSum += f
        fluxFrac.append(percentSum / sumFlux)
    print ("FluxFrac:", fluxFrac)
    print ("Flux Fraction done!")
    print ("----------------------------------------------------------------------")
    print ("R50 and R90 computations")
    print ()
    
    #Find the radius at 50% enclosed
    radius50 = 1
    for i in range(0, len(fluxFrac)):
        if (fluxFrac[i] < .50):
            radius50 += 1       #Radius 50% enclosed
    print ("50 percent enclosed:", radius50)
    
    #Find the radius at 90% encolsed
    radius90 = 1
    for i in range(0, len(fluxFrac)):
        if (fluxFrac[i] < .90):
            radius90 += 1       #Radius 90% enclosed
    print ("90 percent enclosed:", radius90)
    print ()
    
    #Gaussian Blur  
    imgBlur = ndimage.gaussian_filter(newImage, sigma=(2,2), order=0)
    plt.imshow(imgBlur, interpolation='nearest')
    plt.savefig('ImageBlur' + str(galID), dpi=600)
    plt.subplots_adjust(right=2.0)
    plt.show()
    
    #Add circle to cropped image, xy= xAxis, yAxis, radius= radius50, radius90
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    circle50 = plt.Circle((xAxis, yAxis), radius=radius50, color='m', fill=False, lw=2)
    circle90 = plt.Circle((xAxis, yAxis), radius=radius90, color='r', fill=False, lw=2)
    plt.imshow(newImage, cmap='gray')
    ax.add_patch(circle50) #50% enclosed
    ax.add_patch(circle90) #90% enclosed
    plt.savefig('Crop-goodsn-14-G141_' + str(galID))
    plt.show()
    
    print ("End radius computation")
    print ("END OF FITS FILE")
    print ("----------------------------------------------------------------------")
    print ()
    return (radius50, radius90)
#End analyzeFits Function
###############################################################################
def starFormation(lineFitDat, zFile):
    datFile = ascii.read(lineFitDat)
    print (datFile)
    print ()
    line = datFile['line']
    flux = datFile['flux']
    
    #Loop through datFile, find flux for H-Alpha line
    for i in range(0, len(line)):
        if (line[i] == 'Ha'):
            hAlphaFlux = flux[i]
    print (hAlphaFlux,"\n")
    
    #Get galaxy redshift from z_peak_grism column
    zFits = ascii.read(zFile)
    redshifts = zFits['z_peak_grism']
    print (redshifts,"\n")
    
    #Find Luminosity Distance
    cosmo = FlatLambdaCDM(H0=68 * u.km / u.s / u.Mpc, Om0=0.3)
    lumDist = cosmo.luminosity_distance(redshifts)
    print ("Luminosity Distance =", lumDist)
###############################################################################
# Below this line call functions to process grism files.
###############################################################################
redshifts = [1.356, 1.748, 1.239]
galaxyID = [23662, 23527, 23654]

ID = []
redShifts = []
R50 = []
R90 = []
for i in range(0, len(galaxyID)): #Loop through ID #'s.
    galID = galaxyID[i]
    zGal = redshifts[i]
    fitsName = 'goodsn-14-G141_' + str(galID) + '.2D.fits'
    finalName = 'Final-' + fitsName
    print ("Name:",fitsName)
    print ("Redshift:", zGal)
    #call functions
    sciHeight = processFits(fitsName)
    radius50, radius90 = analyzeFits(finalName, fitsName, zGal, sciHeight)
    
    #Prepare items for ascii table
    ID.append(galID)
    redShifts.append(zGal)
    R50.append(radius50)
    R90.append(radius90)    
print ()
###############################################################################
#Call StarFormation function. Finds H-Alpha flux & Galaxy redshift
lineFitDat = 'goodsn-14-G141_23662.linefit.dat'
zFile = 'goodsn-14-G141_23662.new_zfit.dat'
starFormation(lineFitDat, zFile)
###############################################################################
#Write data to ascii table
data = Table([ID, redShifts, R50, R90], names=['Galaxy ID', 'Redshifts', 'R50', 'R90'])
ascii.write(data, 'HSTData.txt', format='fixed_width')
#file.close()
