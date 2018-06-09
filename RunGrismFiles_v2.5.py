# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 22:07:40 2017

@author: Matthew Peek
Last Modified: 13 July 2017

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
import csv
import math
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
def processFits(fitsName, finalName):
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
    print ("Blank Fits Width:", xAxis)
    print ("Blank Fits Height:", yAxis)
    hdu = fits.PrimaryHDU()
    hdu.data = np.zeros((xAxis, yAxis))
    hdu.writeto('BlankFits.fits', overwrite=True)
    print ("I finished making the new fits!")
    print ()
###############################################################################    
# Iterate through cleaned fits and write data into new fits file
    values = []
    f = fits.open('BlankFits.fits')
    newData = f[0].data
        
    for x in range(0, header['NAXIS2']):
        for y in range(0, header['NAXIS1']):
            finalImage = final[x,y]
            newData[x][y] = finalImage  #Use newData[i][j] for flux calculation
            values.append(newData[x][y])
    #fits.writeto('FINAL-SDSS-J001453.19+091217.6-G141_00' + str(galID) + '.2D.fits', newData, f[0].header, overwrite=True)
    fits.writeto(finalName, newData, f[0].header, overwrite=True)
    f.close() 
    return sciHeight
#End processFits Function
###############################################################################
###############################################################################    
def analyzeFits(finalName, fitsName, zGal, sciHeight, hAlphaFlux):
    hdulist = fits.open(fitsName)
    f = fits.open(finalName)
    newData = f[0].data
    
    xAxisMin = 0
    xAxisMax = 0
    wavelength, header = hdulist[9].data, hdulist[9].header
    hAlphaWavelength = 6563 #Angstroms
    OIIIWavelength = 5008 #Angstroms
    print ("zGal:",zGal)
    #print ("Wavelength",wavelength)
    
    if (zGal < 1.6):    #If zGal is less then z=1.6, check H-Alpha lines
        hAlphaGal = hAlphaWavelength * (zGal+1)
        print ("H-Alpha:", hAlphaGal)
        print ()
    
        for m in range(0, len(wavelength)):
            if (abs(wavelength[m] - hAlphaGal) < 800):
                if (xAxisMin == 0):
                    xAxisMin = m
                xAxisMax = m
            print ("H-Alpha Wavelength: ", wavelength[m], " Index: ", m)
        print ("XMAX ", xAxisMin)
        print ("YMAX ", xAxisMax)
        print()
        newImage = newData[:,xAxisMin:xAxisMax]
    
    elif (zGal > 1.6):  #If zGal is greater than z=1.6, check O-III lines
        OIIIGal = (zGal+1) * OIIIWavelength
        print ("O-III:", OIIIGal)
        print ()
    
        for m in range(0, len(wavelength)):
            if (abs(wavelength[m] - OIIIGal) < 800):
                if (xAxisMin == 0):
                    xAxisMin = m
                xAxisMax = m
                print ("OIII Wavelength: ", wavelength[m], " Index: ", m)
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
    r_in = np.linspace(1, 11, 11)   #Inner radii per pixel
    r_out = np.linspace(2, 12, 11)  #Outer radii per pixel
    
    #print ( "r_in =", r_in)
    #print ("r_out =", r_out)
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
        
    #print ("Radius:", radArray)
    print ("Flux:", fluxArray)
    print ()
    
    sumFlux = sum(fluxArray)
    print (sumFlux)
    
    #Plot flux as radius increases
    plt.clf()
    plt.plot(radArray, fluxArray)
    plt.show()
    
    fits.writeto('CROP-SDSS-J001453.19+091217.6-G141_00' + str(galID) + '.2D.fits', newImage, f[0].header, overwrite=True)
    f.close()
    
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
    
    
    #Same as below but without gaussian blur
    #plt.clf()
    plt.imshow(newImage)
    plt.savefig('noImageBlur' + str(galID) + '.png')
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.show()
        
    #Gaussian Blur
    #plt.clf()
    imgBlur = ndimage.gaussian_filter(newImage, sigma=(2,2), order=0)
    plt.imshow(imgBlur)
    plt.savefig('ImageBlur' + str(galID) + '.png', dpi=100)
    plt.subplots_adjust(right=2.0)
    plt.subplots_adjust(top=1.0)
    plt.show()
    
    #Add circle to cropped image, xy= xAxis, yAxis, radius= radius50, radius90
    """
    #Flux and Redshift text added to image.
    #plt.clf()
    fig = plt.figure()
    left, width = 1, 10
    bottom, height = 1, 10
    right = left + width
    top = bottom + height
    galRedshift = 'Redshift:', zGal
    hFlux = 'H-Alpha Flux:', hAlphaFlux
    """
    
    #Gray scale galaxy cutout with gaussian blur and annulus
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    circle50 = plt.Circle((xAxis, yAxis), radius=radius50, color='m', fill=False, lw=2)
    circle90 = plt.Circle((xAxis, yAxis), radius=radius90, color='r', fill=False, lw=2)
    plt.imshow(imgBlur, cmap='gray')
    ax.add_patch(circle50) #50% enclosed
    ax.add_patch(circle90) #90% enclosed
    #ax.text(left, top, galRedshift, horizontalalignment='left', verticalalignment='top', fontsize=12, color='red')
    #ax.text(left, bottom, hFlux, horizontalalignment='left', verticalalignment='bottom', fontsize=14, color='green')
    plt.savefig('CROP-SDSS-J001453.19+091217.6-G141_00' + str(galID) + '.png', dpi=100)
    plt.show()
    
    """
    #Same image as above but with no gaussian blur
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    noBlurCircle50 = plt.Circle((xAxis, yAxis), radius=radius50, color='m', fill=False, lw=2)
    noBlurCircle90 = plt.Circle((xAxis, yAxis), radius=radius90, color='r', fill=False, lw=2)
    plt.imshow(newImage, cmap='gray')
    ax.add_patch(noBlurCircle50) #50% enclosed
    ax.add_patch(noBlurCircle90) #90% enclosed
    #ax.text(left, top, galRedshift, horizontalalignment='left', verticalalignment='top', fontsize=16, color='red')
    #ax.text(left, bottom, hFlux, horizontalalignment='left', verticalalignment='bottom', fontsize=16, color='green')
    plt.savefig('CROP-NB-SDSS-J001453.19+091217.6-G141_00' + str(galID) + '.png')
    plt.show()
    """
    
    print ("End radius computation")
    print ()
    return (radius50, radius90)
#End analyzeFits Function
###############################################################################
#Begin Function to calculate star formation rates
def starFormation(hAlphaFlux, zGal):
    """
    datFile = ascii.read(lineFitDat)
    print (datFile)
    print ()
    line = datFile['line']
    flux = datFile['flux']
    
    #Loop through datFile, find flux for H-Alpha line
    for i in range(0, len(line)):
        if (line[i] == 'Ha'):
            hAlphaFlux = flux[i] * math.pow(10, -17)
        elif (line[i] != 'Ha'):
            hAlphaFlux = 0
    """
    #Find Luminosity Distance
    mpc2cm = 3.085677581e+24
    cosmo = FlatLambdaCDM(H0=68 * u.km / u.s / u.Mpc, Om0=0.3)
    lumDist = cosmo.luminosity_distance(zGal)
    angDiameter = cosmo.angular_diameter_distance(zGal) #Angular Diameter Distance function, param(redshift)
    print ("Luminosity Distance =", lumDist)
    print ("Angular Diameter Distance =", angDiameter)
    
    #Convert angular diameter from MPC to Kpc, get radius, find area = pi * radius^2
    angDiameterKpc = angDiameter[0].value * 10**-3
    galRadius = angDiameterKpc / 2
    area = math.pi * math.pow(galRadius, 2)
    
    #Convert luminosity distance from Mpc to Cm and find star formation luminosity of galaxy
    lum = hAlphaFlux*4*math.pi*(lumDist*mpc2cm)**2.
    sfr = (7.9*10**-42)*lum
    sfr = sfr[0].value
    
    #Star formation surface density
    sfrSurfaceDens = sfr / area
    
    print("Star Formation Rate:", sfr)
    print ("Galaxy Area:", area)
    print ("Star Formation Surface Density:", sfrSurfaceDens,"\n")
    print ("END OF FITS FILE")
    print ("----------------------------------------------------------------------")
    print ("----------------------------------------------------------------------")
    
    return (sfr, area, sfrSurfaceDens)
    
#End starFormation Function
###############################################################################
# Below this line call functions to process grism files.
###############################################################################
#zGals = [1.28855, 1.43851, 1.45279]
#galaxyID = [236, 359]
#fluxHa = []
#z = []
"""
targetFile = ascii.read('Filtered_Objects.dat', data_start=1)
objID = targetFile['col2']
for i in range(0, len(objID)):
    galaxyID.append(objID[i])
"""
with open('Best_redux_Target1.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    galaxyID = []
    orient = []
    orientName = []
    galRedshift = []
    counter = -1
     
    for row in readCSV:
        counter += 1
        if (counter > 0): 
            galaxyID.append(row[0])
            orient.append(row[6])
            galRedshift.append(row[2])
            
            if (row[6] == '0'):
                orientName.append('-08-064-') #Orientation 0
                
            elif (row[6] == '1'):
                orientName.append('-17-084-') #Orientation 1
                
            elif (row[6] == '2'): 
                orientName.append('-') #Orientation 2
                    
print (galaxyID)
print (orient)
print (orientName)
#galaxyID = [303]
#orient = [2]
#orientName = ['-']
ID = []
redshift = []
starForm = []
galArea = []
sfrDens = []
R50 = []
R90 = []
for i in range(0, len(galaxyID)): #Loop through Galaxy ID #'s.    
    print ("BEGIN NEW FITS IMAGE PROCESS")
    print ()    
    
    galID = galaxyID[i]
    frame = orientName[i]
    
    #zfits.dat file
    zFits = ascii.read('SDSS-J001453.19+091217.6' + str(frame) + 'G141_' + str(galID).zfill(5) + '.zfit.dat')     
    zGal = zFits['z_max_spec']
    print (zGal)
    
    #Hubble fits image
    fitsName = 'SDSS-J001453.19+091217.6' + str(frame) + 'G141_' + str(galID).zfill(5) + '.2D.fits'
    finalName = 'FINAL-' + fitsName
    print ("Fits Name:",fitsName)
    print ("Final Name:", finalName)
    
    #linefit.dat file
    lineFitDat = 'SDSS-J001453.19+091217.6' + str(frame) + 'G141_' + str(galID).zfill(5) + '.linefit.dat'
    print (lineFitDat)
    datFile = ascii.read(lineFitDat)
    line = datFile['line']
    flux = datFile['flux']
    #Loop though lines and find H-Alpha flux if it exists.
    for i in range(0, len(line)):
        if (line[i] == 'Ha'):
            hAlphaFlux = flux[i] * math.pow(10, -17)
        elif (line[i] != 'Ha'):
            hAlphaFlux = 0
        
    #call functions
    sciHeight = processFits(fitsName, finalName)
    radius50, radius90 = analyzeFits(finalName, fitsName, zGal, sciHeight, hAlphaFlux)
    sfr, area, sfrSurfaceDens = starFormation(hAlphaFlux, zGal)
    
    #Prepare items for ascii table
    ID.append(galID)
    redshift.append(zGal)
    R50.append(radius50)
    R90.append(radius90)
    starForm.append(sfr)
    galArea.append(area)
    sfrDens.append(sfrSurfaceDens)
print ()
###############################################################################
    
###############################################################################
#Write data to ascii table
data = (Table([ID, redshift, R50, R90, starForm, galArea, sfrDens], 
       names=['Galaxy ID', 'zGals', 'R50', 'R90', 'Star Formation Rate', 'Galaxy Area', 'Sfr Surface Density']))

ascii.write(data, 'HSTData.dat', format='fixed_width', overwrite=True)
print ("HSTData.dat file has been written")
