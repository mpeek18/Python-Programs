# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 22:07:40 2017

@author: Matthew Peek
Last Modified: 2 January 2018
Field 8

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
import matplotlib.mlab as mlab
from matplotlib import pyplot as plt
#from matplotlib.colors import LogNorm
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
    
        for m in range(0, len(wavelength)):
            if (abs(wavelength[m] - hAlphaGal) < 800):
                if (xAxisMin == 0):
                    xAxisMin = m
                xAxisMax = m
            #print ("H-Alpha Wavelength: ", wavelength[m], " Index: ", m)
        print ("XMAX ", xAxisMin)
        print ("YMAX ", xAxisMax)
        print()
        newImage = newData[:,xAxisMin:xAxisMax]
    
    elif (zGal > 1.6):  #If zGal is greater than z=1.6, check O-III lines
        OIIIGal = (zGal+1) * OIIIWavelength
        print ("O-III:", OIIIGal)
    
        for m in range(0, len(wavelength)):
            if (abs(wavelength[m] - OIIIGal) < 800):
                if (xAxisMin == 0):
                    xAxisMin = m
                xAxisMax = m
                #print ("OIII Wavelength: ", wavelength[m], " Index: ", m)
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
    r_in = np.linspace(1, 11, 110)   #Inner radii per pixel
    r_out = np.linspace(2, 12, 120)  #Outer radii per pixel
    
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
    #print ("Flux:", fluxArray)
    #print ()
    
    sumFlux = sum(fluxArray)
    #print (sumFlux)
    
    #Plot flux as radius increases
    plt.clf()
    plt.plot(radArray, fluxArray)
    plt.show()
    
    fits.writeto('CROP-SDSS-J120342.24+102831.8-G141_00' + str(galID) + '.2D.fits', newImage, f[0].header, overwrite=True)
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
    #print ("FluxFrac:", fluxFrac)
    print ("Flux Fraction done!")
    print ("----------------------------------------------------------------------")
    print ("R50 and R90 computations")
    print ()
    
    
    #Find the radius at 50% enclosed
    radius50 = 1.0
    for i in range(0, len(fluxFrac)):
        if (fluxFrac[i] < .50):
            radius50 += 1/10       #Radius 50% enclosed
    print ("50 percent enclosed:", radius50)
    
    #Find the radius at 90% encolsed
    radius90 = 1.0
    for i in range(0, len(fluxFrac)):
        if (fluxFrac[i] < .90):
            radius90 += 1/10       #Radius 90% enclosed
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
    #plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    circle50 = plt.Circle((xAxis, yAxis), radius=radius50, color='m', fill=False, lw=2)
    circle90 = plt.Circle((xAxis, yAxis), radius=radius90, color='r', fill=False, lw=2)
    plt.imshow(imgBlur, cmap='gray')
    ax.add_patch(circle50) #50% enclosed
    ax.add_patch(circle90) #90% enclosed
    #ax.text(left, top, galRedshift, horizontalalignment='left', verticalalignment='top', fontsize=12, color='red')
    #ax.text(left, bottom, hFlux, horizontalalignment='left', verticalalignment='bottom', fontsize=14, color='green')
    plt.savefig('CROP-J120342.24+102831.8-G141_00' + str(galID) + '.png', dpi=100)
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
    plt.savefig('CROP-NB-SDSS-J120342.24+102831.8-G141_00' + str(galID) + '.png')
    plt.show()
    """
    
    print ("End radius computation")
    print ()
    return (radius50, radius90)
#End analyzeFits Function
###############################################################################
#Begin Function to calculate star formation rates
def starFormation(hAlphaFlux, zGal, radius90):
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
    #Find Luminosity Distance & Galaxy Diameter
    mpc2cm = 3.085677581e+24
    cosmo = FlatLambdaCDM(H0=68 * u.km / u.s / u.Mpc, Om0=0.3)
    lumDist = cosmo.luminosity_distance(zGal)
    angDiameter = cosmo.angular_diameter_distance(zGal) #Angular Diameter Distance function, param(redshift)
    angDiameterDist = angDiameter[0].value
    print ("Luminosity Distance =", lumDist)
    lumDist = lumDist[0].value
    print ("Angular Diameter Distance =", angDiameter)

    #Hubble resolutin = .06 arcsec per pixel
    resolution = radius90 * .06
    
    #Convert resolution from arcseconds to degrees 
    degree = resolution / 3600
    
    #Convert resolution from degrees to radians
    radian = (degree * (2 * math.pi)) / 360
   
    #Find radius of galaxy = angular_diameter_distance * radius90, (in radians)
    #Convert from Mpc to Kpc
    #Find area of radius90 = pi * radius^2
    radiusCircle = angDiameterDist * radian
    radiusCircleKpc = radiusCircle * 10**3
    areaKpc = math.pi * radiusCircleKpc**2
    
    #Convert luminosity distance from Mpc to Cm and find star formation luminosity of galaxy
    lum = hAlphaFlux*4*math.pi*(lumDist*mpc2cm)**2.
    sfr = (7.9*10**-42)*lum
    #sfr = sfr[0].value
    
    #Star formation surface density
    sfrSurfaceDens = sfr / areaKpc
    
    print ("Star Formation Rate:", sfr)
    print ("Galaxy Area Kpc:", areaKpc)
    print ("Star Formation Surface Density:", sfrSurfaceDens,"\n")
    print ("END OF FITS FILE")
    print ("----------------------------------------------------------------------")
    print ("----------------------------------------------------------------------")
    
    return (sfr, areaKpc, sfrSurfaceDens, lumDist)
    
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
csvFile = open('Best_redux_Target7.csv')
with open('Best_redux_Target7.csv') as csvFile:
    readCSV = csv.reader(csvFile, delimiter=',')
    galaxyID = []
    orient = []
    orientName = []
    galRedshift = []
    zQual = []
    fluxList = []
    counter = -1
     
    for row in readCSV:
        counter += 1
        if (counter > 0): 
            galaxyID.append(row[0])
            fluxList.append(float(row[4]))
            orient.append(row[6])
            galRedshift.append(float(row[2]))
            zQual.append(row[7])
            
            if (row[6] == '0'):
                orientName.append('-02-122-') #Orientation 0
                
            elif (row[6] == '1'):
                orientName.append('-11-275-') #Orientation 1
                
            elif (row[6] == '2'): 
                orientName.append('-') #Orientation 2
csvFile.close()

#################################################################################
#Read quasar absorber file, read rows galaxy ID's, galaxy redshifts, and mgII
with open('Abs_data_Target7.csv') as csvFile:
    readCSV = csv.reader(csvFile, delimiter=',')
    galQsoID = []
    galZQso = []
    mgII = []
    count = -1
    
    for row in readCSV:
        count += 1
        if (count > 0):
            galZQso.append(row[0])
            mgII.append(row[1])
            galQsoID.append(row[16])
            
print (galZQso,'\n')
print (mgII,'\n')
print (galQsoID,'\n')
#################################################################################

#galaxyID = [303]
#orient = [2]
#orientName = ['-']
ID = []
redshift = []
starForm = []
galAreaKpc = []
sfrDens = []
distance = []
R50 = []
R90 = []
errorLineFile = []
errorZFile = []
countLineFileError = 0
countzFileError = 0
for i in range(0, len(galaxyID)): #Loop through Galaxy ID #'s.    
    print ("BEGIN NEW FITS IMAGE PROCESS")
    print ()    
    
    if (zQual[i] == 'likely' or zQual[i] == 'good' or zQual[i] == 'probable' and galRedshift[i] < 1.6):
        galID = galaxyID[i]
        frame = orientName[i]
        
        #zfits.dat file
        try:
            zFits = ascii.read('SDSS-J120342.24+102831.8' + str(frame) + 'G141_' + str(galID).zfill(5) + '.zfit.dat')     
            zGal = zFits['z_max_spec']
            print (zGal)
            
        except IOError:
            countzFileError += 1
            errorZFile.append(galID)
            """
            for i in range(0, len(galaxyID)):
                if (galaxyID[i] == galID):
                    print ("GalRedshift:", galRedshift[i])
                    zGal = galRedshift[i]
            """
            print ("zfit.dat file does not exist!")
            
            
        #Hubble fits image
        fitsName = 'SDSS-J120342.24+102831.8' + str(frame) + 'G141_' + str(galID).zfill(5) + '.2D.fits'
        finalName = 'FINAL-' + fitsName
        print ("Fits Name:",fitsName)
        print ("Final Name:", finalName)
        
        #linefit.dat file
        try:
            lineFitDat = 'SDSS-J120342.24+102831.8' + str(frame) + 'G141_' + str(galID).zfill(5) + '.linefit.dat'    
            print (lineFitDat)
            datFile = ascii.read(lineFitDat)
            line = datFile['line']
            flux = datFile['flux']
            #Loop though lines and find H-Alpha flux if it exists.
            for i in range(0, len(line)):
                if (line[i] == 'Ha'):
                    hAlphaFlux = flux[i] * math.pow(10, -17)
                elif (line[i] != 'Ha'):
                    hAlphaFlux = -1
                    
        except IOError:
            countLineFileError += 1
            errorLineFile.append(galID)
            for i in range(0, len(fluxList)):
                if (fluxList[i] >= 0):
                    hAlphaFlux = fluxList[i] * math.pow(10, -17)
                else:
                    hAlphaFlux = -1
            
        #call functions
        sciHeight = processFits(fitsName, finalName)
        radius50, radius90 = analyzeFits(finalName, fitsName, zGal, sciHeight, hAlphaFlux)
        sfr, areaKpc, sfrSurfaceDens, lumDist = starFormation(hAlphaFlux, zGal, radius90)
        
        #Prepare items for ascii table
        ID.append(galID)
        redshift.append(zGal)
        R50.append(radius50)
        R90.append(radius90)
        starForm.append(sfr)
        galAreaKpc.append(areaKpc)
        sfrDens.append(sfrSurfaceDens)
        distance.append(lumDist)
################################################################################
#List of galaxy ID's with associated absorption, from absorber csv file.
newGalID = ['322', '426', '332', '470', '418']

print ("newGalID array length:", len(newGalID))
print ("ID array length:", len(ID))

#Find how many sfrDens are zero and non-zero.
count = 0
countZero = 0
for i in range(0, len(sfrDens)):
    if (sfrDens[i] != 0):
        count +=1
    else:
        countZero += 1
print ("Non-zero sfrDens total:", count)
print ("Zero sfrDens total:", countZero)
print ("Total:", count + countZero)

###############################################################################

#Read in RaDecData.dat file and get galaxy id's and angular distance calculations
radecGal = [] #galaxy id's in RaDecData.dat file
angularDist = []

raDecFile = ascii.read('RaDecData.dat', delimiter="|")
galIDFile = raDecFile['col2']
angDist = raDecFile['col5']

#Write assigned columns to new lists
for i in range(0, len(galIDFile)):    
    radecGal.append(galIDFile[i])    
    angularDist.append(angDist[i])

###############################################################################

#Begin looping through sfrDensity array comparing values in ID list with galaxy ID absorber list.
#Find common galaxies that are associated with absorbers and galaxies NOT associated with absorbers.
#Create lists for both galaxy areas and distance for absorber criteria.
#Create new list, Yes if galaxy is absorber, No if not absorber
sfrDensAbsorb = []
galIDAbsorb = []
sfrDensNoAbsorb = []
galIDNoAbsorb = []
areaKpcAbsorb = []
areaKpcNoAbsorb = []
distAbsorb = []
distNoAbsorb = []
SFRAbsorb = []
SFRNoAbsorb = []
zDistAbsorb = []
zDistNonAbsorb = []
zQualAbsorb = []
zQualNonAbsorb = []

#New lists for total numbers to write to absorber ascii table
totalID = []
totalSFR = []
totalSFRDens = []
totalZQual = []
totalZDist = []
galAbsorption = []

for i in range(0, len(sfrDens)):
    if (sfrDens[i] >= 0):   #If sfr surface density is >= 0
        totalID.append(ID[i])
        totalSFR.append(starForm[i])
        totalSFRDens.append(float(sfrDens[i]))
        totalZQual.append(zQual[i])
        totalZDist.append(redshift[i])
        
        if (ID[i] in newGalID): #If galaxy ID in ID array is also in newGalID array
            sfrDensAbsorb.append(sfrDens[i])
            galIDAbsorb.append(ID[i])
            areaKpcAbsorb.append(galAreaKpc[i])
            distAbsorb.append(distance[i])
            SFRAbsorb.append(starForm[i])
            zDistAbsorb.append(redshift[i])
            zQualAbsorb.append(zQual[i])
            absorb = 'Yes'
            
        else:   #If galaxy ID in ID array is not also in newGalID array
            sfrDensNoAbsorb.append(sfrDens[i])
            galIDNoAbsorb.append(ID[i])
            areaKpcNoAbsorb.append(galAreaKpc[i])
            distNoAbsorb.append(distance[i])
            SFRNoAbsorb.append(starForm[i])
            zDistNonAbsorb.append(redshift[i])
            zQualNonAbsorb.append(zQual[i])
            absorb = 'No'
        
        galAbsorption.append(absorb)

print ("sfr density absorption:", len(sfrDensAbsorb))
print ("sfr density no absorption:", len(sfrDensNoAbsorb))
print ("sfr density array:", len(sfrDens))
print ()
print (sfrDensAbsorb,'\n')
print (sfrDensNoAbsorb,'\n')
print (galIDAbsorb,'\n')
print (galIDNoAbsorb)


totalAngDist = []
for i in range(0, len(totalID)):
    for j in range(0, len(radecGal)):
        if (totalID[i] == radecGal[j]):
            totalAngDist.append(float(angularDist[j]))
print (totalAngDist)
print (len(totalAngDist))
###############################################################################
#Non-existant files
print ("Non-existant linefit.dat files:", countLineFileError)
print ("Non-existant zFits.dat files:", countzFileError)
print (errorLineFile)
print (errorZFile)

###############################################################################

#Histogram plot star formation surface density for absorbers and non-absorbers.
#plt.clf()

#binArray local variable for histogram plots. Sets length and bin size.
binArray = np.linspace(0, 2, 10)

plt.hist(sfrDensAbsorb, bins=binArray, normed=True, histtype='step', label= 'MgII Detection (%i)' %len(sfrDensAbsorb))
plt.hist(sfrDensNoAbsorb, bins=binArray, normed=True, histtype='step', label = 'MgII Non-Detection (%i)' %len(sfrDensNoAbsorb))
plt.xlabel("SFR Surface Density $(M_{Sun} yr^{-1} Kpc^{-2})$")
plt.ylabel("Normalized Number")
plt.legend()
plt.savefig('Histogram_SfrDensity')
plt.show()

"""
SFR Units: M_Sun yr^-1
SFR Surface Density Untits: M_Sun yr^-1 KPC^-2
"""

#Histogram plot star formation rate for absorbers and non-absorbers.
#plt.clf()
fig, ax = plt.subplots()
ax.hist(sfrDensAbsorb, bins=binArray, normed=True, histtype='step', label= 'MgII Detection (%i)' %len(sfrDensAbsorb))
ax.hist(sfrDensNoAbsorb, bins=binArray, normed=True, histtype='step', label = 'MgII Non-Detection (%i)' %len(sfrDensNoAbsorb))
ax.set_ylim([0, 8])
plt.xlabel("SFR Surface Density $(M_{Sun} yr^{-1} Kpc^{-2})$")
plt.ylabel("Normalized Number")
plt.legend()
plt.savefig('Histogram_SfrDensity_ScaleY')
plt.show()

#Histogram plot star formation rates for absorbers and non-absorbers.
#plt.clf()
plt.hist(SFRAbsorb, bins=binArray, normed=True, histtype='step', label= 'MgII Detection (%i)' % len(SFRAbsorb))
plt.hist(SFRNoAbsorb, bins=binArray, normed=True, histtype='step', label = 'MgII Non-Detection (%i)' % len(SFRNoAbsorb))
plt.xlabel("SFR $(M_{Sun} yr^{-1})$")
plt.ylabel("Normalized Number")
plt.legend()
plt.savefig('Histogram_Sfr')
plt.show()

#Scatter plot galaxy area vs. star formation rate, absorber and non-absorber
#plt.clf()
plt.scatter(areaKpcAbsorb, SFRAbsorb, c='blue', label='Absorbers')
plt.scatter(areaKpcNoAbsorb, SFRNoAbsorb, c='r', marker='^', label='Non-Absorbers')
plt.xlabel("Galaxy Area (Kpc)")
plt.ylabel("Star Formation Rate")
plt.legend(loc='upper right')
plt.savefig('Absorber_Area_vs_SFR')
plt.show()

#Scatter plot galaxy area vs. star formation rate surface density, absorber and non-absorber
#plt.clf()
plt.scatter(areaKpcAbsorb, sfrDensAbsorb, c='blue', label='Absorbers')
plt.scatter(areaKpcNoAbsorb, sfrDensNoAbsorb, c='r', marker='^', label='Non-Absorbers')
plt.xlabel("Galaxy Area (Kpc)")
plt.ylabel("Star Formation Surface Density")
plt.legend(loc='upper right')
plt.savefig('Absorber_Area_vs_sfrDensity')
plt.show()

#Scatter plot galaxy distance vs. star formation rate, absorber and non-absorber
#plt.clf()
plt.scatter(distAbsorb, SFRAbsorb, c='blue', label='Absorbers')
plt.scatter(distNoAbsorb, SFRNoAbsorb, c='r', marker='^', label='Non-Absorbers')
plt.xlabel("Distance(Mpc)")
plt.ylabel("Star Formation Rate")
plt.legend(loc='upper right')
plt.savefig('Absorber_Distance_vs_SFR')
plt.show()

#Scatter plot galaxy distance vs. star formation rate surface density, absorber and non-absorber
#plt.clf()
plt.scatter(distAbsorb, sfrDensAbsorb, c='blue', label='Absorbers')
plt.scatter(distNoAbsorb, sfrDensNoAbsorb, c='r', marker='^', label='Non-Absorbers')
plt.xlabel("Distance (Mpc)")
plt.ylabel("Star Formation Surface Density")
plt.legend(loc='upper right')
plt.savefig('Absorber_Distance_vs_sfrDensity')
plt.show()

#Scatter plot galaxy distance vs. area, absorber and non-absorber
#plt.clf()
plt.scatter(distAbsorb, areaKpcAbsorb, c='blue', label='Absorbers')
plt.scatter(distNoAbsorb, areaKpcNoAbsorb, c='r', marker='^', label='Non-Absorbers')
plt.xlabel("Distance (Mpc)")
plt.ylabel("Galaxy Area (Kpc)")
plt.legend(loc='upper left')
plt.savefig('Absorber_Distance_vs_Area')
plt.show()

#Scatter plot total galaxy angular distance vs. total sfr densities. Line of best fit included.
plt.scatter(totalAngDist, totalSFRDens)
plt.plot(np.unique(totalAngDist), np.poly1d(np.polyfit((totalAngDist), totalSFRDens, 3))(np.unique(totalAngDist)), c='r')
plt.xlabel('Angular Distance (arcsec)')
plt.ylabel('Star Formation Rate Density')
plt.savefig('TotalSFRDens_vs_AngularDist')
plt.show()
print ()
###############################################################################
    
###############################################################################
#Write data to ascii table

#Table for combining sfrDensity data.
sfrAbsorbData = (Table([sfrDensAbsorb], names=['SFR Density Absorber']))
ascii.write(sfrAbsorbData, 'Field8_SFRDensAbsorb.dat', format='fixed_width', overwrite=True)

#Table for combining sfrDensity data.
sfrNoAbsorbData = (Table([sfrDensNoAbsorb], names=['SFR Density No Absorber']))
ascii.write(sfrNoAbsorbData, 'Field8_NoSFRDensAbsorb.dat', format='fixed_width', overwrite=True)

#Table for combining sfr data.
sfrAbsorbData = (Table([SFRAbsorb], names=['SFR Absorber']))
ascii.write(sfrAbsorbData, 'Field8_SFRAbsorb.dat', format='fixed_width', overwrite=True)

#Table for combining sfr data.
sfrNoAbsorbData = (Table([SFRNoAbsorb], names=['SFR Non Absorber']))
ascii.write(sfrNoAbsorbData, 'Field8_SFRNoAbsorb.dat', format='fixed_width', overwrite=True)

#Following tables are general output data of processed files.
data = (Table([ID, redshift, R50, R90, starForm, galAreaKpc, sfrDens], 
       names=['Galaxy ID', 'zGals', 'R50', 'R90', 'Star Formation Rate',
              'Galaxy Area (Kpc)', 'Sfr Surface Density']))
ascii.write(data, 'HSTData.dat', format='fixed_width', overwrite=True)
print ("HSTData.dat file has been written")

#Write absorber data to ascii table
absorberData = (Table([totalID, totalZDist, totalZQual, totalSFR, totalSFRDens, galAbsorption, totalAngDist],
                names=['Galaxy ID', 'Z Dist', 'Z Qual', 'Star Formation Rate',
                       'SFR Surface Density', 'Absorber', 'Angular Distance']))
ascii.write(absorberData, 'Absorption_Data.dat', format='fixed_width', overwrite=True)
print ("Absorption_Data.dat file has been written")
