# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 13:32:18 2017

@author: Matthew Peek
Last Modified: 14 June 2017
"""

import math
from astropy.io import ascii
import astropy.io.fits as fits
from astropy.table import Table

fileName = 'SDSS-J001453.19+091217.6.cat'

# Begin function to convert hours, minutes, seconds coordinates
# to degrees.
def hms2deg(rahex,dechex):
	
	rastring = str(rahex)
	rah = float(rastring[:2])  # first two characters
	ra1 = rastring[2:] # all but the first two characters
	ram = float(ra1[:2])
	ras = float(ra1[2:])		

	decstring = str(dechex)
	decs = float(decstring[5:])
	decm1 = decstring[:5]
	decm = float(decm1[3:])
	decd = float(decstring[:3])		
	decsign = str(decstring[:1])
	print ("DECSTRING", decstring)
	if (decsign == '-'):
		decd = -1.0 * float(decd)
	else:
		if (decsign != '+'):
			print ("format of hexigesimal coords must have +/- sign input for declination")
		
		
	ra = 15.0 * (float(rah) + (float(ram) + float(ras)/60.0)/60.0)
	dec = float(decd) + (float(decm) + float(decs)/60.0)/60.0

	if (decsign == '-') and dec>0.:
		dec *= -1
			
	radeg=(float(ra))
	decdeg=(float(dec))
		
	return (radeg, decdeg)

# End hms2deg conversion function
#################################################################################

sdssCat = ascii.read(fileName)
coords = fileName[6:-4].split('+')
print (coords)
raDeg, decDeg = hms2deg(coords[0], '+' + coords[1])
print (raDeg," ", decDeg)

# Loop through the main SDSS Catalog and write arrays of Galaxy ID's,
# RA, Dec, and Magnitude.
catRA = []
catDec = []
catID = []
ra1 = raDeg #RA Quasar
dec1 = decDeg #Dec Quasar

galID = sdssCat['NUMBER']
galRA = sdssCat['X_WORLD']
galDec = sdssCat['Y_WORLD']
magnitude = sdssCat['MAG_AUTO']
for i in range(0, len(sdssCat)):
    gals = sdssCat[i]
    ID = galID[i]
    ra2 = galRA[i]
    dec2 = galDec[i]
    
    catID.append(ID)
    catRA.append(ra2)
    catDec.append(dec2)
print (sdssCat)
print ("-------------------------------------------------------------------------", "\n")
print (catID)
print ("-------------------------------------------------------------------------", "\n")

# Calculate the angular distance of galaxy objects from the quasar.
angDist = []
array = []
for i in range(0, len(catID)):
    ra2 = catRA[i]
    dec2 = catDec[i]
    
    # Angular distance equation
    a = math.sin(abs(dec1 - dec2) / 2)**2 #sin(|dec1 - dec2| / 2)^2
    b = math.cos(dec1) * math.cos(dec2) * math.sin(abs(ra1 - ra2) / 2)**2 #cos(dec1) * sin(|ra1 - ra2| / 2)^2
    d = 2 * math.asin(math.sqrt(a + b)) #2 * arcsin(sqrt(a + b))
    angDist.append(d)
print (angDist)
galDist = zip(catID, angDist)
print ("-------------------------------------------------------------------------", "\n")
print (list(galDist))
print ("-------------------------------------------------------------------------")
print ()

Dist = []
galaxyID = []
mag = []
for i in range(0, len(angDist)):
    arcSec = (angDist[i]) * 3600
    if (arcSec < 30):
        Dist.append(arcSec)
        galaxyID.append(catID[i])
        mag.append(magnitude[i])
print (Dist)
print (len(Dist))
print ()
print (galaxyID, "\n")
print (mag)
print ("-------------------------------------------------------------------------")
print ()


# Iterate through galaxy ID's find zfit.dat files that DO exist,
# write new ascii table with cooresponding Galaxy ID's, Magnitudes,
# Distance, and Redshifts.
redshifts = []
filterGals = []
galDist = []
galMag = []
fluxHa = []
for i in range(0, len(galaxyID)):
    ID = galaxyID[i]
    
    try:
        zFits = ascii.read('SDSS-J001453.19+091217.6-G141_00' + str(ID) + '.zfit.dat')
        lineFit = ascii.read('SDSS-J001453.19+091217.6-G141_00' + str(ID) + '.linefit.dat')
        redshift = zFits['z_peak_spec']
        line = lineFit['line']
        flux = lineFit['flux']
        
        for i in range(0, len(line)):
            if (line[i] == 'Ha' and flux[i] > 0.00):
                if (redshift > 0.7 and redshift < 1.6): #Check if redshift = 0.7 < Z < 1.6
                    if (mag[i] > 18.0 and mag[i] < 26):
                        redshifts.append(redshift)
                        filterGals.append(ID)
                        galDist.append(Dist[i])
                        galMag.append(mag[i])
                        fluxHa.append(flux[i])
        
    except IOError:
        print ("Galaxy ID: " + str(ID) + " zfit.dat and/or linefit.dat file does not exist!")

print (filterGals)
print ()        
print (redshifts)
print ()
print (fluxHa)
print ()
print (len(fluxHa))
print (len(filterGals))
print (len(galDist))
print (len(galMag))
print (len(redshifts))
print ("-------------------------------------------------------------------------")
print ()


# Iterate through galaxy ID's and check for linefit.dat files that DO exist
# for the given ID.
# Placed this check inside above try/catch in order to write prospective 
# targets to ascii file.
lineID = []
for i in range(0, len(galaxyID)):
    newID = galaxyID[i]
    
    try:
        linefit = ascii.read('SDSS-J001453.19+091217.6-G141_00' + str(newID) + '.linefit.dat')
        lineID.append(newID)
        
    except IOError:
        print ("Galaxy ID: " + str(newID) + " linefit.dat file does not exist!")
        
print ()
print (lineID, "\n")

# Write the filtered galaxy data to ascii table for,
# processing in "RunGrismFiles" program.
redshiftData = Table([filterGals, galMag, galDist, redshifts, fluxHa], names=['ID', 'Magnitude', 'Distance', 'Redshift', 'H-Alpha'])
ascii.write(redshiftData, 'Filtered_Objects.dat', format='fixed_width', overwrite=True)

newFile = ascii.read('Filtered_Objects.dat', data_start=1)
objID = newFile['col2']
zShift = newFile['col5']
fluxHa = newFile['col6']
print (objID)
print (zShift)
print (fluxHa)
