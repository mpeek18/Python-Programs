# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 18:13:29 2017

@author: Matthew Peek
Program to calculate galaxy angular distance from quasar
Field 2
Last Modified: 16 Sept. 2017
"""

import csv
import math
from math import radians
from astropy.io import ascii
import astropy.io.fits as fits
from astropy.table import Table

fileName = 'SDSS-J082946.90+185222.0.cat'

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
	#print ("DECSTRING", decstring)
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
print ("Quasar Coord", raDeg," ", decDeg)

# Loop through the main SDSS Catalog and write arrays of Galaxy ID's,
# RA, Dec, and Magnitude.
catRA = []
catDec = []
catID = []

ra1 = radians(raDeg) #RA Quasar
dec1 = radians(decDeg) #Dec Quasar

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
#print ("I'm here",catRA[0], catDec[0])
print (sdssCat)
print ("-------------------------------------------------------------------------", "\n")
print (catID)
print ("-------------------------------------------------------------------------", "\n")

#Calculate the angular distance of galaxy objects from the quasar.
angDist = []
array = []
for i in range(0, len(catID)):
    ra2 = radians(catRA[i])
    dec2 = radians(catDec[i])
    
    # Angular distance equation
    angSep = math.acos(math.sin(dec2)*math.sin(dec1) + math.cos(dec2) * math.cos(dec1) * math.cos(ra2-ra1))*(180/math.pi)*3600
    angDist.append(angSep)
    #if (catID[i] == 350):
        #print ("350",angDist[i])
print (angDist)
print ("CATID", len(catID))
print ("ANGDIST", len(angDist))
print ("-------------------------------------------------------------------------")

#Write galaxy ID, galaxy RA, and galaxy Dec to ascii table to read into RunGrismFiles program
RaDecData = Table([catID, catRA, catDec, angDist], 
            names=['Galaxy ID', 'Galaxy RA', 'Galaxy Dec', 'Angular Distance'])
ascii.write(RaDecData, 'RaDecData.dat', format='fixed_width', overwrite=True)
print ("RaDecData.dat file written")
"""
file = ascii.read('RaDecData.dat', delimiter="|")
print (file['col2'])
print (file['col3'])
print (file['col4'])
print (file['col5'])
"""
