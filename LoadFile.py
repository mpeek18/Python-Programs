# -*- coding: utf-8 -*-
"""
Created on Fri May 26 13:35:49 2017

@author: Matthew Peek
Last Modified: 28 May 2017
"""

import math
from astropy.io import ascii
import astropy.io.fits as fits
from astropy.table import Table

sdssCat = ascii.read('SDSS-J001453.19+091217.6.cat')

catRA = []
catDec = []
catID = []
ra1 = 001453.19 #RA Quasar
dec1 = 091217.6 #Dec Quasar
galID = sdssCat['NUMBER']
galRA = sdssCat['X_WORLD']
galDec = sdssCat['Y_WORLD']
for i in range(0, len(sdssCat)):
    gals = sdssCat[i]
    ID = galID[i]
    ra2 = galRA[i]
    dec2 = galDec[i]
    
    catID.append(ID)
    catRA.append(ra2)
    catDec.append(dec2)
print (sdssCat)
print ()
print (catID)
print ()

angDist = []
array = []
for i in range(0, len(catID)):
    ra2 = catRA[i]
    dec2 = catDec[i]
    
    a = math.sin(abs(dec1 - dec2) / 2)**2
    b = math.cos(dec1) * math.cos(dec2) * math.sin(abs(ra1 - ra2) / 2)**2
    d = 2 * math.asin(math.sqrt(a + b))
    angDist.append(d)
print (angDist)
ID_AngDist = zip(catID, angDist)
print ()
print (list(ID_AngDist))
print ("-------------------------------------------------------------------------")
print ()


#Get RA and Dec
sdssRaDec = ascii.read('SDSS-J001453.19+091217.6.radec')
ra = sdssRaDec['col1']
dec = sdssRaDec['col2']

RA = []
Dec = []
for i in range(0, len(ra)):
    RA.append(ra[i])
    
for i in range(0, len(dec)):
    Dec.append(dec[i])
print ("RA DEC", "\n")
print (sdssRaDec)
print ("-------------------------------------------------------------------------")
print ()

RADist = []
sdssDist = ascii.read('SDSS-J001453.19+091217.6_dist.dat')
ID = sdssDist['col1']
ra = sdssDist['col4']

for i in range(1, len(ra)):
    if (ra[i] == 3.7216885519):
        galID = i 
print ("Index:", galID)
print (sdssDist)
print ("-------------------------------------------------------------------------")
print ()

galArray = []
RA = []
DEC = []
sdssSortCat = ascii.read('SDSS-J001453.19+091217.6_sorted.cat')
galID = sdssSortCat['col1']
galRA = sdssSortCat['col4']
galDec = sdssSortCat['col5']
for i in range(63, len(galID)):
    sortGals = galID[i]
    sortRA = galRA[i]
    sortDec = galDec[i]
    galArray.append(sortGals)
    RA.append(sortRA)
    Dec.append(sortDec)
print ("Sorted Cat: ", galArray, RA, Dec)
print ()
print (len(galArray))
print ("-------------------------------------------------------------------------")
print ()

"""
ANGULAR DISTANCE
d = 2 * arcsin(sqrt(a + b))
a = sin(abs(dec1 - dec2) / 2)**2
b = cos(dec1) * cos(dec2) * sin(abs(ra1 - ra2) / 2)**2
"""

#Calculate angular distance of galaxy ID's
result = []
ra1 = 001453.19 #RA Quasar
dec1 = 091217.6 #Dec Quasar
for i in range(0, len(galArray)):
    ra2 = RA[i]
    dec2 = Dec[i]
    a = math.sin(abs(dec1 - dec2) / 2)**2
    b = math.cos(dec1) * math.cos(dec2) * math.sin(abs(ra1 - ra2) / 2)**2
    d = 2 * math.asin(math.sqrt(a + b))
    result.append(d)
print ("Galaxy ID, Angular Distance ", galArray, result, "\n")

result.sort()
print ("Sorted Distance: ", result)
print ()

zipped = zip(galArray, result)
print (list(zipped))

data = Table([galArray, result], names=['Galaxy ID', 'Distance'])
ascii.write(data, 'Test.dat', format='fixed_width', overwrite=True)

#hms2dec = raH + (raM / 60) + (raS / (60 * 60))

"""
# Write your hms2dec and dms2dec functions here
import numpy as np
def hms2dec(hour, minute, second):
  ra2dec = 15*(hour + minute/60 + second/(60*60))
  return ra2dec

def dms2dec(hour, minute, second):
  if (hour < 0):
    hour = hour * -1
    dec2dec = -1*(hour + minute/60 + second/(60*60))
  else:
    dec2dec = 15*(hour + minute/60 + second/(60*60))
  return dec2dec


a = np.sin(np.abs(dec1Rad - dec2Rad)/2)**2
b = np.cos(dec1Rad)*np.cos(dec2Rad)*np.sin(np.abs(ra1 - ra2)/2)**2
d = 2*np.arcsin(np.sqrt(a + b))
"""