# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 10:18:25 2018

@author: Matthew Peek
Last Modified: 27 June 2018
"""

import csv
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt

galID = []
ma90 = []
ma50 = []
galFit50 = []
sextractor90 = []
mag = []
sfr = []
redshift = []

#Read in data from GalaxyMeasurements2 file
csvFile = open('GalaxyMeasurements2.csv')
with open('GalaxyMeasurements2.csv') as csvFile:
    readFile = csv.reader(csvFile, delimiter=',')
    
    #Loop though file and append data from each row to list
    count = -1   
    for row in readFile:
        print (row)
        count += 1
        if (count != 0):  
            galID.append(row[0])
            sextractor90.append(float(row[1]))
            ma90.append(float(row[2]))
            galFit50.append(float(row[3]))
            ma50.append(float(row[4]))
            mag.append(float(row[5]))
            sfr.append(float(row[6]))
            redshift.append(float(row[7]))
csvFile.close()   
###################################################################################
#Begin Function to calculate star formation rates
"""
New for this week, read in Nathan's galfit50, find sfr density by 
sfrDens = mySFR/Nathan's galfit50.

Will need to bring in my sfrAbsorber and sfrNoAbsorber
"""
def starFormation(sfRate, galfit50):
    sfrSurfaceDens = sfRate / galfit50
    #print ("SFR Surface Density:", sfrSurfaceDens)
    return (sfrSurfaceDens)
#End starFormation Function
#################################################################################
"""
Call starFormation() function
"""
sfrDens = []
for i in range(0, len(galID)):
    sfrDensity = starFormation(sfr[i], galFit50[i])
    sfrDens.append(sfrDensity)
print ("SFR Surface Density:", sfrDens, '\n')
print ("Length galID:", len(galID))
print ("Length SFR Surface Density:", len(sfrDens))
#End
##################################################################################
"""
Sort sfr densities by absorbers and non-absorbers, write to lists for use
in histogram plots.
"""
#Galaxy ID's listed in abs_data_Target8.csv file
absorberGalID = ['331', '341', '333', '336']
sfrDensAbsorb = []
sfrDensNoAbsorb = []
for i in range(0, len(galID)):
    if (galID[i] in absorberGalID):
        sfrDensAbsorb.append(sfrDens[i])
    else:
        sfrDensNoAbsorb.append(sfrDens[i])
print (len(sfrDensAbsorb))
print (len(sfrDensNoAbsorb))

#Plot sfr density
sfrDensBinArray = np.linspace(0, 4, 12)

plt.hist(sfrDensAbsorb, density=True, histtype='step', label= 'MgII Detection (%i)' %len(sfrDensAbsorb))
plt.hist(sfrDensNoAbsorb, bins=sfrDensBinArray, density=True, histtype='step', label = 'MgII Non-Detection (%i)' %len(sfrDensNoAbsorb))
plt.xlabel("SFR Surface Density $(M_{Sun} yr^{-1} Kpc^{-2})$")
plt.ylabel("Normalized Number")
plt.legend()
plt.savefig('Combined_Histogram_SfrDensity')
plt.show()
#End
##################################################################################
for i in range(0, len(galID)):
    luminosity = mag[i]
    zGal = redshift[i]

"""
#Create scatter plots for Nathan's 90 and 50% and Matthew's 90 and 50%.
x = np.linspace(0, 15, 15)
y = np.linspace(0, 15, 15)
plt.plot(x, y)
plt.scatter(sextractor90, ma90, c='r', label='90% Enclosed')
plt.scatter(galFit50, ma50, c='b', label='50% Enclosed')
print ()

fig, ax = plt.subplots()
im = ax.scatter(galFit50, ma50, c=sfr)
ax.set_xlim([0, 12])
ax.set_ylim([0, 12])
fig.colorbar(im, ax=ax)
im.set_clim(-4, 4)
print ("SFR List:", sfr)

#plt.subplots_adjust(top=1.5)
#plt.subplots_adjust(right=1.5)
plt.axis([0,35, 0,35])
plt.legend(loc='upper left')
plt.xlabel('Sextractor Measurements')
plt.ylabel('Matthews Measurements')
plt.savefig('Emission_Comparison_Plot')
plt.show()
"""
#Statistics
print ("Nathan's sextractor 90% / my 90% statistics")
print (stats.ks_2samp(sextractor90, ma90),'\n')

print ("Nathan's galfit 50% / my 50% statistics")
print (stats.ks_2samp(galFit50, ma50),'\n')


x = np.linspace(0, 15, 15)
y = np.linspace(0, 15, 15)
plt.plot(x, y)
plt.scatter(sextractor90, ma90, c=sfr)
plt.axis([0, 30, 0, 15])
plt.xlabel('sextractor90')
plt.ylabel('ma90')
plt.colorbar()
plt.savefig('sextractor90_vs_ma90_c=sfr')
plt.show()

x = np.linspace(0, 15, 15)
y = np.linspace(0, 15, 15)
plt.plot(x, y)
plt.scatter(sextractor90, ma90, c=sfrDens)
plt.axis([0, 30, 0, 15])
plt.xlabel('sextractor90')
plt.ylabel('ma90')
plt.colorbar()
plt.savefig('sextractor90_vs_ma90_c=sfrDens')
plt.show()

plt.scatter(galFit50, ma50, c=sfr)
plt.axis([0, 12, 0, 12])
plt.xlabel('galFit50')
plt.ylabel('ma50')
plt.colorbar()
plt.savefig('galFit50_vs_ma50_c=sfr')
plt.show()

plt.scatter(galFit50, ma50, c=sfrDens)
plt.axis([0, 12, 0, 12])
plt.xlabel('galFit50')
plt.ylabel('ma50')
plt.colorbar()
plt.savefig('galFit50_vs_ma50_c=sfrDens')
plt.show()


#Go through Nathan and Matthew's 90% lists and find difference between each value.
#Append difference to list named deltaSize90.
deltaSize90 = []
for i in range(0, len(sextractor90)):
    size = sextractor90[i] - ma90[i]
    deltaSize90.append(size)
print (deltaSize90)
print (len(deltaSize90),'\n')

#Go through Nathan and Matthew's 50% lists and find difference between each value.
#Append difference to list named deltaSize50.
deltaSize50 = []
for i in range(0, len(galFit50)):
    size = galFit50[i] - ma50[i]
    deltaSize50.append(size)
print (deltaSize50)
print (len(deltaSize50),'\n')

#Create scatter plot comparing deltaSize with mag column from file.
plt.scatter(mag, deltaSize90, c='r', label='90% Enclosed')
plt.scatter(mag, deltaSize50, c='b', label='50% Enclosed')
plt.legend(loc='upper right')
#plt.axis([0,25,0,25])
plt.xlabel('Magnitude')
plt.ylabel('Difference in Size')
plt.savefig('Emission_mag_Difference_Plot')
plt.show()

#Create scatter plot comparing deltaSize with redshift column from file.
plt.scatter(redshift, deltaSize90, c='r', label='90% Enclosed')
plt.scatter(redshift, deltaSize50, c='b', label='50% Enclosed')
plt.legend(loc='upper right')
plt.xlabel('Redshift')
plt.ylabel('Difference in Size')
plt.savefig('Emission_redshift_Difference_Plot')
plt.show()