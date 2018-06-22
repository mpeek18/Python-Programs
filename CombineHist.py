# -*- coding: utf-8 -*-
"""
Created on Sat Sep 30 20:26:38 2017

@author: Matthew Peek
Last Modified: 22 June 2018
"""

import numpy as np
from astropy.io import ascii
from matplotlib import pyplot as plt
from scipy import stats

#Get field 1 absorber and non-absorber data.
field1DensAbsorb = []
field1DensNoAbsorb = []
field1Absorb = []
field1NoAbsorb = []
field1R90Absorb = []
field1R90NoAbsorb = []
field1R50Absorb = []
field1R50NoAbsorb = []
field1zAbsorb = []
field1zNoAbsorb = []
field1AreaAbsorb = []
field1AreaNoAbsorb = []
field1angDist = []
field1angDistAbsorber = []
field1angDistNoAbsorber = []

field1FileDensAbsorb = ascii.read('Field1_SFRDensAbsorb.dat')
field1FileDensNoAbsorb = ascii.read('Field1_NoSFRDensAbsorb.dat')
field1FileSFRNoAbsorb = ascii.read('Field1_SFRNoAbsorb.dat')
field1FileSFRAbsorb = ascii.read('Field1_SFRAbsorb.dat')
field1FileR90Absorb = ascii.read('Field1_r90Absorb.dat')
field1FileR90NoAbsorb = ascii.read('Field1_r90NoAbsorb.dat')
field1FilezAbsorb = ascii.read('Field1_zDistAbsorb.dat')
field1FilezNoAbsorb = ascii.read('Field1_zDistNoAbsorb.dat')
field1GalAreaAbsorb = ascii.read('Field1_galAreaAbsorb.dat')
field1GalAreaNoAbsorb = ascii.read('Field1_galAreaNoAbsorb.dat')
field1FileR50Absorb = ascii.read('Field1_r50Absorb.dat')
field1FileR50NoAbsorb = ascii.read('Field1_r50NoAbsorb.dat')
field1angularDist = ascii.read('Field1_AngularDistance.dat')
field1angDistAbsorb = ascii.read('Field1_AbsorbAngularDistance.dat')
field1angDistNoAbsorb = ascii.read('Field1_NoAbsorbAngularDistance.dat')

sfrDensAbsorb = field1FileDensAbsorb['col2']
sfrDensNoAbsorb = field1FileDensNoAbsorb['col2']
sfrAbsorb = field1FileSFRAbsorb['col2']
sfrNoAbsorb = field1FileSFRNoAbsorb['col2']
R90Absorb = field1FileR90Absorb['col2']
R90NoAbsorb = field1FileR90NoAbsorb['col2']
zDistAbsorb = field1FilezAbsorb['col2']
zDistNoAbsorb = field1FilezNoAbsorb['col2']
galArea1Absorb = field1GalAreaAbsorb['col2']
galArea1NoAbsorb = field1GalAreaNoAbsorb['col2']
R50Absorb = field1FileR50Absorb['col2']
R50NoAbsorb = field1FileR50NoAbsorb['col2']
angDist = field1angularDist['col3']
angDistAbsorb = field1angDistAbsorb['col2']
angDistNoAbsorb = field1angDistNoAbsorb['col2']

#Append data to new lists.
for i in range(1, len(sfrDensAbsorb)):
    field1DensAbsorb.append(float(sfrDensAbsorb[i]))
for i in range(1, len(sfrDensNoAbsorb)):
    field1DensNoAbsorb.append(float(sfrDensNoAbsorb[i]))
for i in range(1, len(sfrAbsorb)):
    field1Absorb.append(float(sfrAbsorb[i]))
for i in range(1, len(sfrNoAbsorb)):
    field1NoAbsorb.append(float(sfrNoAbsorb[i]))
for i in range(1, len(R90Absorb)):
    field1R90Absorb.append(float(R90Absorb[i]))
for i in range(1, len(R90NoAbsorb)):
    field1R90NoAbsorb.append(float(R90NoAbsorb[i]))
for i in range(1, len(zDistAbsorb)):
    field1zAbsorb.append(float(zDistAbsorb[i]))
for i in range(1, len(zDistNoAbsorb)):
    field1zNoAbsorb.append(float(zDistNoAbsorb[i]))
for i in range(1, len(galArea1Absorb)):
    field1AreaAbsorb.append(float(galArea1Absorb[i]))
for i in range(1, len(galArea1NoAbsorb)):
    field1AreaNoAbsorb.append(float(galArea1NoAbsorb[i]))
for i in range(1, len(R50Absorb)):
    field1R50Absorb.append(float(R50Absorb[i]))
for i in range(1, len(R50NoAbsorb)):
    field1R50NoAbsorb.append(float(R50NoAbsorb[i]))
for i in range(1, len(angDist)):
    field1angDist.append(float(angDist[i]))
for i in range(1, len(angDistAbsorb)):
    field1angDistAbsorber.append(float(angDistAbsorb[i]))
for i in range(1, len(angDistNoAbsorb)):
    field1angDistNoAbsorber.append(float(angDistNoAbsorb[i]))
    
print (field1DensAbsorb,'\n')
print (field1DensNoAbsorb,'\n')
print (field1Absorb,'\n')
print (field1NoAbsorb,'\n')
print (field1R90Absorb,'\n')
print (field1R90NoAbsorb,'\n')
print (field1zAbsorb,'\n')
print (field1zNoAbsorb,'\n')
print (field1AreaAbsorb,'\n')
print (field1AreaNoAbsorb,'\n')
print (field1R50Absorb,'\n')
print (field1R50NoAbsorb,'\n')
print (field1angDist,'\n')
print (field1angDistAbsorber,'\n')
print (field1angDistNoAbsorber,'\n')
print ("End field 1")
print ("-------------------------------------------------------------------------")


#Get field 2 absorber and non-absorber data.
field2DensAbsorb = []
field2DensNoAbsorb = []
field2Absorb = []
field2NoAbsorb = []
field2R90Absorb = []
field2R90NoAbsorb = []
field2zAbsorb = []
field2zNoAbsorb = []
field2AreaAbsorb = []
field2AreaNoAbsorb = []
field2R50Absorb = []
field2R50NoAbsorb = []
field2angDist = []
field2angDistAbsorber = []
field2angDistNoAbsorber = []

field2FileDensAbsorb = ascii.read('Field2_SFRDensAbsorb.dat')
field2FileDensNoAbsorb = ascii.read('Field2_NoSFRDensAbsorb.dat')
field2FileSFRAbsorb = ascii.read('Field2_SFRAbsorb.dat')
field2FileSFRNoAbsorb = ascii.read('Field2_SFRNoAbsorb.dat')
field2FileR90Absorb = ascii.read('Field2_r90Absorb.dat')
field2FileR90NoAbsorb = ascii.read('Field2_r90NoAbsorb.dat')
field2FilezAbsorb = ascii.read('Field2_zDistAbsorb.dat')
field2FilezNoAbsorb = ascii.read('Field2_zDistNoAbsorb.dat')
field2GalAreaAbsorb = ascii.read('Field2_galAreaAbsorb.dat')
field2GalAreaNoAbsorb = ascii.read('Field2_galAreaNoAbsorb.dat')
field2FileR50Absorb = ascii.read('Field2_r50Absorb.dat')
field2FileR50NoAbsorb = ascii.read('Field2_r50NoAbsorb.dat')
field2angularDist = ascii.read('Field2_AngularDistance.dat')
field2angDistAbsorb = ascii.read('Field2_AbsorbAngularDistance.dat')
field2angDistNoAbsorb = ascii.read('Field2_NoAbsorbAngularDistance.dat')


sfrDensAbsorb2 = field2FileDensAbsorb['col2']
sfrDensNoAbsorb2 = field2FileDensNoAbsorb['col2']
sfrAbsorb2 = field2FileSFRAbsorb['col2'] 
sfrNoAbsorb2 = field2FileSFRNoAbsorb['col2']
R90Absorb2 = field2FileR90Absorb['col2']
R90NoAbsorb2 = field2FileR90NoAbsorb['col2']
zDistAbsorb2 = field2FilezAbsorb['col2']
zDistNoAbsorb2 = field2FilezNoAbsorb['col2']
galArea2Absorb = field2GalAreaAbsorb['col2']
galArea2NoAbsorb = field2GalAreaNoAbsorb['col2']
R50Absorb2 = field2FileR50Absorb['col2']
R50NoAbsorb2 = field2FileR50NoAbsorb['col2']
angDist = field2angularDist['col3']
angDistAbsorb = field2angDistAbsorb['col2']
angDistNoAbsorb = field2angDistNoAbsorb['col2']


#Append data to new lists.
for i in range(1, len(sfrDensAbsorb2)):
    field2DensAbsorb.append(float(sfrDensAbsorb2[i]))
for i in range(1, len(sfrDensNoAbsorb2)):
    field2DensNoAbsorb.append(float(sfrDensNoAbsorb2[i]))
for i in range(1, len(sfrAbsorb2)):
    field2Absorb.append(float(sfrAbsorb2[i]))
for i in range(1, len(sfrNoAbsorb2)):
    field2NoAbsorb.append(float(sfrNoAbsorb2[i]))
for i in range(1, len(R90Absorb2)):
    field2R90Absorb.append(float(R90Absorb2[i]))
for i in range(1, len(R90NoAbsorb2)):
    field2R90NoAbsorb.append(float(R90NoAbsorb2[i]))
for i in range(1, len(zDistAbsorb2)):
    field2zAbsorb.append(float(zDistAbsorb2[i]))
for i in range(1, len(zDistNoAbsorb2)):
    field2zNoAbsorb.append(float(zDistNoAbsorb2[i]))
for i in range(1, len(galArea2Absorb)):
    field2AreaAbsorb.append(float(galArea2Absorb[i]))
for i in range(1, len(galArea2NoAbsorb)):
    field2AreaNoAbsorb.append(float(galArea2NoAbsorb[i]))
for i in range(1, len(R50Absorb2)):
    field2R50Absorb.append(float(R50Absorb2[i]))
for i in range(1, len(R50NoAbsorb2)):
    field2R50NoAbsorb.append(float(R50NoAbsorb2[i]))
for i in range(1, len(angDist)):
    field2angDist.append(float(angDist[i]))
for i in range(1, len(angDistAbsorb)):
    field2angDistAbsorber.append(float(angDistAbsorb[i]))
for i in range(1, len(angDistNoAbsorb)):
    field2angDistNoAbsorber.append(float(angDistNoAbsorb[i]))

    
print (field2DensAbsorb,'\n')
print (field2DensNoAbsorb,'\n')
print (field2Absorb,'\n')
print (field2NoAbsorb,'\n')
print (field2R90Absorb,'\n')
print (field2R90NoAbsorb,'\n')
print (field2zAbsorb,'\n')
print (field2zNoAbsorb,'\n')
print (field2AreaAbsorb,'\n')
print (field2AreaNoAbsorb,'\n')
print (field2R50Absorb,'\n')
print (field2R50NoAbsorb,'\n')
print (field2angDist,'\n')
print (field2angDistAbsorber,'\n')
print (field2angDistNoAbsorber,'\n')
print ("End field 2")
print ("-------------------------------------------------------------------------")


#Get field 3 absorber and non-absorber data.
field3DensAbsorb = []
field3DensNoAbsorb = []
field3Absorb = []
field3NoAbsorb = []
field3R90Absorb = []
field3R90NoAbsorb = []
field3zAbsorb = []
field3zNoAbsorb = []
field3AreaAbsorb = []
field3AreaNoAbsorb = []
field3R50Absorb = []
field3R50NoAbsorb = []
field3angDist = []
field3angDistAbsorber = []
field3angDistNoAbsorber = []

field3FileDensAbsorb = ascii.read('Field3_SFRDensAbsorb.dat')
field3FileDensNoAbsorb = ascii.read('Field3_NoSFRDensAbsorb.dat')
field3FileSFRAbsorb = ascii.read('Field3_SFRAbsorb.dat')
field3FileSFRNoAbsorb = ascii.read('Field3_SFRNoAbsorb.dat')
field3FileR90Absorb = ascii.read('Field3_r90Absorb.dat')
field3FileR90NoAbsorb = ascii.read('Field3_r90NoAbsorb.dat')
field3FilezAbsorb = ascii.read('Field3_zDistAbsorb.dat')
field3FilezNoAbsorb = ascii.read('Field3_zDistNoAbsorb.dat')
field3GalAreaAbsorb = ascii.read('Field3_galAreaAbsorb.dat')
field3GalAreaNoAbsorb = ascii.read('Field3_galAreaNoAbsorb.dat')
field3FileR50Absorb = ascii.read('Field3_r50Absorb.dat')
field3FileR50NoAbsorb = ascii.read('Field3_r50NoAbsorb.dat')
field3angularDist = ascii.read('Field3_AngularDistance.dat')
field3angDistAbsorb = ascii.read('Field3_AbsorbAngularDistance.dat')
field3angDistNoAbsorb = ascii.read('Field3_NoAbsorbAngularDistance.dat')


sfrDensAbsorb3 = field3FileDensAbsorb['col2']
sfrDensNoAbsorb3 = field3FileDensNoAbsorb['col2']
sfrAbsorb3 = field3FileSFRAbsorb['col2'] 
sfrNoAbsorb3 = field3FileSFRNoAbsorb['col2']
R90Absorb3 = field3FileR90Absorb['col2']
R90NoAbsorb3 = field3FileR90NoAbsorb['col2']
zDistAbsorb3 = field3FilezAbsorb['col2']
zDistNoAbsorb3 = field3FilezNoAbsorb['col2']
galArea3Absorb = field3GalAreaAbsorb['col2']
galArea3NoAbsorb = field3GalAreaNoAbsorb['col2']
R50Absorb3 = field3FileR50Absorb['col2']
R50NoAbsorb3 = field3FileR50NoAbsorb['col2']
angDist = field3angularDist['col3']
angDistAbsorb = field3angDistAbsorb['col2']
angDistNoAbsorb = field3angDistNoAbsorb['col2']


#Append data to new lists.
for i in range(1, len(sfrDensAbsorb3)):
    field3DensAbsorb.append(float(sfrDensAbsorb3[i]))
for i in range(1, len(sfrDensNoAbsorb3)):
    field3DensNoAbsorb.append(float(sfrDensNoAbsorb3[i]))
for i in range(1, len(sfrAbsorb3)):
    field3Absorb.append(float(sfrAbsorb3[i]))
for i in range(1, len(sfrNoAbsorb3)):
    field3NoAbsorb.append(float(sfrNoAbsorb3[i]))
for i in range(1, len(R90Absorb3)):
    field3R90Absorb.append(float(R90Absorb3[i]))
for i in range(1, len(R90NoAbsorb3)):
    field3R90NoAbsorb.append(float(R90NoAbsorb[i]))
for i in range(1, len(zDistAbsorb3)):
    field3zAbsorb.append(float(zDistAbsorb3[i]))
for i in range(1, len(zDistNoAbsorb3)):
    field3zNoAbsorb.append(float(zDistNoAbsorb3[i]))
for i in range(1, len(galArea3Absorb)):
    field3AreaAbsorb.append(float(galArea3Absorb[i]))
for i in range(1, len(galArea3NoAbsorb)):
    field3AreaNoAbsorb.append(float(galArea3NoAbsorb[i]))
for i in range(1, len(R50Absorb3)):
    field3R50Absorb.append(float(R50Absorb3[i]))
for i in range(1, len(R50NoAbsorb3)):
    field3R50NoAbsorb.append(float(R50NoAbsorb3[i]))
for i in range(1, len(angDist)):
    field3angDist.append(float(angDist[i]))
for i in range(1, len(angDistAbsorb)):
    field3angDistAbsorber.append(float(angDistAbsorb[i]))
for i in range(1, len(angDistNoAbsorb)):
    field3angDistNoAbsorber.append(float(angDistNoAbsorb[i]))

    
print (field3DensAbsorb,'\n')
print (field3DensNoAbsorb,'\n')
print (field3Absorb,'\n')
print (field3NoAbsorb,'\n')
print (field3R90Absorb,'\n')
print (field3R90NoAbsorb,'\n')
print (field3zAbsorb,'\n')
print (field3zNoAbsorb,'\n')
print (field3AreaAbsorb,'\n')
print (field3AreaNoAbsorb,'\n')
print (field3R50Absorb,'\n')
print (field3R50NoAbsorb,'\n')
print (field3angDist,'\n')
print (field3angDistAbsorber,'\n')
print (field3angDistNoAbsorber,'\n')
print ("End field 3")
print ("-------------------------------------------------------------------------")


#Get field 4 absorber and non-absorber data.
field4DensAbsorb = []
field4DensNoAbsorb = []
field4Absorb = []
field4NoAbsorb = []
field4R90Absorb = []
field4R90NoAbsorb = []
field4zAbsorb = []
field4zNoAbsorb = []
field4AreaAbsorb = []
field4AreaNoAbsorb = []
field4R50Absorb = []
field4R50NoAbsorb = []
field4angDist = []
field4angDistAbsorber = []
field4angDistNoAbsorber = []

field4FileDensAbsorb = ascii.read('Field4_SFRDensAbsorb.dat')
field4FileDensNoAbsorb = ascii.read('Field4_NoSFRDensAbsorb.dat')
field4FileSFRNoAbsorb = ascii.read('Field4_SFRNoAbsorb.dat')
field4FileSFRAbsorb = ascii.read('Field4_SFRAbsorb.dat')
field4FileR90Absorb = ascii.read('Field4_r90Absorb.dat')
field4FileR90NoAbsorb = ascii.read('Field4_r90NoAbsorb.dat')
field4FilezAbsorb = ascii.read('Field4_zDistAbsorb.dat')
field4FilezNoAbsorb = ascii.read('Field4_zDistNoAbsorb.dat')
field4GalAreaAbsorb = ascii.read('Field4_galAreaAbsorb.dat')
field4GalAreaNoAbsorb = ascii.read('Field4_galAreaNoAbsorb.dat')
field4FileR50Absorb = ascii.read('Field4_r50Absorb.dat')
field4FileR50NoAbsorb = ascii.read('Field4_r50NoAbsorb.dat')
field4angularDist = ascii.read('Field4_AngularDistance.dat')
field4angDistAbsorb = ascii.read('Field4_AbsorbAngularDistance.dat')
field4angDistNoAbsorb = ascii.read('Field4_NoAbsorbAngularDistance.dat')


sfrDensAbsorb4 = field4FileDensAbsorb['col2']
sfrDensNoAbsorb4 = field4FileDensNoAbsorb['col2']
sfrAbsorb4 = field4FileSFRAbsorb['col2']
sfrNoAbsorb4 = field4FileSFRNoAbsorb['col2']
R90Absorb4 = field4FileR90Absorb['col2']
R90NoAbsorb4 = field4FileR90NoAbsorb['col2']
zDistAbsorb4 = field4FilezAbsorb['col2']
zDistNoAbsorb4 = field4FilezNoAbsorb['col2']
galArea4Absorb = field4GalAreaAbsorb['col2']
galArea4NoAbsorb = field4GalAreaNoAbsorb['col2']
R50Absorb4 = field4FileR50Absorb['col2']
R50NoAbsorb4 = field4FileR50NoAbsorb['col2']
angDist = field4angularDist['col3']
angDistAbsorb = field4angDistAbsorb['col2']
angDistNoAbsorb = field4angDistNoAbsorb['col2']


#Append data to new lists.
for i in range(1, len(sfrDensAbsorb4)):
    field4DensAbsorb.append(float(sfrDensAbsorb4[i]))
for i in range(1, len(sfrDensNoAbsorb4)):
    field4DensNoAbsorb.append(float(sfrDensNoAbsorb4[i]))
for i in range(1, len(sfrAbsorb4)):
    field4Absorb.append(float(sfrAbsorb4[i]))
for i in range(1, len(sfrNoAbsorb4)):
    field4NoAbsorb.append(float(sfrNoAbsorb4[i]))
for i in range(1, len(R90Absorb4)):
    field4R90Absorb.append(float(R90Absorb4[i]))
for i in range(1, len(R90NoAbsorb4)):
    field4R90NoAbsorb.append(float(R90NoAbsorb4[i]))
for i in range(1, len(zDistAbsorb4)):
    field4zAbsorb.append(float(zDistAbsorb4[i]))
for i in range(1, len(zDistNoAbsorb4)):
    field4zNoAbsorb.append(float(zDistNoAbsorb4[i]))
for i in range(1, len(galArea4Absorb)):
    field4AreaAbsorb.append(float(galArea4Absorb[i]))
for i in range(1, len(galArea4NoAbsorb)):
    field4AreaNoAbsorb.append(float(galArea4NoAbsorb[i]))
for i in range(1, len(R50Absorb4)):
    field4R50Absorb.append(float(R50Absorb4[i]))
for i in range(1, len(R50NoAbsorb4)):
    field4R50NoAbsorb.append(float(R50NoAbsorb4[i]))
for i in range(1, len(angDist)):
    field4angDist.append(float(angDist[i]))
for i in range(1, len(angDistAbsorb)):
    field4angDistAbsorber.append(float(angDistAbsorb[i]))
for i in range(1, len(angDistNoAbsorb)):
    field4angDistNoAbsorber.append(float(angDistNoAbsorb[i]))

    
print (field4DensAbsorb,'\n')
print (field4DensNoAbsorb,'\n')
print (field4Absorb,'\n')
print (field4NoAbsorb,'\n')
print (field4R90Absorb,'\n')
print (field4R90NoAbsorb,'\n')
print (field4zAbsorb,'\n')
print (field4zNoAbsorb,'\n')
print (field4AreaAbsorb,'\n')
print (field4AreaNoAbsorb,'\n')
print (field4R50Absorb,'\n')
print (field4R50NoAbsorb,'\n')
print (field4angDist,'\n')
print (field4angDistAbsorber,'\n')
print (field4angDistNoAbsorber,'\n')
print ("End field 4")
print ("-------------------------------------------------------------------------")


#Get field 5 absorber and non-absorber data.
field5DensAbsorb = []
field5DensNoAbsorb = []
field5Absorb = []
field5NoAbsorb = []
field5R90Absorb = []
field5R90NoAbsorb = []
field5zAbsorb = []
field5zNoAbsorb = []
field5AreaAbsorb = []
field5AreaNoAbsorb = []
field5R50Absorb = []
field5R50NoAbsorb = []
field5angDist = []
field5angDistAbsorber = []
field5angDistNoAbsorber = []

field5FileDensAbsorb = ascii.read('Field5_SFRDensAbsorb.dat')
field5FileDensNoAbsorb = ascii.read('Field5_NoSFRDensAbsorb.dat')
field5FileSFRNoAbsorb = ascii.read('Field5_SFRNoAbsorb.dat')
field5FileSFRAbsorb = ascii.read('Field5_SFRAbsorb.dat')
field5FileR90Absorb = ascii.read('Field5_r90Absorb.dat')
field5FileR90NoAbsorb = ascii.read('Field5_r90NoAbsorb.dat')
field5FilezAbsorb = ascii.read('Field5_zDistAbsorb.dat')
field5FilezNoAbsorb = ascii.read('Field5_zDistNoAbsorb.dat')
field5GalAreaAbsorb = ascii.read('Field5_galAreaAbsorb.dat')
field5GalAreaNoAbsorb = ascii.read('Field5_galAreaNoAbsorb.dat')
field5FileR50Absorb = ascii.read('Field5_r50Absorb.dat')
field5FileR50NoAbsorb = ascii.read('Field5_r50NoAbsorb.dat')
field5angularDist = ascii.read('Field5_AngularDistance.dat')
field5angDistAbsorb = ascii.read('Field5_AbsorbAngularDistance.dat')
field5angDistNoAbsorb = ascii.read('Field5_NoAbsorbAngularDistance.dat')


sfrDensAbsorb5 = field5FileDensAbsorb['col2']
sfrDensNoAbsorb5 = field5FileDensNoAbsorb['col2']
sfrAbsorb5 = field5FileSFRAbsorb['col2']
sfrNoAbsorb5 = field5FileSFRNoAbsorb['col2']
R90Absorb5 = field5FileR90Absorb['col2']
R90NoAbsorb5 = field5FileR90NoAbsorb['col2']
zDistAbsorb5 = field5FilezAbsorb['col2']
zDistNoAbsorb5 = field5FilezNoAbsorb['col2']
galArea5Absorb = field5GalAreaAbsorb['col2']
galArea5NoAbsorb = field5GalAreaNoAbsorb['col2']
R50Absorb5 = field5FileR50Absorb['col2']
R50NoAbsorb5 = field5FileR50NoAbsorb['col2']
angDist = field5angularDist['col3']
angDistAbsorb = field5angDistAbsorb['col2']
angDistNoAbsorb = field5angDistNoAbsorb['col2']


#Append data to new lists.
for i in range(1, len(sfrDensAbsorb5)):
    field5DensAbsorb.append(float(sfrDensAbsorb5[i]))
for i in range(1, len(sfrDensNoAbsorb5)):
    field5DensNoAbsorb.append(float(sfrDensNoAbsorb5[i]))
for i in range(1, len(sfrAbsorb5)):
    field5Absorb.append(float(sfrAbsorb5[i]))
for i in range(1, len(sfrNoAbsorb5)):
    field5NoAbsorb.append(float(sfrNoAbsorb5[i]))
for i in range(1, len(R90Absorb5)):
    field5R90Absorb.append(float(R90Absorb5[i]))
for i in range(1, len(R90NoAbsorb5)):
    field5R90NoAbsorb.append(float(R90NoAbsorb5[i]))
for i in range(1, len(zDistAbsorb5)):
    field5zAbsorb.append(float(zDistAbsorb5[i]))
for i in range(1, len(zDistNoAbsorb5)):
    field5zNoAbsorb.append(float(zDistNoAbsorb5[i]))
for i in range(1, len(galArea5Absorb)):
    field5AreaAbsorb.append(float(galArea5Absorb[i]))
for i in range(1, len(galArea5NoAbsorb)):
    field5AreaNoAbsorb.append(float(galArea5NoAbsorb[i]))
for i in range(1, len(R50Absorb5)):
    field5R50Absorb.append(float(R50Absorb5[i]))
for i in range(1, len(R50NoAbsorb5)):
    field5R50NoAbsorb.append(float(R50NoAbsorb5[i]))
for i in range(1, len(angDist)):
    field5angDist.append(float(angDist[i]))
for i in range(1, len(angDistAbsorb)):
    field5angDistAbsorber.append(float(angDistAbsorb[i]))
for i in range(1, len(angDistNoAbsorb)):
    field5angDistNoAbsorber.append(float(angDistNoAbsorb[i]))

    
print (field5DensAbsorb,'\n')
print (field5DensNoAbsorb,'\n')
print (field5Absorb,'\n')
print (field5NoAbsorb,'\n')
print (field5R90Absorb,'\n')
print (field5R90NoAbsorb,'\n')
print (field5zAbsorb,'\n')
print (field5zNoAbsorb,'\n')
print (field5AreaAbsorb,'\n')
print (field5AreaNoAbsorb,'\n')
print (field5R50Absorb,'\n')
print (field5R50NoAbsorb,'\n')
print (field5angDist,'\n')
print (field5angDistAbsorber,'\n')
print (field5angDistNoAbsorber,'\n')
print ("End field 5")
print ("-------------------------------------------------------------------------")


#Get field 7 absorber and non-absorber data.
field7DensAbsorb = []
field7DensNoAbsorb = []
field7Absorb = []
field7NoAbsorb = []
field7R90Absorb = []
field7R90NoAbsorb = []
field7zAbsorb = []
field7zNoAbsorb = []
field7AreaAbsorb = []
field7AreaNoAbsorb = []
field7R50Absorb = []
field7R50NoAbsorb = []
field7angDist = []
field7angDistAbsorber = []
field7angDistNoAbsorber = []

field7FileDensAbsorb = ascii.read('Field7_SFRDensAbsorb.dat')
field7FileDensNoAbsorb = ascii.read('Field7_NoSFRDensAbsorb.dat')
field7FileSFRNoAbsorb = ascii.read('Field7_SFRNoAbsorb.dat')
field7FileSFRAbsorb = ascii.read('Field7_SFRAbsorb.dat')
field7FileR90Absorb = ascii.read('Field7_r90Absorb.dat')
field7FileR90NoAbsorb = ascii.read('Field7_r90NoAbsorb.dat')
field7FilezAbsorb = ascii.read('Field7_zDistAbsorb.dat')
field7FilezNoAbsorb = ascii.read('Field7_zDistNoAbsorb.dat')
field7GalAreaAbsorb = ascii.read('Field7_galAreaAbsorb.dat')
field7GalAreaNoAbsorb = ascii.read('Field7_galAreaNoAbsorb.dat')
field7FileR50Absorb = ascii.read('Field7_r50Absorb.dat')
field7FileR50NoAbsorb = ascii.read('Field7_r50NoAbsorb.dat')
field7angularDist = ascii.read('Field7_AngularDistance.dat')
field7angDistAbsorb = ascii.read('Field7_AbsorbAngularDistance.dat')
field7angDistNoAbsorb = ascii.read('Field7_NoAbsorbAngularDistance.dat')


sfrDensAbsorb7 = field7FileDensAbsorb['col2']
sfrDensNoAbsorb7 = field7FileDensNoAbsorb['col2']
sfrAbsorb7 = field7FileSFRAbsorb['col2']
sfrNoAbsorb7 = field7FileSFRNoAbsorb['col2']
R90Absorb7 = field7FileR90Absorb['col2']
R90NoAbsorb7 = field7FileR90NoAbsorb['col2']
zDistAbsorb7 = field7FilezAbsorb['col2']
zDistNoAbsorb7 = field7FilezNoAbsorb['col2']
galArea7Absorb = field7GalAreaAbsorb['col2']
galArea7NoAbsorb = field7GalAreaNoAbsorb['col2']
R50Absorb7 = field7FileR50Absorb['col2']
R50NoAbsorb7 = field7FileR50NoAbsorb['col2']
angDist = field7angularDist['col3']
angDistAbsorb = field7angDistAbsorb['col2']
angDistNoAbsorb = field7angDistNoAbsorb['col2']


#Append data to new lists.
for i in range(1, len(sfrDensAbsorb7)):
    field7DensAbsorb.append(float(sfrDensAbsorb7[i]))
for i in range(1, len(sfrDensNoAbsorb7)):
    field7DensNoAbsorb.append(float(sfrDensNoAbsorb7[i]))
for i in range(1, len(sfrAbsorb7)):
    field7Absorb.append(float(sfrAbsorb7[i]))
for i in range(1, len(sfrNoAbsorb7)):
    field7NoAbsorb.append(float(sfrNoAbsorb7[i]))
for i in range(1, len(R90Absorb7)):
    field7R90Absorb.append(float(R90Absorb7[i]))
for i in range(1, len(R90NoAbsorb7)):
    field7R90NoAbsorb.append(float(R90NoAbsorb7[i]))
for i in range(1, len(zDistAbsorb7)):
    field7zAbsorb.append(float(zDistAbsorb7[i]))
for i in range(1, len(zDistNoAbsorb7)):
    field7zNoAbsorb.append(float(zDistNoAbsorb7[i]))
for i in range(1, len(galArea7Absorb)):
    field7AreaAbsorb.append(float(galArea7Absorb[i]))
for i in range(1, len(galArea7NoAbsorb)):
    field7AreaNoAbsorb.append(float(galArea7NoAbsorb[i]))
for i in range(1, len(R50Absorb7)):
    field7R50Absorb.append(float(R50Absorb7[i]))
for i in range(1, len(R50NoAbsorb7)):
    field7R50NoAbsorb.append(float(R50NoAbsorb7[i]))
for i in range(1, len(angDist)):
    field7angDist.append(float(angDist[i]))
for i in range(1, len(angDistAbsorb)):
    field7angDistAbsorber.append(float(angDistAbsorb[i]))
for i in range(1, len(angDistNoAbsorb)):
    field7angDistNoAbsorber.append(float(angDistNoAbsorb[i]))

    
print (field7DensAbsorb,'\n')
print (field7DensNoAbsorb,'\n')
print (field7Absorb,'\n')
print (field7NoAbsorb,'\n')
print (field7R90Absorb,'\n')
print (field7R90NoAbsorb,'\n')
print (field7zAbsorb,'\n')
print (field7zNoAbsorb,'\n')
print (field7AreaAbsorb,'\n')
print (field7AreaNoAbsorb,'\n')
print (field7R50Absorb,'\n')
print (field7R50NoAbsorb,'\n')
print (field7angDist,'\n')
print (field7angDistAbsorber,'\n')
print (field7angDistNoAbsorber,'\n')
print ("End field 7")
print ("-------------------------------------------------------------------------")


#Get field 8 absorber and non-absorber data.
field8DensAbsorb = []
field8DensNoAbsorb = []
field8Absorb = []
field8NoAbsorb = []
field8R90Absorb = []
field8R90NoAbsorb = []
field8zAbsorb = []
field8zNoAbsorb = []
field8AreaAbsorb = []
field8AreaNoAbsorb = []
field8R50Absorb = []
field8R50NoAbsorb = []
field8angDist = []
field8angDistAbsorber = []
field8angDistNoAbsorber = []

field8FileDensAbsorb = ascii.read('Field8_SFRDensAbsorb.dat')
field8FileDensNoAbsorb = ascii.read('Field8_NoSFRDensAbsorb.dat')
field8FileSFRNoAbsorb = ascii.read('Field8_SFRNoAbsorb.dat')
field8FileSFRAbsorb = ascii.read('Field8_SFRAbsorb.dat')
field8FileR90Absorb = ascii.read('Field8_r90Absorb.dat')
field8FileR90NoAbsorb = ascii.read('Field8_r90NoAbsorb.dat')
field8FilezAbsorb = ascii.read('Field8_zDistAbsorb.dat')
field8FilezNoAbsorb = ascii.read('Field8_zDistNoAbsorb.dat')
field8GalAreaAbsorb = ascii.read('Field8_galAreaAbsorb.dat')
field8GalAreaNoAbsorb = ascii.read('Field8_galAreaNoAbsorb.dat')
field8FileR50Absorb = ascii.read('Field8_r50Absorb.dat')
field8FileR50NoAbsorb = ascii.read('Field8_r50NoAbsorb.dat')
field8angularDist = ascii.read('Field8_AngularDistance.dat')
field8angDistAbsorb = ascii.read('Field8_AbsorbAngularDistance.dat')
field8angDistNoAbsorb = ascii.read('Field8_NoAbsorbAngularDistance.dat')


sfrDensAbsorb8 = field8FileDensAbsorb['col2']
sfrDensNoAbsorb8 = field8FileDensNoAbsorb['col2']
sfrAbsorb8 = field8FileSFRAbsorb['col2']
sfrNoAbsorb8 = field8FileSFRNoAbsorb['col2']
R90Absorb8 = field8FileR90Absorb['col2']
R90NoAbsorb8 = field8FileR90NoAbsorb['col2']
zDistAbsorb8 = field8FilezAbsorb['col2']
zDistNoAbsorb8 = field8FilezNoAbsorb['col2']
galArea8Absorb = field8GalAreaAbsorb['col2']
galArea8NoAbsorb = field8GalAreaNoAbsorb['col2']
R50Absorb8 = field8FileR50Absorb['col2']
R50NoAbsorb8 = field8FileR50NoAbsorb['col2']
angDist = field8angularDist['col3']
angDistAbsorb = field8angDistAbsorb['col2']
angDistNoAbsorb = field8angDistNoAbsorb['col2']


#Append data to new lists.
for i in range(1, len(sfrDensAbsorb8)):
    field8DensAbsorb.append(float(sfrDensAbsorb8[i]))
for i in range(1, len(sfrDensNoAbsorb8)):
    field8DensNoAbsorb.append(float(sfrDensNoAbsorb8[i]))
for i in range(1, len(sfrAbsorb8)):
    field8Absorb.append(float(sfrAbsorb8[i]))
for i in range(1, len(sfrNoAbsorb8)):
    field8NoAbsorb.append(float(sfrNoAbsorb8[i]))
for i in range(1, len(R90Absorb8)):
    field8R90Absorb.append(float(R90Absorb8[i]))
for i in range(1, len(R90NoAbsorb8)):
    field8R90NoAbsorb.append(float(R90NoAbsorb8[i]))
for i in range(1, len(zDistAbsorb8)):
    field8zAbsorb.append(float(zDistAbsorb8[i]))
for i in range(1, len(zDistNoAbsorb8)):
    field8zNoAbsorb.append(float(zDistNoAbsorb8[i]))
for i in range(1, len(galArea8Absorb)):
    field8AreaAbsorb.append(float(galArea8Absorb[i]))
for i in range(1, len(galArea8NoAbsorb)):
    field8AreaNoAbsorb.append(float(galArea8NoAbsorb[i]))
for i in range(1, len(R50Absorb8)):
    field8R50Absorb.append(float(R50Absorb8[i]))
for i in range(1, len(R50NoAbsorb8)):
    field8R50NoAbsorb.append(float(R50NoAbsorb8[i]))
for i in range(1, len(angDist)):
    field8angDist.append(float(angDist[i]))
for i in range(1, len(angDistAbsorb)):
    field8angDistAbsorber.append(float(angDistAbsorb[i]))
for i in range(1, len(angDistNoAbsorb)):
    field8angDistNoAbsorber.append(float(angDistNoAbsorb[i]))

    
print (field8DensAbsorb,'\n')
print (field8DensNoAbsorb,'\n')
print (field8Absorb,'\n')
print (field8NoAbsorb,'\n')
print (field8R90Absorb,'\n')
print (field8R90NoAbsorb,'\n')
print (field8zAbsorb,'\n')
print (field8zNoAbsorb,'\n')
print (field8AreaAbsorb,'\n')
print (field8AreaNoAbsorb,'\n')
print (field8R50Absorb,'\n')
print (field8R50NoAbsorb,'\n')
print (field8angDist,'\n')
print (field8angDistAbsorber,'\n')
print (field8angDistNoAbsorber,'\n')
print ("End field 8")
print ("-------------------------------------------------------------------------")


#Get field 9 absorber and non-absorber data.
field9DensAbsorb = []
field9DensNoAbsorb = []
field9Absorb = []
field9NoAbsorb = []
field9R90Absorb = []
field9R90NoAbsorb = []
field9zAbsorb = []
field9zNoAbsorb = []
field9AreaAbsorb = []
field9AreaNoAbsorb = []
field9R50Absorb = []
field9R50NoAbsorb = []
field9angDist = []
field9angDistAbsorber = []
field9angDistNoAbsorber = []

field9FileDensAbsorb = ascii.read('Field9_SFRDensAbsorb.dat')
field9FileDensNoAbsorb = ascii.read('Field9_NoSFRDensAbsorb.dat')
field9FileSFRNoAbsorb = ascii.read('Field9_SFRNoAbsorb.dat')
field9FileSFRAbsorb = ascii.read('Field9_SFRAbsorb.dat')
field9FileR90Absorb = ascii.read('Field9_r90Absorb.dat')
field9FileR90NoAbsorb = ascii.read('Field9_r90NoAbsorb.dat')
field9FilezAbsorb = ascii.read('Field9_zDistAbsorb.dat')
field9FilezNoAbsorb = ascii.read('Field9_zDistNoAbsorb.dat')
field9GalAreaAbsorb = ascii.read('Field9_galAreaAbsorb.dat')
field9GalAreaNoAbsorb = ascii.read('Field9_galAreaNoAbsorb.dat')
field9FileR50Absorb = ascii.read('Field9_r50Absorb.dat')
field9FileR50NoAbsorb = ascii.read('Field9_r50NoAbsorb.dat')
field9angularDist = ascii.read('Field9_AngularDistance.dat')
field9angDistAbsorb = ascii.read('Field9_AbsorbAngularDistance.dat')
field9angDistNoAbsorb = ascii.read('Field9_NoAbsorbAngularDistance.dat')


sfrDensAbsorb9 = field9FileDensAbsorb['col2']
sfrDensNoAbsorb9 = field9FileDensNoAbsorb['col2']
sfrAbsorb9 = field9FileSFRAbsorb['col2']
sfrNoAbsorb9 = field9FileSFRNoAbsorb['col2']
R90Absorb9 = field9FileR90Absorb['col2']
R90NoAbsorb9 = field9FileR90NoAbsorb['col2']
zDistAbsorb9 = field9FilezAbsorb['col2']
zDistNoAbsorb9 = field9FilezNoAbsorb['col2']
galArea9Absorb = field9GalAreaAbsorb['col2']
galArea9NoAbsorb = field9GalAreaNoAbsorb['col2']
R50Absorb9 = field9FileR50Absorb['col2']
R50NoAbsorb9 = field9FileR50NoAbsorb['col2']
angDist = field9angularDist['col3']
angDistAbsorb = field9angDistAbsorb['col2']
angDistNoAbsorb = field9angDistNoAbsorb['col2']


#Append data to new lists.
for i in range(1, len(sfrDensAbsorb9)):
    field9DensAbsorb.append(float(sfrDensAbsorb9[i]))
for i in range(1, len(sfrDensNoAbsorb9)):
    field9DensNoAbsorb.append(float(sfrDensNoAbsorb9[i]))
for i in range(1, len(sfrAbsorb9)):
    field9Absorb.append(float(sfrAbsorb9[i]))
for i in range(1, len(sfrNoAbsorb9)):
    field9NoAbsorb.append(float(sfrNoAbsorb9[i]))
for i in range(1, len(R90Absorb9)):
    field9R90Absorb.append(float(R90Absorb9[i]))
for i in range(1, len(R90NoAbsorb9)):
    field9R90NoAbsorb.append(float(R90NoAbsorb9[i]))
for i in range(1, len(zDistAbsorb9)):
    field9zAbsorb.append(float(zDistAbsorb9[i]))
for i in range(1, len(zDistNoAbsorb9)):
    field9zNoAbsorb.append(float(zDistNoAbsorb9[i]))
for i in range(1, len(galArea9Absorb)):
    field9AreaAbsorb.append(float(galArea9Absorb[i]))
for i in range(1, len(galArea9NoAbsorb)):
    field9AreaNoAbsorb.append(float(galArea9NoAbsorb[i]))
for i in range(1, len(R50Absorb9)):
    field9R50Absorb.append(float(R50Absorb9[i]))
for i in range(1, len(R50NoAbsorb9)):
    field9R50NoAbsorb.append(float(R50NoAbsorb9[i]))
for i in range(1, len(angDist)):
    field9angDist.append(float(angDist[i]))
for i in range(1, len(angDistAbsorb)):
    field9angDistAbsorber.append(float(angDistAbsorb[i]))
for i in range(1, len(angDistNoAbsorb)):
    field9angDistNoAbsorber.append(float(angDistNoAbsorb[i]))

    
print (field9DensAbsorb,'\n')
print (field9DensNoAbsorb,'\n')
print (field9Absorb,'\n')
print (field9NoAbsorb,'\n')
print (field9R90Absorb,'\n')
print (field9R90NoAbsorb,'\n')
print (field9zAbsorb,'\n')
print (field9zNoAbsorb,'\n')
print (field9AreaAbsorb,'\n')
print (field9AreaNoAbsorb,'\n')
print (field9R50Absorb,'\n')
print (field9R50NoAbsorb,'\n')
print (field9angDist,'\n')
print (field9angDistAbsorber,'\n')
print (field9angDistNoAbsorber,'\n')
print ("End field 9")
print ("-------------------------------------------------------------------------")


#Combine fields 1, 2, 3, 4, 5, 7, 8, and 9 lists into combined lists, absorbers and non-absorbers.
totalSFRDensAbsorb = (field1DensAbsorb + field2DensAbsorb + field3DensAbsorb + field4DensAbsorb + 
                     field5DensAbsorb + field7DensAbsorb + field8DensAbsorb + field9DensAbsorb)
print (len(totalSFRDensAbsorb))
print ("Combined SFR Density Absorb:", totalSFRDensAbsorb,'\n')

totalSFRDensNoAbsorb = (field1DensNoAbsorb + field2DensNoAbsorb + field3DensNoAbsorb + 
                        field4DensNoAbsorb + field5DensNoAbsorb + field7DensNoAbsorb + 
                        field8DensNoAbsorb + field9DensNoAbsorb)
print (len(totalSFRDensNoAbsorb))
print ("Combined SFR Density Non Absorb:", totalSFRDensNoAbsorb,'\n')

totalSFRDens = (totalSFRDensAbsorb + totalSFRDensNoAbsorb)
print (len(totalSFRDens))
print ("Combined SFR Density:", totalSFRDens,'\n')

totalSFRAbsorb = (field1Absorb + field2Absorb + field3Absorb + field4Absorb + 
                  field5Absorb + field7Absorb + field8Absorb +field9Absorb)
print (len(totalSFRAbsorb))
print ("Combined SFR Absorb:", totalSFRAbsorb,'\n')

totalSFRNoAbsorb = (field1NoAbsorb + field2NoAbsorb + field3NoAbsorb + field4NoAbsorb + 
                    field5NoAbsorb + field7NoAbsorb + field8NoAbsorb + field9NoAbsorb)
print (len(totalSFRNoAbsorb))
print ("Combined SFR Non Absorb:", totalSFRNoAbsorb,'\n')

totalSFR = (totalSFRAbsorb + totalSFRNoAbsorb)
print (len(totalSFR))
print ("Combined SFR:", totalSFR,'\n')

totalR90Absorb = (field1R90Absorb + field2R90Absorb + field3R90Absorb + field4R90Absorb + 
                  field5R90Absorb + field7R90Absorb + field8R90Absorb + field9R90Absorb)
print (len(totalR90Absorb))
print ("Combined R90 Absorbers:", totalR90Absorb,'\n')

totalR90NoAbsorb = (field1R90NoAbsorb + field2R90NoAbsorb +field3R90NoAbsorb + 
                    field4R90NoAbsorb + field5R90NoAbsorb + field7R90NoAbsorb + 
                    field8R90NoAbsorb + field9R90NoAbsorb) 
print (len(totalR90NoAbsorb))
print ("Combined R90 Non-Absorbers:", totalR90NoAbsorb,'\n')

totalR50Absorb = (field1R50Absorb + field2R50Absorb + field3R50Absorb + field4R50Absorb +
                  field5R50Absorb + field7R50Absorb + field8R50Absorb + field9R50Absorb)
print (len(totalR50Absorb))
print ("Combined R50 Absorbers:", totalR50Absorb)

totalR50NoAbsorb = (field1R50NoAbsorb + field2R50NoAbsorb + field3R50NoAbsorb + 
                    field4R50NoAbsorb + field5R50NoAbsorb + field7R50NoAbsorb + 
                    field8R50NoAbsorb + field9R50NoAbsorb)
print (len(totalR50NoAbsorb))
print ("Combined R50 Non-Absorbers:", totalR50NoAbsorb)

totalzAbsorb = (field1zAbsorb + field2zAbsorb + field3zAbsorb + field4zAbsorb + 
                field5zAbsorb + field7zAbsorb + field8zAbsorb + field9zAbsorb)
print (len(totalzAbsorb))
print ("Combined zDist Absorbers:", totalzAbsorb, '\n')

totalzNoAbsorb = (field1zNoAbsorb + field2zNoAbsorb + field3zNoAbsorb + field4zNoAbsorb + 
                  field5zNoAbsorb + field7zNoAbsorb + field8zNoAbsorb + field9zNoAbsorb)
print (len(totalzNoAbsorb),'\n')
print ("Combined zDist Non-Absorbers:", totalzNoAbsorb, '\n')

totalGalAreaAbsorb = (field1AreaAbsorb + field2AreaAbsorb + field3AreaAbsorb + 
                      field4AreaAbsorb + field5AreaAbsorb + field7AreaAbsorb + 
                      field8AreaAbsorb + field9AreaAbsorb)
print (len(totalGalAreaAbsorb))
print ("Combined galaxy area absorbers:", totalGalAreaAbsorb,'\n')

totalGalAreaNoAbsorb = (field1AreaNoAbsorb + field2AreaNoAbsorb + field3AreaNoAbsorb +
                        field4AreaNoAbsorb + field5AreaNoAbsorb + field7AreaNoAbsorb +
                        field8AreaNoAbsorb + field9AreaNoAbsorb)
print (len(totalGalAreaNoAbsorb))
print ("Combined galaxy area Non-Absorbers:", totalGalAreaNoAbsorb,'\n')

totalAngDist = (field1angDist + field2angDist + field3angDist + field4angDist +
                field5angDist + field7angDist + field8angDist + field9angDist)
print (len(totalAngDist))
print ("Combined Angular Distances:", totalAngDist,'\n')

totalAngDistAbsorb = (field1angDistAbsorber + field2angDistAbsorber + field3angDistAbsorber +
                      field4angDistAbsorber + field5angDistAbsorber + field7angDistAbsorber +
                      field8angDistAbsorber + field9angDistAbsorber)
print (len(totalAngDistAbsorb))
print ("Combined Absorber Angular Distances:", totalAngDistAbsorb,'\n')

totalAngDistNoAbsorb = (field1angDistNoAbsorber + field2angDistNoAbsorber + field3angDistNoAbsorber + 
                        field4angDistNoAbsorber + field5angDistNoAbsorber + field7angDistNoAbsorber +
                        field8angDistNoAbsorber + field9angDistNoAbsorber)
print (len(totalAngDistNoAbsorb))
print ("Combined Non-Absorber Angular Distances:", totalAngDistNoAbsorb,'\n')

print ("--------------------------------------------------------------------------")
print ("Statistics")

#Statistics
print ("Total SFR Density Absorbers / Non-Absorbers Statistics")
print (stats.ks_2samp(totalSFRDensAbsorb, totalSFRDensNoAbsorb),'\n')
#print (totalSFRDensAbsorb,'\n')

print ("Total SFR Absorbers / Non-Absorbers Statistics")
print (stats.ks_2samp(totalSFRAbsorb, totalSFRNoAbsorb),'\n')
#print (totalSFRAbsorb,'\n')

print ("Total R90 Absorbers / R90 Non-Absorbers Statistics")
print (stats.ks_2samp(totalR90Absorb, totalR90NoAbsorb),'\n')

print ("Total R50 Absorbers / R50 Non-Absorbers Statistics")
print (stats.ks_2samp(totalR50Absorb, totalR50NoAbsorb),'\n')
print ("--------------------------------------------------------------------------")

#Plot histogram of combined sfr density data.
#binArray local variable for histogram plots. Sets length and bin size.
sfrDensBinArray = np.linspace(0, 3.2, 10)
sfrBinArray = np.linspace(0, 40, 10)

plt.hist(totalSFRDensAbsorb, bins=sfrDensBinArray, density=True, histtype='step', label= 'MgII Detection (%i)' %len(totalSFRDensAbsorb))
plt.hist(totalSFRDensNoAbsorb, bins=sfrDensBinArray, density=True, histtype='step', label = 'MgII Non-Detection (%i)' %len(totalSFRDensNoAbsorb))
plt.xlabel("SFR Surface Density $(M_{Sun} yr^{-1} Kpc^{-2})$")
plt.ylabel("Normalized Number")
plt.legend()
plt.savefig('Histogram_Combined_SfrDensity')
plt.show()

plt.hist(totalSFRDensAbsorb, bins = np.linspace(0, 2, 10), density=True, histtype='step', label= 'MgII Detection (%i)' %len(totalSFRDensAbsorb))
plt.hist(totalSFRDensNoAbsorb, bins = np.linspace(0, 2, 10), density=True, histtype='step', label = 'MgII Non-Detection (%i)' %len(totalSFRDensNoAbsorb))
plt.xlabel("SFR Surface Density $(M_{Sun} yr^{-1} Kpc^{-2})$")
plt.ylabel("Normalized Number")
plt.legend()
plt.savefig('Histogram_Combined_SfrDensity_NoBinArray')
plt.show()

#print (totalSFRDensAbsorb,'\n')
#print (totalSFRDensNoAbsorb,'\n')

fig, ax = plt.subplots()
ax.hist(totalSFRDensAbsorb, bins=sfrDensBinArray, density=True, histtype='step', label= 'MgII Detection (%i)' %len(totalSFRDensAbsorb))
ax.hist(totalSFRDensNoAbsorb, bins=sfrDensBinArray, density=True, histtype='step', label = 'MgII Non-Detection (%i)' %len(totalSFRDensNoAbsorb))
ax.set_ylim([0, 1])
plt.xlabel("SFR Surface Density $(M_{Sun} yr^{-1} Kpc^{-2})$")
plt.ylabel("Normalized Number")
plt.legend()
plt.savefig('Histogram_Combined_SfrDensity_ScaleY')
plt.show()

plt.hist(totalSFRAbsorb, bins=sfrBinArray, density=True, histtype='step', label= 'MgII Detection (%i)' % len(totalSFRAbsorb))
plt.hist(totalSFRNoAbsorb, bins=sfrBinArray, density=True, histtype='step', label = 'MgII Non-Detection (%i)' % len(totalSFRNoAbsorb))
plt.xlabel("SFR $(M_{Sun} yr^{-1})$")
plt.ylabel("Normalized Number")
plt.legend()
plt.savefig('Histogram_Combined_Sfr')
plt.show()

plt.hist(totalSFRAbsorb, density=True, histtype='step', label= 'MgII Detection (%i)' % len(totalSFRAbsorb))
plt.hist(totalSFRNoAbsorb, density=True, histtype='step', label = 'MgII Non-Detection (%i)' % len(totalSFRNoAbsorb))
plt.xlabel("SFR $(M_{Sun} yr^{-1})$")
plt.ylabel("Normalized Number")
plt.legend()
plt.savefig('Histogram_Combined_Sfr_NoBinArray')
plt.show()

plt.hist(totalR90Absorb, bins=np.linspace(0, 14, 10), density=True, histtype='step', label='MgII Detection(%i)' % len(totalR90Absorb))
plt.hist(totalR90NoAbsorb, bins=np.linspace(0, 13, 10), density=True, histtype='step', label='MgII Non-Detection (%i)' % len(totalR90NoAbsorb))
plt.xlabel("Radius 90")
plt.legend(loc='upper left')
plt.savefig('Histogram_Combined_R90')
plt.show()

plt.hist(totalR50Absorb, bins=np.linspace(0, 7, 10), density=True, histtype='step', label='MgII Detection(%i)' % len(totalR50Absorb))
plt.hist(totalR50NoAbsorb, bins=np.linspace(0, 8, 10), density=True, histtype='step', label='MgII Non-Detection(%i)' % len(totalR50NoAbsorb))
plt.xlabel("Radius 50")
plt.legend(loc='upper right')
plt.savefig('Histogram_Combined_R50')
plt.show()
######################################################################################

#totalzAbsorb.sort()
#totalGalAreaAbsorb.sort()
#totalzNoAbsorb.sort()
#totalGalAreaNoAbsorb.sort()

#Begin scatter plot comparing redshifts vs galaxy areas
plt.scatter(totalzAbsorb, totalGalAreaAbsorb, label='MgII Detection(%i)' % len(totalzAbsorb))
plt.scatter(totalzNoAbsorb, totalGalAreaNoAbsorb, label='MgII Non-Detection(%i)' % len(totalzNoAbsorb))
plt.plot(np.unique(totalzAbsorb), np.poly1d(np.polyfit((totalzAbsorb), totalGalAreaAbsorb, 2))(np.unique(totalzAbsorb)), c='b')
plt.plot(np.unique(totalzNoAbsorb), np.poly1d(np.polyfit((totalzNoAbsorb), totalGalAreaNoAbsorb, 2))(np.unique(totalzNoAbsorb)), c='orange')
plt.xlabel("Redshift Distance")
plt.ylabel("Galaxy Area Kpc")
plt.legend(loc='upper left')
plt.savefig('Redshift_vs_GalaxyArea')
plt.show()
########################################################################################

#totalzAbsorb.sort()
#totalR90Absorb.sort()
#totalzNoAbsorb.sort()
#totalR90NoAbsorb.sort()

#Begin scatter plot comparing redshifts vs radius 50
plt.scatter(totalzAbsorb, totalR50Absorb, label='MgII Detection(%i)' % len(totalzAbsorb))
plt.scatter(totalzNoAbsorb, totalR50NoAbsorb, label='MgII Non-Detection(%i)' % len(totalzNoAbsorb))
plt.plot(np.unique(totalzAbsorb), np.poly1d(np.polyfit((totalzAbsorb), totalR50Absorb, 2))(np.unique(totalzAbsorb)), c='b')
plt.plot(np.unique(totalzNoAbsorb), np.poly1d(np.polyfit((totalzNoAbsorb), totalR50NoAbsorb, 2))(np.unique(totalzNoAbsorb)), c='orange')
plt.xlabel("Redshift Distance")
plt.ylabel("50% Emission Enclosed")
plt.legend(loc='upper left')
plt.savefig('Redshift_vs_Radius50')
plt.show()

#Go through totalzAbsorb list and remove galaxies with redshifts less than 0.8
for i in range(0, len(totalzAbsorb)):
    if (totalzAbsorb[i] < 0.8 and totalR90Absorb[i] < 4):
        print ("Redshift < 0.8:", totalzAbsorb[i])
        print ("R90 < 4:", totalR90Absorb[i],'\n')

#Go through totalzNoAbsorb list and remove galaxies with redshifts less than 0.8
for i in range(0, len(totalzNoAbsorb)):
    if (totalzNoAbsorb[i] < 0.8 and totalR90NoAbsorb[i] < 6):
        print ("Redshift < 0.8:", totalzNoAbsorb[i])
        print ("R90 < 6:", totalR90NoAbsorb[i])
        

#totalzAbsorb.remove(0.69801)
#totalR90Absorb.remove(2.3)
#totalzAbsorb.remove(0.73986)
#totalR90Absorb.remove(3.2)
#totalzAbsorb.remove(0.74474)
#totalR90Absorb.remove(3.5)

#totalzNoAbsorb.remove(0.68748)
#totalR90NoAbsorb.remove(3.5)
#totalzNoAbsorb.remove(0.69801)
#totalR90NoAbsorb.remove(6.5)
#totalzNoAbsorb.remove(0.7878)
#totalR90NoAbsorb.remove(6.7)

#Begin scatter plot comparing redshifts vs radius 90
plt.scatter(totalzAbsorb, totalR90Absorb, label='MgII Detection(%i)' % len(totalzAbsorb))
plt.scatter(totalzNoAbsorb, totalR90NoAbsorb, label='MgII Non-Detection(%i)' % len(totalzNoAbsorb))
plt.plot(np.unique(totalzAbsorb), np.poly1d(np.polyfit((totalzAbsorb), totalR90Absorb, 2))(np.unique(totalzAbsorb)), c='b')
plt.plot(np.unique(totalzNoAbsorb), np.poly1d(np.polyfit((totalzNoAbsorb), totalR90NoAbsorb, 2))(np.unique(totalzNoAbsorb)), c='orange')
plt.xlabel("Redshift Distance")
plt.ylabel("90% Emission Enclosed")
plt.legend(loc='lower right')
plt.savefig('Redshift_vs_Radius90')
plt.show()

#Plot total Angular Distance vs. total SFR Densities
plt.scatter(totalAngDist, totalSFRDens)
plt.plot(np.unique(totalAngDist), np.poly1d(np.polyfit((totalAngDist), totalSFRDens, 2))(np.unique(totalAngDist)), c='r')
plt.xlabel("Angular Distance from QSO (kpc)")
plt.ylabel("SFR Densities")
plt.savefig('totalAngDist_vs_totalSFRDens')
plt.show()

#Plot total Angular Distance vs. total SFR
plt.scatter(totalAngDist, totalSFR)
plt.plot(np.unique(totalAngDist), np.poly1d(np.polyfit((totalAngDist), totalSFR, 2))(np.unique(totalAngDist)), c='r')
plt.xlabel("Angular Distance from QSO (kpc)")
plt.ylabel("Star Formation Rates")
plt.savefig('totalAngDist_vs_TotalSFR')
plt.show()

#Plot absorbers and non-absorbers angular distances vs. SFR Density absorbers and non-absorbers
plt.scatter(totalAngDistAbsorb, totalSFRDensAbsorb, label='MgII Detection(%i)' % len(totalAngDistAbsorb), c='r', marker='^')
plt.scatter(totalAngDistNoAbsorb, totalSFRDensNoAbsorb, label='MgII Non-Detection(%i)' % len(totalAngDistNoAbsorb), c='b')
plt.xlabel("Angular Distance from QSO (kpc)")
plt.ylabel("SFR Densities")
plt.legend(loc='upper left')
plt.savefig('totalAngDistAbsorbers_vs_totalSFRDensAbsorbers')
plt.show()

#Plot absorbers and non-absorbers angular distances vs. SFR absorbers and non-absorbers
plt.scatter(totalAngDistAbsorb, totalSFRAbsorb, label='MgII Detection(%i)' % len(totalAngDistAbsorb), c='r', marker='^')
plt.scatter(totalAngDistNoAbsorb, totalSFRNoAbsorb, label='MgII Non-Detection(%i)' % len(totalAngDistAbsorb), c='b')
plt.xlabel("Angular Distance from QSO (kpc)")
plt.ylabel("Star Formation Rates")
plt.legend(loc='upper left')
plt.savefig('totalAngDistAbsorbers_vs_totalSFRAbsorbers')
plt.show()