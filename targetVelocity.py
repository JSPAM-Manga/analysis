
from astropy.io import fits
import numpy as np
#from PIL import Image


import matplotlib.pyplot as plt
import numpy as np


def readTargetInfo(tdirectory, galaxyID):
    
    vdone = []
    fn = tdirectory +  "target_info.txt"
    f = open(fn,"r")
    for l in f:
        v = l.strip().split()
        if len(v) > 2:
            if v[1] == galaxyID:
                vdone = v

    return vdone

def readManga(fn, printHeader, printData):

    hdulist = fits.open(fn)
    ii = 4
    if printHeader:
        print (repr(hdulist[ii].header))
    

    # grab the critical header information
    naxis1 = hdulist[ii].header['NAXIS1']
    naxis2 = hdulist[ii].header['NAXIS2']
    
    dx =  hdulist[ii].header['PC1_1']
    dy =  hdulist[ii].header['PC2_2']
    xc =  hdulist[ii].header['CRVAL1']
    yc =  hdulist[ii].header['CRVAL2']
    
    nxc = hdulist[ii].header['CRPIX1']
    nyc = hdulist[ii].header['CRPIX2']
    

    # get the datablocks for the file
    a = hdulist[ii].data


    # get the data for the stars
    vstars = hdulist[17].data

    w = naxis1
    h = naxis2
    vgas = []
    for i in range(w):
        vgasRow = []
        for j in range(h):

            # position of the current pixel
            ra =  (i - nxc) * dx + xc
            dec = (j - nyc) * dy + yc 
            vstar = vstars[i,j]  #stellar velocity - floating point
            
            vblock = []  # string values of gas velocity
            vb2 = []     # floating point values of gas velocity
            for k in range(8):
                v = a[k, i, j]  # gas velocity for the kth plane
                vblock.append(str(v))
                vb2.append(v)
            vgasAverage = np.average(vb2)
            vgasRow.append(vgasAverage)
            ss = str(i) + ", " + str(j) + ", "+ str(ra) + ", " + str(dec) + ", " + str(vstar) + ", " + str(vgasAverage)+ ", " +", ".join(vblock)
            if vstar != 0 and printData  == 1:
                print ss
            elif printData == 2:
                print ss
            elif printData == 3:
                if vstar <> 0:
                    print vstar, vgasAverage
                
        vgas.append(vgasRow)

    return vstars, vgas, w, h


def dumpVelocities(vstars, vgas, w, h):
    for i in range(w):
        for j in range(h):    
            vb = vgas[:,i,j]
            vgasAverage = np.average(vb)
            vblock = []
            for vv in vb:
                vblock.append(str(vv))
            ss = str(i) + ", " + str(j) +  ", " + str(vstar) + ", " + str(vgasAverage)+ ", " +", ".join(vblock)
            print ss
            

######

fn = "astro/mangadap-8481-6103-default.fits"
printHeader=0
printData = 3

vstars, vgas, w, h = readManga(fn, printHeader, printData)

#dumpVelocities(vstars, vgas, w, h)

#print vgas
tdirectory = "./"
galaxyID = "1237678620102623480"

galaxyPositions = readTargetInfo(tdirectory, galaxyID)
print galaxyPositions


vmin = 100000
vmax = -10000

vg = np.array(vgas)
for i in range(w):
    for j in range(h):
        v = vg[i,j]
        if v< vmin:
            vmin = v
        if v > vmax:
            vmax = v

vmin = -800
vmax = 800
for i in range(w):
    for j in range(h):
        v = vg[i,j]
        v =  (v-vmin) / (vmax -vmin)
        if v>1 or v < 0:
            v = 0
        vg[i,j] = v
        
        print vg[i,j]

print "dumped"
plt.imshow(vg, cmap='hot', interpolation='nearest')
plt.show()
