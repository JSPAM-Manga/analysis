import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import math
import os

class sdssImage:
    def __init__(self, target_info_line):
        self.valid = False

        self.iscale = 0.2
        self.iheight = 512
        self.iwidth = 512

        self.nx = None
        self.ny = None
        self.ixcenter = iwidth/2.
        self.iycenter = iheight/2.
        self.raFactor = 1. / math.cos( float(self.dec) * math.pi/180.)
        self.image = None
        self.ipath = 'img/'+self.name+".jpg"

        v = target_info_line.strip().split("\t")
        if len(v) > 2 :
            self.name = v[1]
            self.ra = float(v[2])
            self.dec = float(v[3])
            self.ra2 = float(v[5])
            self.dec2 = float(v[6])
            self.deltaRA = (self.ra2 - self.ra) * self.raFactor
            self.deltaDec = (self.dec2 - self.dec)

            # separation between galaxies in pixels on the SDSS image
            # convert ra and dec offsets into arcseconds, then scale into pixels
            self.dxPixels =  -self.deltaRA * 3600 / self.iscale
            self.dyPixels =  -self.deltaDec * 3600 / self.iscale

            self.valid = True

        if valid :
            self.grabSDSSimage()
            self.readImage()

    def grabSDSSimage(self):
        if not os.path.exists('img'):
            os.makedirs('img')
        if not os.path.isfile(self.ipath):
            cmd = "wget -O "+self.ipath+
                " \"http://skyservice.pha.jhu.edu/DR7/ImgCutout/getjpeg.aspx"+
                "?ra="+str(self.ra)+
                "&dec="+str(self.dec)+
                "&scale="+str(self.iscale)+
                "&width="+str(self.iheight)+
                "&height="+str(self.iwidth)+
                "&opt=L"+
                "&query="+
                "&Label=on\""
            os.system(s)

    def readImage(self):
        self.image = mpimg.imread(self.ipath)

    def radec2pixel(self,ra,dec):
        x = ra * self.dxPixels
        y = dec * self.dyPixels
        return x,y

    def pixel2radec(self):
        x = ra / self.dxPixels
        y = dec / self.dyPixels
        return x,y

    def plotField(self,alpha=1.):
        plt.imshow(self.image, alpha=alpha)
