#!/bin/env python

# Compute the pupil illumination through a slit, taking the time
# average of many realizations of atmospheric turbulence.


import pyfits
from numpy import *
from matplotlib import pylab
import matplotlib
from numpy.fft import *

r0     =   0.2    # meters
height = 10000    # meters

waver0 =   0.5e-6  # Wavelength for r0
npix   = 1024
wave   =   2.5e-6  # Wavelength in microns
slitwidth = 0.2
slitheight = 10

Dpri = 8.4
gap = 0.3
dt  = 0.1 # sec
windx       =  25       # m/sec
windy       =   0       # m/sec

h = array([ 25 , 275, 425, 1250, 4000, 8000, 13000   ])
f = array([ .126, .087, .067, .350, .227, .068, .075 ])
angles = [1,2,4,8,15,30,60,120]

segments =[[0.0, 4.2, 0.0, 0.0],
           [0.0, 4.2, 8.7, 0.0],
           [0.0, 4.2, 8.7, 60.0],
           [0.0, 4.2, 8.7, 120.0],
           [0.0, 4.2, 8.7, 180.0],
           [0.0, 4.2, 8.7, 240.0],
           [0.0, 4.2, 8.7, 300.0]]

deltaphi=30

screendir='Screens'
screenfiles=['screen4096.1.fits','screen4096.2.fits']


pupilscale = 12  # pixels / meter

# Scale r0 to observation wavelength and to pixels
r00 = r0 * pow(wave/waver0, 6./5.)
r0pix = r00 * pupilscale
print "pupilscale", pupilscale, "pixels/meter"
print "r0 at ",wave*1e6 , "microns:", r00, "meter; ", r0pix, "pixels"


pupil=zeros((npix,npix))
slit=zeros((npix,npix))
Ilyot=zeros((npix,npix))

y,x=indices(pupil.shape,dtype=float32)
y = (y - npix/2) / pupilscale
x = (x - npix/2) / pupilscale

# Make pupil mask

# Initialize array of boolean zeros (probably a better way)
mask = (zeros((npix,npix)) == 1)

for segment in segments:
    [rin, rout, r0, phi] = segment
    x0 = r0 * cos(radians(phi+deltaphi))
    y0 = r0 * sin(radians(phi+deltaphi))
    print x0, y0
    xseg = x - x0
    yseg = y - y0
    r2 = xseg*xseg + yseg*yseg
    mask = mask + (r2 < rout*rout) * (r2 > rin*rin)
    
pixsize = wave * 206265. * pupilscale / npix   # arcsec / pix
pupil[mask] = 1.
slit[(abs(x)*pupilscale*pixsize < slitwidth/2.)*(abs(y)*pupilscale*pixsize < slitheight/2.0)] = 1.
     
print "Pixsize ", pixsize
print "npix ", npix

#
# Compute PSF with no atmosphere
#
hdu=pyfits.PrimaryHDU(pupil)
hdu.writeto('pupil.fits',clobber=True)

# Slit function
hdu=pyfits.PrimaryHDU(slit)
hdu.writeto('slit.fits',clobber=True)

y,x=indices((npix,npix))

for screenfile in screenfiles:
    # Read in an atmospheric phase screen
    atmo = pyfits.getdata(screendir+'/'+screenfile)
    
    # Scale amplitude to r0
    ir0pix = int(r0pix)
    atmo[0]
    print atmo[0]
    print atmo[ir0pix]
    d = atmo[0:-ir0pix] - atmo[ir0pix:]
    print (d*d).mean()
    scale = sqrt(6.88 / (d*d).mean())
    atmo *=scale
    
    dytot = 0

    pupnx=len(atmo[0])
    pupny=len(atmo)

    dx = windx * dt * pupilscale
    dy = npix
    print (pupnx - npix)/dx

    for j in range((pupny-npix)/dy):
        for tick in range((pupnx-npix)/dx):
            time = tick * dt


            # Extract the right part of the phase screen

            dxtot = int(windx * time * pupilscale)
            dytot = int(windy * time * pupilscale) + j * npix
            phase=atmo[dytot:dytot+npix,dxtot:dxtot+npix].copy()

            # Telescope pupil
            Etel = pupil * exp( complex(0,1) * phase)
            
            # Slit plane
            Eslit = fftshift(fft2(fftshift(Etel))) * slit
            Islit = abs(Eslit*Eslit)
            
            # Just before lyot stop
            Elyot = fftshift(fft2(fftshift(Eslit)))
            Ilyot += abs(Elyot*Elyot)
            print screenfile,j,tick


hdu=pyfits.PrimaryHDU(Ilyot)
hdu.writeto('pupilintens.fits',clobber=True)




