#!/bin/env python

# Compute the image of the pupil caused by slit diffraction
# We compute the PSF 


import pyfits
from numpy import *
from matplotlib import pylab
import matplotlib
from numpy.fft import *

from optparse import OptionParser


r0     =   0.2    # meters
height = 10000    # meters

waver0 =   0.5e-6  # Wavelength for r0
npix   = 1024
slitheight = 10

Dpri = 8.4
gap = 0.3

parser = OptionParser()
parser.add_option("-w", type="float", dest="wave"       , default=2.4e-6)
parser.add_option("-s", type="float", dest="slitwidth"  , default=0.5)
parser.add_option("-o", type="float", dest="outbkg"     , default=1.0)
parser.add_option("-t", type="float", dest="telbkg"     , default=1.0)

(options,args) = parser.parse_args()

wave=options.wave
slitwidth=options.slitwidth
outbkg=options.outbkg
telbkg=options.telbkg

segments =[[0.0, 4.2, 0.0, 0.0],
           [0.0, 4.2, 8.7, 0.0],
           [0.0, 4.2, 8.7, 60.0],
           [0.0, 4.2, 8.7, 120.0],
           [0.0, 4.2, 8.7, 180.0],
           [0.0, 4.2, 8.7, 240.0],
           [0.0, 4.2, 8.7, 300.0]]

deltaphi=30
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

print "Pupilsum ", pupil.sum()
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

PSFslit = abs(pow(fftshift(fft2(fftshift(slit))),2.))
hdu=pyfits.PrimaryHDU(PSFslit)
hdu.writeto('PSFslit.fits',clobber=True)

pupilconv = abs(fftshift(fft2(fft2(fftshift(PSFslit))*fft2(fftshift(pupil)))))/pow(npix,6)
hdu=pyfits.PrimaryHDU(pupilconv)
hdu.writeto('pupilconv.fits',clobber=True)

bkgconv = abs(fftshift(fft2(fft2(fftshift(PSFslit))*fft2(fftshift((1-pupil)*outbkg + pupil*telbkg)))))/pow(npix,6)


hdu=pyfits.PrimaryHDU(bkgconv)
hdu.writeto('bkgconv.fits',clobber=True)

bkgnodiffraction =(abs(fftshift(fft2(fft2(fftshift(PSFslit))*fft2(fftshift(pupil*telbkg)))))/pow(npix,6)).sum()
signodiffraction =(abs(fftshift(fft2(fft2(fftshift(PSFslit))*fft2(fftshift(pupil)))))/pow(npix,6)).sum()

print "Ratios"
print bkgconv.sum() / bkgnodiffraction
print pupilconv.sum() / signodiffraction

snrnodiffraction = signodiffraction / sqrt(bkgnodiffraction)

print "lyotrad", "sig", "bkg", "snr", "signodiffraction", "bkgnodiffraction", "snrnodiffraction", "snr/snrnodiffraction"

for lyotradius in [3.9, 3.95, 4.0, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4]:
    lyotmask = (zeros((npix,npix)) == 1)
    for segment in segments:
        [rin, rout, r0, phi] = segment
        
        x0 = r0 * cos(radians(phi+deltaphi))
        y0 = r0 * sin(radians(phi+deltaphi))
        xseg = x - x0
        yseg = y - y0
        r2 = xseg*xseg + yseg*yseg
        lyotmask = lyotmask + (r2 < lyotradius*lyotradius)

    
    bkg = (bkgconv   * lyotmask).sum()
    sig = (pupilconv * lyotmask).sum()
    snr = sig/sqrt(bkg)


    print lyotradius, sig, bkg, snr, signodiffraction, bkgnodiffraction, snrnodiffraction, snr/snrnodiffraction

