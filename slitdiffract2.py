#!/bin/env python

# 2014-02 BMcLeod
# Compute the image of the pupil caused by slit diffraction
# In this version we can calculate off-axis


import warnings
import astropy.io.fits as pyfits
import numpy as np
from numpy import *
from matplotlib import pylab
import matplotlib
from numpy.fft import *

from optparse import OptionParser

# Least Recently Used cache
from pylru import lrudecorator

npix         = 1024       # Size of array used to compute FFT
slitheight   = 10         # Long dimension of slit in arcsec -- meant to be comfortably larger than the diffraction limit.

Dpri         = 8.4        # Primary mirror diamter in meters
rseg         = Dpri/2

m2shift      = 0.0469     # Amount M2 pupil shifts per arcminute, referenced to entrance pupil in meters
m2skirtwidth = 0.0005*8.4 # Width of high-emissivity annulus surrounding M2 segments, referenced to entrance pupil
m2cyloff     = 0.126*8.4 # Apparent offset of top m2 cylindrical volume
m2topcyloff  = 0.446*8.4 # Apparent offset of top m2 cylindrical volume

#
# Get the options off the command line
#
parser = OptionParser()
parser.add_option("-w", type="float", dest="wave"       , default=2.4e-6, help="default=%default")
parser.add_option("-s", type="float", dest="slitwidth"  , default=0.5,    help="default=%default")
parser.add_option("-o", type="float", dest="outbkg"     , default=0.05,   help="default=%default")
parser.add_option("-t", type="float", dest="telbkg"     , default=0.1,    help="default=%default")
parser.add_option("-f", type="float", dest="fldang"     , default=0.0,    help="default=%default")
parser.add_option("-m", type="string",dest="m2type"     , default="base", help="default=%default")
parser.add_option("-l", type="string",dest="looptype"   , default="lyot", help="default=%default")
parser.add_option("-r", type="string",dest="root"       , default="em_",  help="default=%default")
parser.add_option("--lyot",type="float", dest="lyotsize", default="4.15", help="default=%default")
parser.add_option("--saveimgs",  action="store_true",  dest="saveimgs"   , default=False, help="default=%default")
parser.add_option("--nonorm",    action="store_true",  dest="nonorm"     , default=False, help="default=%default")

(options,args) = parser.parse_args()

wave      = options.wave
slitwidth = options.slitwidth
outbkg    = options.outbkg
telbkg    = options.telbkg
fldang    = options.fldang
m2type    = options.m2type
looptype  = options.looptype
root      = options.root
lyotsize  = options.lyotsize
saveimgs  = options.saveimgs
nonorm    = options.nonorm

#
# Set up loops
#
if looptype=="lyot":
    fldlist=[fldang]
    # Effect radius of Lyot stop segment diameter in meters
    lyotlist=[3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4.0, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4]
elif looptype=="field":
    # Field radii range from 0 to 5 arcmin
    fldlist=arange(0,5.001,0.5)
    lyotlist=[lyotsize]
else:
    # No loop
    fldlist=[fldang]
    lyotlist=[lyotsize]

deltaphi=30
pupilscale = 24  # pixels / meter

# The total mirror emmisivity is assigned to M1 and M2 in equal proportions
#
m1refl = 1 - telbkg/2
m2refl = 1 - telbkg/2

pupil = zeros((npix,npix))
slit  = zeros((npix,npix))
Ilyot = zeros((npix,npix))


#
# Compute pupil plane coordinates in meters
#
yindx,xindx=indices((npix,npix),dtype=float32)
y = ( yindx - npix/2 ) / pupilscale
x = ( xindx - npix/2 ) / pupilscale

# Makes a circle, approximating the flux on the boundary of the circle using algorithm from IRAF/apphot design doc
# Inside the circle has value 1, outside has 0
# lrudecorator invokes a cache so we don't have to recompute it if it is in the cache.
#
@lrudecorator(100)
def circle (npix, r0, x0=0, y0=0, scale=1):
    y = yindx - npix/2 - y0 * scale
    x = xindx - npix/2 - x0 * scale
    r = sqrt(x*x+y*y)
    c=1-(r-(r0*scale-0.5))
    c[c<0] = 0
    c[c>1] = 1
    return c

# cornerdist: radius of the corner
# ang: angle of hexagon in degrees. 0 means corners on y axis
# dx, dy: offset of the center
# Uses global y,x arrays
def hex(cornerdist, ang=0, dx=0, dy=0):
    xx = (x+dx) * cos(radians(ang)) + (y+dy) * sin(radians(ang))
    yy = (y+dy) * cos(radians(ang)) - (x+dx) * sin(radians(ang))
    hbool = (abs(xx) < cornerdist * sqrt(3) / 2) * (abs(yy) < cornerdist - tan(radians(30)) * abs(xx))
    h = np.zeros((npix,npix))
    h[hbool] = 1.0
    return h
#
# Makes a GMT pupil shape approximated by 7 circles.
# (real pupil has elliptical off-axis segments and various obstructions
#
def makegmtpupil(rseg=Dpri/2, dx=0, dy=0, segdist=8.7):

    segments =[[0.0, rseg, 0.0,     0.0],
               [0.0, rseg, segdist, 0.0],
               [0.0, rseg, segdist, 60.0],
               [0.0, rseg, segdist, 120.0],
               [0.0, rseg, segdist, 180.0],
               [0.0, rseg, segdist, 240.0],
               [0.0, rseg, segdist, 300.0]]

    pupil = zeros((npix,npix))
    for segment in segments:
        [rin, rout, r0, phi] = segment
        if rseg!=0:
            rout=rseg
        x0 = r0 * cos(radians(phi+deltaphi)) + dx
        y0 = r0 * sin(radians(phi+deltaphi)) + dy
        pupil += circle(npix,rout,x0,y0,pupilscale)
    return pupil

#
# Makes a set of ds9 regions in the shape of the GMT pupil
#
def makegmtpupil_ds9(rseg=Dpri/2, dx=0, dy=0, segdist=8.7):
    segments =[[0.0, rseg, 0.0,     0.0],
               [0.0, rseg, segdist, 0.0],
               [0.0, rseg, segdist, 60.0],
               [0.0, rseg, segdist, 120.0],
               [0.0, rseg, segdist, 180.0],
               [0.0, rseg, segdist, 240.0],
               [0.0, rseg, segdist, 300.0]]

    str="physical\n"
    for segment in segments:
        [rin, rout, r0, phi] = segment
        if rseg!=0:
            rout=rseg
        x0 = r0 * cos(radians(phi+deltaphi)) + dx
        y0 = r0 * sin(radians(phi+deltaphi)) + dy
        str=str+"circle(%f,%f,%f)\n" % (x0 * pupilscale + npix/2, y0*pupilscale+npix/2, rout*pupilscale)
    return str


#
# Save a FITS file
#
def savefits(root,name):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        hdu=pyfits.PrimaryHDU(globals()[name])
        hdu.writeto(root+name+'.fits',clobber=True)

######
# OK, lets get to work now
######

#
# Print the header of our output file
#
print "m2type\tfldang\tlyotrad\tthrput\tbkg\tsnratio"
print "------\t------\t-------\t------\t---\t-------"

#
# Define the M1 pupil shape
#
m1pupil = makegmtpupil()

#
# The area outside the M1 pupil has emissivity = 1
# The area within  the M1 pupil has emissivity = 1 - m1refl
#
m1bkg = m1pupil * (1-m1refl) +  (1-m1pupil)

for fldang in fldlist:

    #
    # Make the correct M2 geometry for one of the various cases under consideration
    #
    if m2type=="base":
        # The baseline design has just the M2 segments with a narrow skirt surrounding them

        # The M2 segments looks like the M1 geometry but shifted depending on the field angle
        m2pupil = makegmtpupil(dx=m2shift*fldang)

        # Add the skirt around each segment
        m2skirt = makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth) - m2pupil

        # Within the segments has emissivity = 1 - m2refl
        # The skirt has emissivity 1
        # ( Outside of that has emissivity 0 )
        m2bkg   = m2skirt + m2pupil * (1-m2refl)

    elif m2type=="cyl":     # Cylindrical reference bodies
        m2pupil = makegmtpupil(dx=m2shift*fldang)
        m2skirt = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth) + makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth,segdist=8.7+m2cyloff), 1.0) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)

    elif m2type=="tallcyl": # Tall Cylindrical reference bodies
        m2pupil = makegmtpupil(dx=m2shift*fldang)
        m2skirt = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth) +
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth,segdist=8.7+m2cyloff) +
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth,segdist=8.7+m2topcyloff),
                             1.0) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)

    elif m2type=="bigm2":  # Monolithic M2
        # M2 is one large circle, sized just large enough to go 10 arcmin off axis without vignetting
        m2pupil = circle(npix,(8.4+2*8.7+10*m2shift)/2,0,0,pupilscale)
        m2bkg   = m2pupil * (1-m2refl)

    elif m2type=="bigref": # Monolithic circular reference body
        # Segmented M2, but surrounded by a large high emissivity circle
        m2pupil = makegmtpupil(dx=m2shift*fldang)
        m2skirt = circle(npix,(8.4+2*8.7+10*m2shift)/2,0,0,pupilscale) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)

    elif m2type=="filledcyl": # AdOptica monolithic reference body
        m2pupil = makegmtpupil(dx=m2shift*fldang)
        m2skirt = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth) +
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth,segdist=8.7+m2cyloff) +
                             circle(npix,9.5,m2shift*fldang,0,pupilscale),
                             1.0) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)

    elif m2type=="bigfilledcyl":  # AdOptical monolithic reference body + Oversize the m2 segments by 10mm radius
        m2pupil = makegmtpupil(dx=m2shift*fldang,rseg=rseg+.010*8.4)
        m2skirt = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth+.010*8.4) +
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth+.010*8.4,segdist=8.7+m2cyloff) +
                             circle(npix,9.5,m2shift*fldang,0,pupilscale),
                             1.0) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)

    elif m2type=="bighex3.7": # All segments on their own piezo hexapods attached to big master hexapod modeled as 3.7m disk
        m2pupil = makegmtpupil(dx=m2shift*fldang)
        m2skirt = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth) +
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth,segdist=8.7+m2cyloff) +
                             circle(npix,3.7/2.*8.4,m2shift*fldang,0,pupilscale),
                             1.0) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)
    elif m2type=="bighex3.6": # All segments on their own piezo hexapods attached to big master hexapod modeled as 3.6m disk
        m2pupil = makegmtpupil(dx=m2shift*fldang)
        m2skirt = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth) +
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth,segdist=8.7+m2cyloff) +
                             circle(npix,3.6/2.*8.4,m2shift*fldang,0,pupilscale),
                             1.0) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)
    elif m2type=="bighex3.0": # All segments on their own piezo hexapods attached to big master hexapod modeled as 3.0m disk
        m2pupil = makegmtpupil(dx=m2shift*fldang)
        m2skirt = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth) +
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth,segdist=8.7+m2cyloff) +
                             circle(npix,3.0/2.*8.4,m2shift*fldang,0,pupilscale),
                             1.0) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)
    elif m2type=="round2017": # Round baffle with diameter 3308.2;  From non-segmented Zemax model, actual M2 beam dia = 3167.36
        m2pupil = makegmtpupil(dx=m2shift*fldang)
        m2skirt = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth) +
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth,segdist=8.7+m2cyloff) +
                             circle(npix,3308.2 / 3167.36 * (8.7*2+Dpri) / 2. , m2shift*fldang,0,pupilscale),
                             1.0) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)

    elif m2type=="hex2017":
        # Hexagon with the corners trimmed off
        m2pupil = makegmtpupil(dx=m2shift*fldang)
        m2skirt = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth) +
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth,segdist=8.7+m2cyloff) +
                             hex(11.73*2/sqrt(3),dx=m2shift*fldang) * hex(12.90*2/sqrt(3),dx=m2shift*fldang,ang=90), 1) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)


    # Here we include the direct view of the sky outside of M2
    #
    m2bkg[m2bkg==0]=outbkg

    # The signal in the pupil plane is the product of the M1 and M2 pupil images
    #
    sig = m1pupil * m2pupil

    # Combine the emissivities of M1 seen in M2, M2 seen directly, and sky seen in M1+M2
    #
    bkg  = m2pupil * (1 - (1-m1bkg)*(1-m2bkg)) + (1-m2pupil) * m2bkg  + sig * outbkg * (1-m1bkg) * (1-m2bkg)

    # Now we compute the effects of slit diffraction
    #
    pixsize = wave * 206265. * pupilscale / npix   # arcsec / pix

    # This makes a rectangle in the shape of the slit
    #
    slit[(abs(x)*pupilscale*pixsize < slitwidth/2.)*(abs(y)*pupilscale*pixsize < slitheight/2.0)] = 1.

    # Take the FFT of the slit to get the convolution kernel for the pupil image
    #
    PSFslit = abs(pow(fftshift(fft2(fftshift(slit))),2.))
    PSFslit /= PSFslit.sum()

    # Convolve the signal with the kernel to get the blurred signal in the instrument pupil plane
    #
    sigconv = abs(fftshift(fft2(fft2(fftshift(PSFslit))*fft2(fftshift(sig)))))/pow(npix,2)

    # Convolve the background image with the kernel to get the blurred background in the instrument pupil plane
    #
    bkgconv = abs(fftshift(fft2(fft2(fftshift(PSFslit))*fft2(fftshift(bkg)))))/pow(npix,2)

    # Define the shape of the cold stop
    #
    lyot=makegmtpupil(rseg=lyotsize)

    # Save the shape of the stop as a ds9 region file
    with open(root+"lyot.reg","w") as f: f.write(makegmtpupil_ds9(rseg=lyotsize))

    # Clip the blurred pupil plane images by the cold stop
    # This is the background and the signal that gets through the instrument cold stop
    bkglyot=bkgconv * lyot
    siglyot=sigconv * lyot

    # Save a bunch of images in FITS format for diagnostic purposes
    if saveimgs:
        fitslist=["m2bkg", "m1bkg", "bkg", "sig", "slit", "PSFslit", "sigconv", "bkgconv", "bkglyot", "siglyot"]
        for f in fitslist:
            savefits(root,f)

    # Compute the normalization values for the S/N
    if nonorm:
        # Don't normalize
        sigref = 1
        bkgref = 1
        snrref = 1
    else:
        # Normalize to what you would get without slit diffraction
        sigref=m1pupil.sum()
        bkgref=m1pupil.sum() * (1 - m1refl*m2refl + outbkg*m1refl*m2refl)
        snrref=sigref / sqrt(bkgref)

    for lyotradius in lyotlist:

        lyotmask = makegmtpupil(rseg=lyotradius)

        # Add up the background going through the cold stop
        bkg = (bkgconv * lyotmask).sum()

        # Add up the signal going through the cold stop
        sig = (sigconv * lyotmask).sum()

        # Calculate the S/N assuming we are background limited
        snr = sig / sqrt(bkg)

        # Write the results to the standard output
        print "%s\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.4f" % ( m2type, fldang, lyotradius, sig/sigref, bkg/bkgref, snr/snrref )
