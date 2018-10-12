#!/bin/env python

# 2014-02 BMcLeod
# Compute the image of the pupil caused by slit diffraction
# In this version we can calculate off-axis
# 2015-10 In this version we compute the effect area for various contributors to the signal and noise

import sys
import warnings
import astropy.io.fits as pyfits
import numpy as np
from numpy import *
from matplotlib import pylab
import matplotlib
from numpy.fft import *

from scipy.ndimage.interpolation import zoom

from optparse import OptionParser

# Least Recently Used cache
from pylru import lrudecorator

npix         = 1024       # Size of array used to compute FFT
slitheight   = 10         # Long dimension of slit in arcsec -- meant to be comfortably larger than the diffraction limit.

Dpri         = 8.4        # Primary mirror diamter in meters
rseg         = Dpri/2

m1m2         = 8.4           # Ratio of beam diameter at m1 vs m2 (approximately)
m2shift      = 0.0469        # Amount M2 pupil shifts per arcminute, referenced to entrance pupil in meters

m2cyloff     = 0.126  * m1m2 # Apparent offset of top m2 cylindrical volume
m2topcyloff  = 0.446  * m1m2 # Apparent offset of top m2 cylindrical volume

validloops = ["lyotsegrad", "lyotcenrad", "fldang", "wave", "slitwid", "m2skirt", "m2baff"]

#
# Get the options off the command line
#
parser = OptionParser()
parser.add_option("--wave",    type="float", dest="wave"      , default=2.4e-6, help="units=meters default=%default")
parser.add_option("--slitwid", type="float", dest="slitwid"   , default=0.5,    help="units=arcsec default=%default")
parser.add_option("--fldang",  type="float", dest="fldang"    , default=0.0,    help="units=arcmin default=%default")

parser.add_option("--m2type",  type="string",dest="m2type"    , default="base", help="options: base,cyl,tallcyl,bigm2,filledcyl,hex  default=%default")

parser.add_option("--loopvar", type="string",dest="loopvar"   , default="lyotsegrad", help="options: " + ",".join(validloops)+" default=%default")
parser.add_option("--loopmin", type="float", dest="loopmin"   , default=3.7,    help="default=%default")
parser.add_option("--loopmax", type="float", dest="loopmax"   , default=4.4,    help="default=%default")
parser.add_option("--loopstp", type="float", dest="loopstp"   , default=0.05,   help="default=%default")

parser.add_option("-r", type="string",dest="root"       , default="em_",  help="default=%default")


parser.add_option("--m2skirt",   type="float", dest="m2skirt",    default="0.0005", help="units=meters default=%default")
parser.add_option("--m2baff",    type="float", dest="m2baff",     default="0.0",    help="units=meters default=%default")
parser.add_option("--lyotsegrad",type="float", dest="lyotsegrad", default="4.15",   help="units=meters default=%default")
parser.add_option("--lyotcenrad",type="float", dest="lyotcenrad", default="1.5",    help="units=meters default=%default")
parser.add_option("--saveimgs",  action="store_true",  dest="saveimgs"   , default="False", help="default=%default")

(options,args) = parser.parse_args()

wave      = options.wave
slitwid   = options.slitwid
fldang    = options.fldang
m2type    = options.m2type
loopvar   = options.loopvar
loopmin   = options.loopmin
loopmax   = options.loopmax
loopstp   = options.loopstp
root      = options.root
m2skirt   = options.m2skirt
m2baff    = options.m2baff
lyotsegrad= options.lyotsegrad
lyotcenrad= options.lyotcenrad
saveimgs  = options.saveimgs

#
# Set up loops
#

if loopvar=="none":
    loopmin=0
    loopmax=0
    none=0
elif not loopvar in validloops:
    sys.exit("Error: loopvar must be one of "+" ".join(validloops))


deltaphi=30
pupilscale = 24  # pixels / meter

#
# The total mirror emmisivity is assigned to M1 and M2 in equal proportions
#
pupil = zeros((npix,npix))
slit  = zeros((npix,npix))
Ilyot = zeros((npix,npix))

#
# Compute pupil plane coordinates in meters
#
yindx,xindx=indices((npix,npix),dtype=float32)
y = ( yindx - npix/2 ) / pupilscale
x = ( xindx - npix/2 ) / pupilscale

#
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
    return np.minimum(pupil,1.0)

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

def SlitDiffract(pupil):
    return abs(fftshift(fft2(fft2(fftshift(PSFslit))*fft2(fftshift(pupil)))))/pow(npix,2)

######
# OK, lets get to work now
######

#
# Print the header of our output file
#
tableheader = """
thrutel: Effective area for light from sky hitting M1 and M2
m1emiss: Effective area for thermal emission from the reflective surface of M1
m2emiss: Effective area for thermal emission from the reflective surface of M2
m2struc: Effective area for thermal emission from the structure surrounding M2
m1struc: Effective area for thermal emission from the structure surrounding M1, reflected by M2
drctsky: Effective area for the direct view of the sky seen around the edges of M2

m2type:  Type of M2 geometry
lyotseg: Radius of each segment of the Lyot stop in meters
lyotcen: Radius of central obscuration of the Lyot stop in meters

m2type\tfldang\twave     \tlyotseg\tlyotcen\tm2skirt\tm2baff\tslitwid\tthrutel\tm1emiss\tm2emiss\tm2struc\tm1struc\tdrctsky
------\t------\t----     \t-------\t-------\t-------\t------\t-------\t-------\t-------\t-------\t-------\t-------\t-------""" % ()

print tableheader

#
# Define the M1 pupil shape
#
m1pupil = makegmtpupil()

#
# The area outside the M1 pupil has emissivity = 1
#

m1struc = 1 - m1pupil



for loopval in np.arange(loopmin, loopmax + 1.e-9, loopstp):

    globals()[loopvar] = loopval

    m2pupil = makegmtpupil(dx=m2shift*fldang)

    m2strucwidth = m2skirt * m1m2 # Width of high-emissivity annulus surrounding M2 segments, referenced to entrance pupil
    #
    # Make the correct M2 geometry for one of the various cases under consideration
    #
    if m2type=="base":
        # The baseline design has just the M2 segments with a narrow skirt surrounding them, i.e. m2baff=0

        # The M2 segments looks like the M1 geometry but shifted depending on the field angle

        # Add the skirt around each segment, and an overall circular obscuration
        m2obscure = np.minimum (
            makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2strucwidth) +
            circle(npix, m2baff*m1m2,m2shift*fldang,0,pupilscale),
            1.0 )



        # Within the segments has emissivity = 1 - m2refl
        # The skirt has emissivity 1
        # ( Outside of that has emissivity 0 )

    elif m2type=="cyl":     # Cylindrical reference bodies
        m2obscure = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2strucwidth) + makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2strucwidth,segdist=8.7+m2cyloff), 1.0)

    elif m2type=="tallcyl": # Tall Cylindrical reference bodies
        m2obscure = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2strucwidth) +
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2strucwidth,segdist=8.7+m2cyloff) +
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2strucwidth,segdist=8.7+m2topcyloff),
                             1.0)

    elif m2type=="bigm2":  # Monolithic M2
        # M2 is one large circle, sized just large enough to go 10 arcmin off axis without vignetting
        m2pupil = circle(npix,(8.4+2*8.7+10*m2shift)/2,0,0,pupilscale)
        m2obscure = m2pupil

    elif m2type=="bigref": # Monolithic circular reference body
        # Segmented M2, but surrounded by a large high emissivity circle
        m2obscure = circle(npix,(8.4+2*8.7+10*m2shift)/2,0,0,pupilscale)

    elif m2type=="filledcyl": # AdOptica monolithic reference body
        m2obscure = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2strucwidth) +
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2strucwidth,segdist=8.7+m2cyloff) +
                             circle(npix,9.5,m2shift*fldang,0,pupilscale),
                             1.0)

    elif m2type=="bigfilledcyl":  # AdOptical monolithic reference body + Oversize the m2 segments by 10mm radius
        m2obscure = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2strucwidth+.010*m1m2) +
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2strucwidth+.010*m1m2,segdist=8.7+m2cyloff) +
                             circle(npix,9.5,m2shift*fldang,0,pupilscale),
                             1.0)
    elif m2type=="hex":
        # Hexagon with the corners trimmed off
        m2obscure = hex(11.73*2/sqrt(3),dx=m2shift*fldang) * hex(12.90*2/sqrt(3),dx=m2shift*fldang,ang=90)
        savefits(root,"m2obscure")

    # Here we include the direct view of the sky outside of M2
    #

    m2struc = m2obscure - m2pupil

    drctsky = 1 - m2obscure

    #
    # Scale down M2 obscuration to figure out how much light is blocked on the way to M1
    m2obscure_demag = np.zeros((npix,npix))
    zoomed = zoom(m2obscure, 1./m1m2)
    (newsizex,newsizey) = zoomed.shape
    left = npix/2 - newsizex/2
    m2obscure_demag[left:left+newsizey,left:left+newsizex] = zoomed

    m1emiss = m1pupil * m2pupil

    thrutel = m1pupil * m2pupil * ( 1 - m2obscure_demag )

    m2emiss = m2pupil

    # Now we compute the effects of slit diffraction
    #
    pixsize = wave * 206265. * pupilscale / npix   # arcsec / pix

    # This makes a rectangle in the shape of the slit
    #
    slit[(abs(x)*pupilscale*pixsize < slitwid/2.)*(abs(y)*pupilscale*pixsize < slitheight/2.0)] = 1.

    # Take the FFT of the slit to get the convolution kernel for the pupil image
    #
    PSFslit = abs(pow(fftshift(fft2(fftshift(slit))),2.))
    PSFslit /= PSFslit.sum()

    # Save the shape of the stop as a ds9 region file
    with open(root+"lyot.reg","w") as f: f.write(makegmtpupil_ds9(rseg=lyotsegrad))

    # Define the shape of the cold stop
    #
    lyot = makegmtpupil(rseg=lyotsegrad)
    lyot -= circle(npix,lyotcenrad*pupilscale)

    # Clip the blurred pupil plane images by the cold stop
    # This is the signal that gets through the instrument cold stop
    thrutel_fp = lyot * SlitDiffract(thrutel)

    # And these are the various sources of background
    m1emiss_fp = lyot * SlitDiffract(m1emiss)
    m2emiss_fp = lyot * SlitDiffract(m2pupil)
    m2struc_fp = lyot * SlitDiffract(m2struc)
    m1struc_fp = lyot * SlitDiffract(m1struc * m2pupil)
    drctsky_fp = lyot * SlitDiffract(drctsky)

    # Add up the effective area passing through the stop

    thrutel_ea = thrutel_fp.sum() / pupilscale**2
    m1emiss_ea = m1emiss_fp.sum() / pupilscale**2
    m2emiss_ea = m2emiss_fp.sum() / pupilscale**2
    m2struc_ea = m2struc_fp.sum() / pupilscale**2
    m1struc_ea = m1struc_fp.sum() / pupilscale**2
    drctsky_ea = drctsky_fp.sum() / pupilscale**2

    # Write the results to the standard output
    print "%s\t%6.3f\t%6.2e\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f" % (
      m2type, fldang, wave, lyotsegrad, lyotcenrad, m2skirt, m2baff, slitwid, thrutel_ea, m1emiss_ea, m2emiss_ea, m2struc_ea, m1struc_ea, drctsky_ea )

    # Save a bunch of images in FITS format for diagnostic purposes, first time through the loop only
    if saveimgs:
        fitslist = ["slit", "PSFslit", "thrutel", "m1emiss", "m2emiss", "m2struc", "m1struc", "drctsky","thrutel_fp", "m2emiss_fp", "m2struc_fp", "m1struc_fp", "drctsky_fp"]
        for f in fitslist:
            savefits(root,f)
        saveimgs = 0
