#!/bin/env python

# 2014-02 BMcLeod
# Compute the image of the pupil caused by slit diffraction
# In this version we can calculate off-axis

import os
import warnings
import pyfits
import numpy as np
from numpy import *
from matplotlib import pylab
import matplotlib
from numpy.fft import *

from optparse import OptionParser

# Least Recently Used cache
from pylru import lrudecorator

npix   = 1024
slitheight = 10

Dpri = 8.4
rseg=Dpri/2

m2shift = 0.0469     # Amount M2 pupil shifts per arcminute, referenced to entrance pupil in meters
m2skirtwidth = 0.0005*8.4 # Width of high-emissivity annulus surrounding M2 segments, referenced to entrance pupil
m2cyloff=0.126*8.4 # Apparent offset of top m2 cylindrical volume
m2topcyloff=0.446*8.4 # Apparent offset of top m2 cylindrical volume

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
parser.add_option("--saveimgs",  action="store_true",  dest="saveimgs"   , default="False", help="default=%default")

(options,args) = parser.parse_args()

wave=options.wave
slitwidth=options.slitwidth
outbkg=options.outbkg
telbkg=options.telbkg
fldang=options.fldang
m2type=options.m2type
looptype=options.looptype
root=options.root
lyotsize=options.lyotsize
saveimgs=options.saveimgs


if looptype=="lyot":
    fldlist=[fldang]
    lyotlist=[3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 4.0, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4]
elif looptype=="field":
    fldlist=arange(0,5.001,0.5)
    lyotlist=[lyotsize]
else:
    fldlist=[fldang]
    lyotlist=[lyotsize]
    
deltaphi=30
pupilscale = 24  # pixels / meter
m1refl = 1 - telbkg/2
m2refl = 1 - telbkg/2

pupil=zeros((npix,npix))
slit=zeros((npix,npix))
Ilyot=zeros((npix,npix))


yindx,xindx=indices((npix,npix),dtype=float32)
y = ( yindx - npix/2 ) / pupilscale
x = ( xindx - npix/2 ) / pupilscale

# Makes a circle, approximating the flux on the boundary of the circle using algorithm from IRAF/apphot design doc
@lrudecorator(100)
def circle (npix, r0, x0=0, y0=0, scale=1):
    y = yindx - npix/2 - y0 * scale
    x = xindx - npix/2 - x0 * scale
    r = sqrt(x*x+y*y)
    c=1-(r-(r0*scale-0.5))
    c[c<0] = 0
    c[c>1] = 1
    return c

def insiderects ( rectlist, scale=1 ):

    tot = zeros(y.shape)
    for rect in rectlist:
        print rect
        (x0, y0, xwid, ywid, rot, xpixsize, ypixsize) = rect
        x1 = x- x0
        y1 = y- y0
        xx =  x1 * cos(radians(rot)) + y1 * sin(radians(rot))
        yy = -x1 * sin(radians(rot)) + y1 * cos(radians(rot))

        xx=-(abs(xx) - xwid/2) * scale + 0.5
        xx[xx<0] = 0
        xx[xx>1] = 1

        yy=-(abs(yy) - ywid/2) * scale + 0.5
        yy[yy<0] = 0
        yy[yy>1] = 1

        tot += xx * yy

    return tot


# Make GMT pupil shape
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

def makephasingapertures(gap, segdist=8.7, phwid=1.5):

    rectlist=[]
    rot = arange(12.) * 60.0 + 30
    rot[6:] += 30
    x0 = cos(radians(rot)) * segdist/2.
    y0 = sin(radians(rot)) * segdist/2.
    x0[6:] *= 7.643 / (segdist/2)
    y0[6:] *= 7.643 / (segdist/2)

    for (x,y,r) in zip(x0,y0,rot):
        rect=[x,y,phwid,phwid,r,pupilscale,pupilscale]
        print rect
        rectlist.append(rect)
        
    return insiderects(rectlist, pupilscale) * makegmtpupil(rseg=(segdist-gap)/2)

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

def savefits(root,name):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        hdu=pyfits.PrimaryHDU(globals()[name])
        hdu.writeto(root+name+'.fits',clobber=True)
    
m1pupil = makegmtpupil()

phz=makephasingapertures(0.9)
savefits("","phz")
exit()

# Area outside M1 pupil has emissivity=1
m1bkg = m1pupil * (1-m1refl) +  (1-m1pupil)

print "m2type\tfldang\tlyotrad\tthrput\tbkg\tsnratio"
print "------\t------\t-------\t------\t---\t-------"

for fldang in fldlist:

    
    if m2type=="base":
        m2pupil = makegmtpupil(dx=m2shift*fldang)
        m2skirt = makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)
    elif m2type=="cyl":
        m2pupil = makegmtpupil(dx=m2shift*fldang)
        m2skirt = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth) + makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth,segdist=8.7+m2cyloff), 1.0) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)
    elif m2type=="tallcyl":
        m2pupil = makegmtpupil(dx=m2shift*fldang)
        m2skirt = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth) + 
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth,segdist=8.7+m2cyloff) + 
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth,segdist=8.7+m2topcyloff),
                             1.0) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)
    elif m2type=="bigm2":
        m2pupil = circle(npix,(8.4+2*8.7+10*m2shift)/2,0,0,pupilscale)
        m2bkg   = m2pupil * (1-m2refl)
    elif m2type=="bigref":
        m2pupil = makegmtpupil(dx=m2shift*fldang)
        m2skirt = circle(npix,(8.4+2*8.7+10*m2shift)/2,0,0,pupilscale) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)
    elif m2type=="filledcyl":
        m2pupil = makegmtpupil(dx=m2shift*fldang)
        m2skirt = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth) + 
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth,segdist=8.7+m2cyloff) +
                             circle(npix,9.5,m2shift*fldang,0,pupilscale), 
                             1.0) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)
    elif m2type=="bigfilledcyl":
        # Oversize the m2 segments by 10mm radius
        m2pupil = makegmtpupil(dx=m2shift*fldang,rseg=rseg+.010*8.4)
        m2skirt = np.minimum(makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth+.010*8.4) + 
                             makegmtpupil(dx=m2shift*fldang,rseg=rseg+m2skirtwidth+.010*8.4,segdist=8.7+m2cyloff) +
                             circle(npix,9.5,m2shift*fldang,0,pupilscale), 
                             1.0) - m2pupil
        m2bkg   = m2skirt + m2pupil * (1-m2refl)


    m2bkg[m2bkg==0]=outbkg
    sig = m1pupil * m2pupil

    
    # Combine the emissivities of M1 seen in M2, M2 seen directly, and sky seen in M1+M2
    bkg  = m2pupil * (1 - (1-m1bkg)*(1-m2bkg)) + (1-m2pupil) * m2bkg  + sig * outbkg * (1-m1bkg) * (1-m2bkg)

    # Slit function
    #
    pixsize = wave * 206265. * pupilscale / npix   # arcsec / pix
    slit[(abs(x)*pupilscale*pixsize < slitwidth/2.)*(abs(y)*pupilscale*pixsize < slitheight/2.0)] = 1.
    PSFslit = abs(pow(fftshift(fft2(fftshift(slit))),2.))
    PSFslit /= PSFslit.sum()

    sigconv = abs(fftshift(fft2(fft2(fftshift(PSFslit))*fft2(fftshift(sig)))))/pow(npix,2)
    bkgconv = abs(fftshift(fft2(fft2(fftshift(PSFslit))*fft2(fftshift(bkg)))))/pow(npix,2)

    lyot=makegmtpupil(rseg=lyotsize)

    with open(root+"lyot.reg","w") as f: f.write(makegmtpupil_ds9(rseg=lyotsize))
    
    bkglyot=bkgconv * lyot
    siglyot=sigconv * lyot

    if saveimgs:
        fitslist=["m2bkg", "m1bkg", "bkg", "sig", "slit", "PSFslit", "sigconv", "bkgconv", "bkglyot", "siglyot"]
        for f in fitslist:
            savefits(root,f)
    
    sigref=m1pupil.sum()
    bkgref=m1pupil.sum() * (1 - m1refl*m2refl + outbkg*m1refl*m2refl) 
    snrref=sigref / sqrt(bkgref)

    for lyotradius in lyotlist:

        lyotmask = makegmtpupil(rseg=lyotradius)

        bkg = (bkgconv * lyotmask).sum()
        sig = (sigconv * lyotmask).sum()
        snr = sig / sqrt(bkg)

        print "%s\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.4f" % ( m2type, fldang, lyotradius, sig/sigref, bkg/bkgref, snr/snrref )
