import pyfits
h=pyfits.open("p.fits")
d=h[0].data

def circle (npix, r0, x0=0, y0=0, scale=1):
    y = yindx - npix/2 - y0 * scale
    x = xindx - npix/2 - x0 * scale
    r = sqrt(x*x+y*y)
    c=1-(r-(r0*scale-0.5))
    c[c<0] = 0
    c[c>1] = 1
    return c

from numpy import *
npix=866
pupilscale=1

npix=886
yindx,xindx=indices((npix,npix),dtype=float32)
c=circle(0,154,265,656)
c.sum()
(c*d).sum()

import pylab
pylab.imshow(c*d)

e=circle(0,154,265,229)
pylab.imshow(e*d)
(e*d).sum()
(e*d).sum()/(c*d).sum()

b=circle(0,154,689,656)
pylab.imshow(b*d)
(b*d).sum()/(c*d).sum()

import readline
readline.write_history_file("hist.txt")
