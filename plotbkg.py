import starbase
import pylab

t=starbase.Starbase("phasingbkg.tab")

pylab.plot(t.wave_Bin,t.InPupil,label="In pupil")
pylab.plot(t.wave_Bin,t.OutPupil,label="Out pupil")
pylab.xlabel("Wavelength [nm]")
pylab.ylabel(r"Background [$photons/sec/arcsec^2/nm/m^2$]")
pylab.title("Phasing camera background components")
pylab.legend(loc=2)
pylab.savefig("phasebkg.png")

pylab.figure()

pylab.plot(t.wave_Bin,t.sky_Mean,label="Sky")
pylab.plot(t.wave_Bin,t.bb_Mean,label="Blackbody (284.75K)")
pylab.xlabel("Wavelength [nm]")
pylab.ylabel(r"Background [$photons/sec/arcsec^2/nm/m^2$]")
pylab.title("Sky and thermal background")
pylab.legend(loc=2)
pylab.savefig("skythermbkg.png")



