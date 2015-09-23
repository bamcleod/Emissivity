import starbase
import pylab


color = {
    "base"      : "red",
    "bigm2"     : "blue",
    "bigref"    : "green",
    "cyl"       : "orange",
    "tallcyl"   : "cyan",
    "filledcyl" : "yellow"
    }

legend = {
    "base"      : "Nominal ASM (no edge sensors)",
    "bigm2"     : "Round oversized M2",
    "bigref"    : "Round ref body",
    "cyl"       : "Cylindrical ref body",
    "tallcyl"   : "Cyl ref body + \"red\" volume",
    "filledcyl" : "Monolithic ref body"
    }

ltype = {
    "0.5": "-",
    "10":  "--"
    }


norm = None
pylab.figure()
for subplot,telbkg in enumerate(["0.05", "0.1", "0.2"]):
    pylab.subplot(2,3,subplot+1)
    for config in ["base", "bigm2", "bigref"] :
        for slit in ["10", "0.5"]:
            data=starbase.Starbase(config+"_tel_s"+slit+"_t"+telbkg+"_o0.03_nonorm.tab")
            if norm==None:
                norm = data.snratio[0];
            if slit=="0.5" and telbkg=="0.05":
                pylab.plot(data.fldang,data.snratio/norm,ltype[slit],color=color[config],label=legend[config])
            else:
                pylab.plot(data.fldang,data.snratio/norm,ltype[slit],color=color[config])
                
    #pylab.legend(loc=3)
    if (subplot==0):
        pylab.ylabel("S/N")
    else:
        pylab.setp(pylab.gca(), yticklabels=[])
    pylab.xlabel("Field angle (arcmin)")
    pylab.title("Emissivity="+telbkg)
    pylab.xlim([0,5])
    pylab.ylim([0.5,1])
    #    pylab.ylim([0,1700])
    pylab.savefig("tel_threecases_nonorm.png")

pylab.figure()
for subplot,telbkg in enumerate(["0.05", "0.1", "0.2"]):
    pylab.subplot(2,3,subplot+1)
    for config in ["base", "cyl", "tallcyl", "filledcyl"] :
        for slit in ["10", "0.5"]:
            data=starbase.Starbase(config+"_tel_s"+slit+"_t"+telbkg+"_o0.03_nonorm.tab")
            if slit=="0.5" and telbkg=="0.05":
                pylab.plot(data.fldang,data.snratio/norm,ltype[slit],color=color[config],label=legend[config])
            else:
                pylab.plot(data.fldang,data.snratio/norm,ltype[slit],color=color[config])
                
    if (subplot==0):
        pylab.ylabel("S/N")
    else:
        pylab.setp(pylab.gca(), yticklabels=[])
    pylab.xlabel("Field angle (arcmin)")
    pylab.title("Emissivity="+telbkg)
    pylab.xlim([0,5])
    pylab.ylim([0.5,1])
    #    pylab.ylim([0,1700])
pylab.savefig("tel_adopticaperf_nonorm.png")

pylab.show()

