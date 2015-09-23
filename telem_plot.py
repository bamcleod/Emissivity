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



pylab.figure()
for subplot,telbkg in enumerate(["0.05", "0.1", "0.2"]):
    pylab.subplot(2,3,subplot+1)
    for config in ["base", "bigm2", "bigref"] :
        for slit in ["10", "0.5"]:
            data=starbase.Starbase(config+"_tel_s"+slit+"_t"+telbkg+"_o0.03.tab")
            if slit=="0.5" and telbkg=="0.05":
                pylab.plot(data.fldang,data.snratio,ltype[slit],color=color[config],label=legend[config])
            else:
                pylab.plot(data.fldang,data.snratio,ltype[slit],color=color[config])
                
    #pylab.legend(loc=3)
    if (subplot==0):
        pylab.ylabel("S/N")
    else:
        pylab.setp(pylab.gca(), yticklabels=[])
    pylab.xlabel("Field angle (arcmin)")
    pylab.title("Emissivity="+telbkg)
    pylab.xlim([0,5])
    pylab.ylim([0.6,1])
    pylab.savefig("tel_threecases.png")

pylab.figure()
for subplot,telbkg in enumerate(["0.05", "0.1", "0.2"]):
    pylab.subplot(2,3,subplot+1)
    for config in ["base", "cyl", "tallcyl", "filledcyl"] :
        for slit in ["10", "0.5"]:
            data=starbase.Starbase(config+"_tel_s"+slit+"_t"+telbkg+"_o0.03.tab")
            if slit=="0.5" and telbkg=="0.05":
                pylab.plot(data.fldang,data.snratio,ltype[slit],color=color[config],label=legend[config])
            else:
                pylab.plot(data.fldang,data.snratio,ltype[slit],color=color[config])
                
    if (subplot==0):
        pylab.ylabel("S/N")
    else:
        pylab.setp(pylab.gca(), yticklabels=[])
    pylab.xlabel("Field angle (arcmin)")
    pylab.title("Emissivity="+telbkg)
    pylab.xlim([0,5])
    pylab.ylim([0.6,1])
    pylab.savefig("tel_adopticaperf.png")

    pylab.show()

