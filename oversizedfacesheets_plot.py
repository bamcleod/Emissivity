import starbase
import pylab


color = {
	"base"      : "red",
	"bigm2"     : "blue",
	"bigref"    : "green",
    "cyl"       : "orange",
    "tallcyl"   : "cyan",
    "filledcyl" : "yellow",
    "bigfilledcyl": "brown"
	}

legend = {
	"base"      : "Nominal ASM (no edge sensors)",
	"bigm2"     : "Round oversized M2",
	"bigref"    : "Round ref body",
    "cyl"       : "Cylindrical ref body",
    "tallcyl"   : "Cyl ref body + \"red\" volume",
    "filledcyl" : "Monolithic ref body",
    "bigfilledcyl": "Monolithic ref, oversized segments"
	}

ltype = {
	"0.5": "-",
	"10":  "--"
	}



for config in ["filledcyl", "bigfilledcyl"] :
	for slit in ["10", "0.5"]:
		data=starbase.Starbase(config+"_s"+slit+"_t0.2_o.1.tab")
		if slit=="0.5":
			pylab.plot(data.fldang,data.snratio,ltype[slit],color=color[config],label=legend[config])
		else:
			pylab.plot(data.fldang,data.snratio,ltype[slit],color=color[config])
				
pylab.legend(loc=3)
pylab.ylabel("S/N")
pylab.xlabel("Field angle (arcmin)")
pylab.xlim([0,5])
pylab.savefig("oversizedfacesheets.png")
pylab.show()

