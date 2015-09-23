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



for config in ["base", "bigm2", "bigref"] :
    for slit in ["10", "0.5"]:
        for fld in ["0.0", "2.5", "5.0"]:
            data=starbase.Starbase(config+"_lyot_s"+slit+"_t0.2_o.1_f"+fld+".tab")
            if slit=="0.5" and fld=="0.0":
                pylab.plot(data.lyotrad,data.snratio,ltype[slit],color=color[config],label=legend[config])
            else:
                pylab.plot(data.lyotrad,data.snratio,ltype[slit],color=color[config])
				
pylab.legend(loc=3)
pylab.ylabel("S/N")
pylab.xlabel("Lyot stop segment radius (m)")
pylab.show()
pylab.savefig("lyotrad_plot.png")

