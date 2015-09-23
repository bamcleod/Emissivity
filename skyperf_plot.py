import starbase
import pylab
import numpy

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

fld="2.5"
for config in ["base", "bigm2", "bigref"] :
    for slit in ["10", "0.5"]:
        skybrt=[]
        sn=[]
        for sky in ["0.01", "0.03", "0.1", "0.3", "1", "3", "10"]:
            data=starbase.Starbase(config+"_sky_s"+slit+"_t0.2_o"+sky+"_f"+fld+".tab")
            skybrt.append(float(sky))
            sn.append(data.snratio[0])
        skybrt=pylab.log(skybrt)/pylab.log(10.)
        if slit=="0.5" and fld=="2.5":
            pylab.plot(skybrt,sn,ltype[slit],color=color[config],label=legend[config])
        else:
            pylab.plot(skybrt,sn,ltype[slit],color=color[config])
				
pylab.legend(loc=3)
pylab.ylabel("S/N")
pylab.xlabel("log(Sky effective emissivity)")

bars=[
    ["H", 1.04, -0.3,  1.5 ],
    ["K", 1.03, -2,   -0.3 ],
    ["L", 1.02, -1.5, -0.5 ],
    ["M", 1.01, -1.0,  0.0 ]
    ]

left=[]
height=[]
width=[]
bottom=[]
h=0.01
for label, y, xmin, xmax in bars:
    left.append(xmin)
    height.append(h)
    width.append(xmax-xmin)
    bottom.append(y-h/2)
    pylab.bar(left,height,width,bottom=bottom,color="yellow")
    pylab.text((xmin+xmax)/2, y, label, ha='center', va='center')

pylab.savefig("skyperf_plot.png")
pylab.show()

