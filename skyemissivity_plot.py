import starbase
import pylab as plt
d=starbase.Starbase("atmo_tel.results")

w=d.wave.reshape((-1,10)).mean(axis=1)
s=d.skyem.reshape((-1,10)).mean(axis=1)


plt.suptitle("Effective Emissivity of Sky")

configs=[
    [ 1, "H", 1400, 1900, -1, 3],
    [ 2, "K", 1900, 2500, -2.5, 0.5],
    [ 3, "L", 2800, 4000, -2, 0],
    [ 4, "M", 4400, 5500, -2, 0]
    ]
    
for isub,label,xmn,xmx,ymn,ymx in configs:
    plt.subplot(2,2,isub)
    plt.plot(w,plt.log(s)/plt.log(10.))
    plt.xlim([xmn,xmx])
    plt.ylim([ymn,ymx])
    plt.text(xmn+50,ymn+0.25, label)
    if isub>2:
        plt.xlabel("Wavelength (nm)")
    if isub%2==1:
        plt.ylabel("log(Emissivity)")

plt.savefig("sky_emissivity.png")

