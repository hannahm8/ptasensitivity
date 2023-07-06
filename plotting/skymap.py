
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic 
import sys
sys.path.append('/home/hannahm/repositories/ptasensitivity/snr/')
sys.path.append('/fred/oz005/users/hmiddlet/ptasensitivity/snr/')




def galacticPlane(shift):
    nsamp = 1000
    #glat = np.zeroes((1, nsamp))
    glong = np.linspace(-180, 180, nsamp)
    ra = []
    dec = []
    ra2 = []
    dec2 = []
    for i in range(0, nsamp):
        l = glong[i]
        c = SkyCoord(l, -5, frame=Galactic,
                     unit=(u.deg, u.deg))
        ra.append(c.icrs.ra.value)
        dec.append(c.icrs.dec.value)
        c2 = SkyCoord(l, 5, frame=Galactic,
                     unit=(u.deg, u.deg))
        ra2.append(c2.icrs.ra.value)
        dec2.append(c2.icrs.dec.value)


    ra = np.array(ra)*np.pi/180 + shiftTo18h
    ra[(ra > np.pi)] = ra[(ra > np.pi)] - 2*np.pi
    dec = np.array(dec)*np.pi/180

    ra2 = np.array(ra2)*np.pi/180 + shiftTo18h
    ra2[(ra2 > np.pi)] = ra2[(ra2 > np.pi)] - 2*np.pi
    dec2 = np.array(dec2)*np.pi/180


    ra_resamp = np.linspace(min(ra), max(ra), nsamp)
    index = np.argsort(ra)
    dec_resamp = np.interp(ra_resamp, ra[index], dec[index])
    index = np.argsort(ra2)
    dec2_resamp = np.interp(ra_resamp, ra2[index], dec2[index])

    plt.fill_between(-ra_resamp, dec_resamp, dec2_resamp, zorder=0, alpha=0.3, color='k')

    ra = []
    dec = []
    ra2 = []
    dec2 = []
    for i in range(0, nsamp):
        l = glong[i]
        c = SkyCoord(l, -10, frame=Galactic,
                     unit=(u.deg, u.deg))
        ra.append(c.icrs.ra.value)
        dec.append(c.icrs.dec.value)
        c2 = SkyCoord(l, 10, frame=Galactic,
                     unit=(u.deg, u.deg))
        ra2.append(c2.icrs.ra.value)
        dec2.append(c2.icrs.dec.value)


    ra = np.array(ra)*np.pi/180 + shiftTo18h
    ra[(ra > np.pi)] = ra[(ra > np.pi)] - 2*np.pi
    dec = np.array(dec)*np.pi/180 

    ra2 = np.array(ra2)*np.pi/180 + shiftTo18h
    ra2[(ra2 > np.pi)] = ra2[(ra2 > np.pi)] - 2*np.pi
    dec2 = np.array(dec2)*np.pi/180


    ra_resamp = np.linspace(min(ra), max(ra), nsamp)
    index = np.argsort(ra)
    dec_resamp = np.interp(ra_resamp, ra[index], dec[index])
    index = np.argsort(ra2)
    dec2_resamp = np.interp(ra_resamp, ra2[index], dec2[index])

    plt.fill_between(-ra_resamp, dec_resamp, dec2_resamp, zorder=0, alpha=0.1, color='k')

    return 


def convert(ra_deg, dec_deg):

    position = SkyCoord(ra_deg*u.deg, dec_deg*u.deg)

    return position.galactic.l.rad, position.galactic.b.rad



psrDataFile = sys.argv[1]
shuffleFile = sys.argv[2]
shuffleLabel = sys.argv[3]

#psrDataFile = '/home/hannahm/repositories/ptasensitivity/data/psrDetails.dat'
data = np.genfromtxt(psrDataFile,names=True)
startTobs = data['IntTime']


# get shuffle times
#shuffleFile = './shuffle_133.dat'
shuffleTobs = np.genfromtxt(shuffleFile,usecols=1)

nPSRs = len(data['RA'])

# size depends on how much time shifted
timeDiff = (shuffleTobs - startTobs)/10

# set up for plot colours
c = []
for i in range(nPSRs):
    if timeDiff[i]>=0: c.append('#F5793A')
    elif timeDiff[i]<0: c.append('#85C8F9')



fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111,projection="mollweide")
ax.grid(True)


# we have RA and dec in degrees
shiftTo18h = +np.pi/2.#-np.pi/2

ras, decs = np.zeros(nPSRs),np.zeros(nPSRs)
for i in range(nPSRs):

    ras[i]  = data['RA'][i] * (np.pi/180) + shiftTo18h
    decs[i] = data['DEC'][i] * (np.pi/180)
    



GC = SkyCoord(0,0,frame=Galactic,unit=(u.deg,u.deg))
raGC = GC.icrs.ra.value *np.pi/180 + shiftTo18h
decGC = GC.icrs.dec.value  *np.pi/180
if raGC > np.pi:
    raGC = raGC - 2*np.pi
plt.scatter(-raGC, decGC, color='yellow', s=200, marker='*')

for i,ra in enumerate(ras):
    if ra > np.pi:
        ras[i] = ra - 2*np.pi

"""
#plt.scatter(0,0)
raPick,decPick = 107.97594058,-68.51322193
testRA = raPick * np.pi/180
if testRA > np.pi: testRA = testRA - 2*np.pi
testDec = decPick * np.pi/180
plt.scatter(testRA,testDec,color='r',marker='x')

"""

galacticPlane(shiftTo18h)

# plot all pulsar positions
plt.scatter(-ras, decs, color='k',marker='x',alpha=0.5)
plt.scatter(-ras, decs, s=abs(timeDiff), alpha=0.9,color=c)

ax.set_xticklabels(['4h','2h','0h','22h','20h','18h','16h','14h','12h','10h','8h'])

plt.tight_layout()

plt.savefig('skymap_{}.png'.format(shuffleLabel))
plt.savefig('skymap_{}.pdf'.format(shuffleLabel))
plt.show()
exit()



