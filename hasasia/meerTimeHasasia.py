import numpy as np
import matplotlib.pyplot as plt
import hasasia.sensitivity as hsen
import hasasia.sim as hsim

import sys
sys.path.append('/home/hannahm/repositories/ptasensitivity/snr/')
sys.path.append('/fred/oz005/users/hmiddlet/ptasensitivity/snr/')
import readInData

"""
import snrFunctions
import readInData

# data file
psrDataFile = '../data/psrDetails.dat'

# read in data and compute angles etc
psrNames, \
psrObsConstants, \
psrStartingObsTimes, \
angles, \
hdValues = readInData.readDataIntoDicts(psrDataFile)

"""

def constructPTA(sig,radRA,radDec,freqs,redNoise=False,redA=None,redGamma=None,c=26):
    T = 10 # in years 
    
    if redNoise==True: 
        print(redA)
        print(redGamma)
        psrs = hsim.sim_pta(timespan=T,\
                            cad=c,\
                            sigma=sig,\
                            phi=radRA,\
                            theta=radDEC,\
                            A_rn=redA,\
                            alpha=redGamma,\
                            freqs=freqs)

    elif redNoise==False:
        print('no red noise or gamma provided')      
        psrs = hsim.sim_pta(timespan=T, cad=c, sigma=sig, phi=radRA, theta=radDEC, freqs=freqs)
    

    spectra = []
    for p in psrs:
        sp = hsen.Spectrum(p, freqs=freqs)
        sp.NcalInv
        spectra.append(sp)

    scGWB = hsen.GWBSensitivityCurve(spectra)

    return scGWB

def convertShuffleData(shuffleFile,psrNames,obsConstants):

    shuffleNames = np.genfromtxt(shuffleFile,usecols=0,dtype=str)
    shuffleTobs = np.genfromtxt(shuffleFile,usecols=1)

    # put in dictionary
    shuffleTimes = {}
    for ipsr,time in zip(shuffleNames, shuffleTobs):    
        shuffleTimes[ipsr] = time

    # work out the new sigmas
    sigmaShuffle = np.zeros(len(psrNames))
    for i,ipsr in enumerate(psrNames):
        sigmaShuffle[i] = psrObsConstants[ipsr] / np.sqrt(shuffleTimes[ipsr])

    return sigmaShuffle


# original t obs and sigmas 
psrDataFile = '/fred/oz005/users/hmiddlet/ptasensitivity/data/psrDetails.dat'
redNoiseFile = '/fred/oz005/users/hmiddlet/ptasensitivity/data/redNoise.dat'

psrNames, \
psrObsConstants, \
psrStartingObsTimes, \
angles, \
hdValues, \
ampRed, gammaRed, \
_ = readInData.readDataIntoDicts(psrDataFile, 'HD', \
                                 redNoiseFile=redNoiseFile)


# read in again to get positions and sigmas
data = np.genfromtxt(psrDataFile,names=True)
radRA  = [  (ra*np.pi)/180. for  ra in data['RA']  ]
radDEC = [ (dec*np.pi)/180. for dec in data['DEC'] ]
sigmas = data['ExpPrecision']*1.E-6


# put red noise into arrays
rA, rG = np.zeros(len(psrNames)), np.zeros(len(psrNames))
for i,ipsr in enumerate(psrNames):
    rA[i] = ampRed[ipsr]
    rG[i] = gammaRed[ipsr]





# get the new times
shuffleFileHD     = './eighthTimeHDMin256/times/shuffle_133.dat'
sigmaShuffleHD    = convertShuffleData(shuffleFileHD,psrNames,psrObsConstants)

shuffleFileDPHD   = './eighthTimeDPHDDiffMin256/times/shuffle_172.dat'
sigmaShuffleDPHD  = convertShuffleData(shuffleFileDPHD,psrNames,psrObsConstants)

shuffleFileEqual  = './eighthTimeEqualMin256/times/shuffle_136.dat'
sigmaShuffleEqual = convertShuffleData(shuffleFileEqual,psrNames,psrObsConstants)

#linestyles = ['solid','dotted','dashed','dashdot']
# for plotting 
freqs = np.logspace(np.log10(5e-10),np.log10(5e-7),500)


plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(figsize=[10,8])

scGWBOriginal = constructPTA(sigmas, radRA, radDEC, freqs, redNoise=True, redA=rA, redGamma=rG)
ax.loglog(freqs,scGWBOriginal.h_c,label='Original',ls='solid')

scGWBShuffleHD = constructPTA(sigmaShuffleHD, \
                              radRA, radDEC, \
                              freqs, redNoise=True, \
                              redA=rA, redGamma=rG)
ax.loglog(freqs,scGWBShuffleHD.h_c,label='HD',ls='dotted')

scGWBShuffleDPHD = constructPTA(sigmaShuffleDPHD, \
                                radRA, radDEC, \
                                freqs, redNoise=True, \
                                redA=rA, redGamma=rG)
ax.loglog(freqs,scGWBShuffleDPHD.h_c,label='DP-2HD',ls='dashed')

scGWBShuffleEqual = constructPTA(sigmaShuffleEqual, \
                                 radRA, radDEC, \
                                 freqs, redNoise=True, \
                                 redA=rA, redGamma=rG)
ax.loglog(freqs,scGWBShuffleEqual.h_c,label='Equal',ls='dashdot')


# inset
#from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)

plt.rcParams.update({'font.size': 12})
ax2 = plt.axes([0,1,1,1])
ip = InsetPosition(ax,[0.14,0.4,0.4,0.55])
ax2.set_axes_locator(ip)
mark_inset(ax,ax2,loc1=2,loc2=4,fc="none",ec="0.5")
ax2.loglog(freqs,scGWBOriginal.h_c,ls='solid')
ax2.loglog(freqs,scGWBShuffleHD.h_c,ls='dotted')
ax2.loglog(freqs,scGWBShuffleDPHD.h_c,ls='dashed')
ax2.loglog(freqs,scGWBShuffleEqual.h_c,ls='dashdot')
ax2.set_xlim(2E-9,7E-9)
ax2.set_ylim(4E-16,6.5E-16)


plt.rcParams.update({'font.size': 15})
#plt.loglog(spectra[0].freqs,spectra[0].h_c)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Characteristic Strain, $h_c$')
ax.legend()
plt.tight_layout()
plt.savefig('eighthsCompare/hasasiaShuffle.png')
plt.savefig('eighthsCompare/hasasiaShuffle.pdf')
plt.show()


