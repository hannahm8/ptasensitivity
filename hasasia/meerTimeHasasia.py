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

def constructPTA(sig,radRA,radDec,freqs,redNoise=False,redA=None,redGamma=None):
    T = 10 # in years 
    c = 26
    
    if redNoise==True: 
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
psrDataFile = '../data/psrDetails.dat'
redNoiseFile = '../data/redNoise.dat'

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

shuffleFileDPHD   = './eighthTimeDPHDMin256/times/shuffle_172.dat'
sigmaShuffleDPHD  = convertShuffleData(shuffleFileDPHD,psrNames,psrObsConstants)

shuffleFileEqual  = './eigthTimeEqual/Min256/times/shuffle_136.dat'
sigmaShuffleEqual = convertShuffleData(shuffleFileEqual,psrNames,psrObsConstants)

#linestyles = ['solid','dotted','dashed','dashdot']
# for plotting 
freqs = np.logspace(np.log10(5e-10),np.log10(5e-7),500)

scGWBOriginal = constructPTA(sigmas, radRA, radDEC, freqs, redNoise=True, redA=rA, redGamma=rG)
plt.loglog(freqs,scGWBOriginal.h_c,label='original',ls='solid')

scGWBShuffle = constructPTA(sigmaShuffleHD, radRA, radDEC, freqs, redNoise=True, redA=rA, redGamma=rG)
plt.loglog(freqs,scGWBShuffle.h_c,label='HD',ls='dotted')

scGWBShuffle = constructPTA(sigmaShuffleDPHD, radRA, radDEC, freqs, redNoise=True, redA=rA, redGamma=rG)
plt.loglog(freqs,scGWBShuffle.h_c,label='DP-2HD',ls='dashed')

scGWBShuffle = constructPTA(sigmaShuffleEqual, radRA, radDEC, freqs, redNoise=True, redA=rA, redGamma=rG)
plt.loglog(freqs,scGWBShuffle.h_c,label='Equal',ls='dashdot')

#plt.loglog(spectra[0].freqs,spectra[0].h_c)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Characteristic Strain, $h_c$')
plt.legend()
#plt.savefig('hasasia_shuffle_with_dipole_hd_diff.png')
plt.show()


