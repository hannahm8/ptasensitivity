import numpy as np
import matplotlib.pyplot as plt
import hasasia.sensitivity as hsen
import hasasia.sim as hsim

import sys
sys.path.append('../scalingRelations/')
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

def constructPTA(sig,radRA,radDec,freqs):
    T = 10 # in years 
    c = 26

    psrs = hsim.sim_pta(timespan=T, cad=c, sigma=sig,phi=radRA,theta=radDEC)
    spectra = []
    for p in psrs:
        sp = hsen.Spectrum(p, freqs=freqs)
        sp.NcalInv
        spectra.append(sp)

    scGWB = hsen.GWBSensitivityCurve(spectra)

    return scGWB


# original t obs and sigmas 
psrDataFile = '../data/psrDetails.dat'
psrNames, \
psrObsConstants, \
psrStartingObsTimes, \
angles, \
hdValues = readInData.readDataIntoDicts(psrDataFile)

# read in again to get positions and sigmas
data = np.genfromtxt(psrDataFile,names=True)
radRA  = [  (ra*np.pi)/180. for  ra in data['RA']  ]
radDEC = [ (dec*np.pi)/180. for dec in data['DEC'] ]
sigmas = data['ExpPrecision']*1.E-6

# get the new times
#shuffleFile = '../scalingRelations/halfTimeShuffle/shuffle_3.dat'
#shuffleFile = '../../runs/shuffle/quarterTimeShuffle/times/shuffle_165.dat'
shuffleFile = '../../runs/shuffle/quarterTimeShuffleDPHD/times/shuffle_87.dat'
shuffleNames = np.genfromtxt(shuffleFile,usecols=0,dtype=str)
shuffleTobs = np.genfromtxt(shuffleFile,usecols=1)
# put in dictionary
shuffleTimes = {}
for ipsr,time in zip(shuffleNames, shuffleTobs):    
    shuffleTimes[ipsr] = time

# work out the new sigmas
sigmaShuffle = np.zeros(len(psrNames))
print(sigmaShuffle)
for i,ipsr in enumerate(psrNames):
    sigmaShuffle[i] = psrObsConstants[ipsr] / np.sqrt(shuffleTimes[ipsr])

print((sigmaShuffle-sigmas))
print(sum(sigmaShuffle-sigmas))
print('original', sigmas)
print('shuffle', sigmaShuffle)


# for plotting 
freqs = np.logspace(np.log10(5e-10),np.log10(5e-7),500)


scGWBOriginal = constructPTA(sigmas, radRA, radDEC, freqs)
plt.loglog(freqs,scGWBOriginal.h_c,label='original')

scGWBShuffle = constructPTA(sigmaShuffle, radRA, radDEC, freqs)
plt.loglog(freqs,scGWBShuffle.h_c,label='shuffle',ls=':')


#plt.loglog(spectra[0].freqs,spectra[0].h_c)
plt.xlabel('Frequency [Hz]')
plt.ylabel('Characteristic Strain, $h_c$')
plt.legend()
plt.savefig('hasasia_shuffle_with_dipole_hd_diff.png')
plt.show()


