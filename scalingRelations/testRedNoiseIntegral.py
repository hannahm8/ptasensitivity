import numpy as np
import matplotlib.pyplot as plt
import readInData

import snrFunctions

import time


def wPSD(sigma,c):  
    deltat = 1./c
    return 2.*sigma*sigma*deltat

def rPSD(f,redAmp,gamma,fref):
    nSecondsInYear = 365.25*24.*60.*60.
    a = (redAmp*redAmp) / (12.*np.pi*np.pi) * nSecondsInYear**3.
    r = a * (f/fref)**-gamma 
    return r

def gwPSD(f,A,alpha,beta,fref):
    b = snrFunctions.get_b(A,fref,alpha)
    return b*(f**-beta)





psrDataFile = '../data/psrDetails.dat'
#psrDataFile = '/home/hannahm/repositories/ptasensitivity/data/trialPSRData.dat'
dataOriginalFormat = np.genfromtxt(psrDataFile, names=True)

# path to red noise file
redNoiseData = '../data/redNoise.dat'

psrNames, \
psrObsConstants, \
psrStartingObsTimes, \
angles, \
angCorrValues, \
redAmps, redGammas = readInData.readDataIntoDicts(psrDataFile,redNoiseFile=redNoiseData)


#redAmps = psrObsConstants.copy()


###### read in red noise 




oneYearInSeconds = (365.25*24.*60.*60.)

T = 11.
TInSeconds = T * oneYearInSeconds

A = 2.E-15
beta = 13./3
alpha = (3.-beta)/2.
fref = 1./oneYearInSeconds

c = 26./oneYearInSeconds

fL=1./TInSeconds
fH=0.5*c
deltat = 1./c
freqs = np.linspace(fL,fH,100)

#def get_integral_with_red_noise(c,fref,sigI,rAI,gamI,sigJ,rAJ,gamJ,A,alpha):

start_time = time.time()

Ts = np.linspace(1,12,10)
avPSD, avPSDR = np.zeros(len(Ts)), np.zeros(len(Ts))
for i,Ti in enumerate(Ts):

    TInSeconds = Ti * oneYearInSeconds
    avPSD[i]  = snrFunctions.avePTASNR(psrNames,psrObsConstants,angCorrValues,psrStartingObsTimes,A,alpha,beta,fref,TInSeconds,c)


    avPSDR[i] = snrFunctions.avePTASNR_incRedNoise(psrNames,psrObsConstants,angCorrValues,psrStartingObsTimes,\
                            redAmps,redGammas,A,alpha,beta,fref,TInSeconds,c)



end_time = time.time()

print('time was ', end_time-start_time)

print('plotting')
plt.plot(Ts,avPSD,color='b',label='original')
plt.plot(Ts,avPSDR,color='r',label='withred')
plt.legend()
plt.show()
plt.clf()
print('w only', avPSD)
print('w+r', avPSDR)



"""


freqs = np.linspace(fL,fH,100)
w = np.zeros(len(freqs))
r = np.zeros(len(freqs))
g = np.zeros(len(freqs))

sigmaI,redAI,gammaI = 1E-6, 0, 4.
sigmaJ,redAJ,gammaJ = 1E-6, 0, 2.

sig=1E-6
redA=1E-14
gamma = 3.

for i,f in enumerate(freqs):
    w[i] = wPSD(sig,c)
    r[i] = rPSD(f,redA,gamma,fref)
    g[i] = gwPSD(f,A,alpha,beta,fref)
    
total = w+g+r

plt.plot(freqs,w,label='w')
plt.plot(freqs,r,label='r')
plt.plot(freqs,g,label='g')
plt.plot(freqs,total,label='total')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()
"""

