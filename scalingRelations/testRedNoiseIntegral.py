import numpy as np
import matplotlib.pyplot as plt
import readInData

import snrFunctions

def wPSD(sigma,c):  
    deltat = 1./c
    return 2.*sigma*deltat

def rPSD(f,redAmp,gamma,fref):
    nSecondsInYear = 365.25*24.*60.*60.
    a = (redAmp*redAmp) / (12.*np.pi*np.pi) * nSecondsInYear**3.
    r = a * (f/fref)**-gamma 
    return r

def gwPSD(f,A,alpha,beta,fref):
    b = snrFunctions.get_b(A,fref,alpha)
    return b*(f**-beta)



def integrand(f,c,fref,sigmaI,redAI,gammaI,sigmaJ,redAJ,gammaJ,A,alpha,beta):
    # white noise
    deltat = 1./c 
    wI = 2.*sigmaI*deltat
    wJ = 2.*sigmaJ*deltat
    
    # red noise
    nSecondsInYear = 365.25*24.*60.*60.
    aI = (redAI*redAI) / (12.*np.pi*np.pi) * nSecondsInYear**3.
    rI = aI * (f/fref)**-gammaI
    aJ = (redAJ*redAJ) / (12.*np.pi*np.pi) * nSecondsInYear**3.
    rJ = aJ * (f/fref)**-gammaJ
    
    # gw signal
    b = snrFunctions.get_b(A,fref,alpha)
    gw = b*f**-beta

    # overall 
    I = (gw*gw) / ((wI+rI+gw)*(wJ+rJ+gw))
    return I

def integrandEqualSigma(f,sig,deltat,beta,A,fref,alpha):  
    b = snrFunctions.get_b(A,fref,alpha)
    I = 1./(1.+(sig*sig*deltat)/(b*f**-beta))**2.
    return I 


psrDataFile = '../data/psrDetails.dat'
#psrDataFile = '/home/hannahm/repositories/ptasensitivity/data/trialPSRData.dat'
dataOriginalFormat = np.genfromtxt(psrDataFile, names=True)

psrNames, \
psrObsConstants, \
psrStartingObsTimes, \
angles, \
angCorrValues = readInData.readDataIntoDicts(psrDataFile)


redAmps = psrObsConstants.copy()
# update
for ipsr in psrNames:
   redAmps[ipsr] = 0

gammas = psrObsConstants.copy()
for ipsr in psrNames: 
    gammas[ipsr] = 1



oneYearInSeconds = (365.25*24.*60.*60.)

T = 11.
TInSeconds = T * oneYearInSeconds

A = 2.E-15
beta = 13./3
alpha = (3.-beta)/2.
fref = 1./oneYearInSeconds

c = 26./oneYearInSeconds


#def get_integral_with_red_noise(c,fref,sigI,rAI,gamI,sigJ,rAJ,gamJ,A,alpha):


avPSD  = snrFunctions.avePTASNR(psrNames,psrObsConstants,angCorrValues,psrStartingObsTimes,A,alpha,beta,fref,TInSeconds,c)


avPSDR = snrFunctions.avePSD_incRedNoise(psrNames,psrObsConstants,angCorrValues,psrStartingObsTimes,\
                            redAmps,gammas,A,alpha,beta,fref,TInSeconds,c)

print('w only', avPSD)
print('w+r', avPSDR)

exit()

fL=1./TInSeconds
fH=0.5*c
deltat = 1./c


freqs = np.linspace(fL,fH,100)
w = np.zeros(len(freqs))
r = np.zeros(len(freqs))
g = np.zeros(len(freqs))

sigmaI,redAI,gammaI = 1E-6, 0, 4.
sigmaJ,redAJ,gammaJ = 1E-6, 0, 2.

sig=1E-6
redA=1E-12
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

