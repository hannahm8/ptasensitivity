import numpy as np
import snrFunctions



# red in th deat
psrDataFile = '../data/psrDetails.dat'
psrNames = np.genfromtxt(psrDataFile,usecols=0,dtype=str)
psrData = np.genfromtxt(psrDataFile,names=True)



# compute obs constants 
obsConstants = [ sig * np.sqrt(intT) for sig, intT in \
                                     zip (psrData['ExpPrecision']*1.E-6, \
                                          psrData['IntTime'])]
psrObsConstants = {}
psrStartingObsTimes = {}
for psr, oc, st in zip(psrNames, obsConstants, psrData['IntTime']):
    psrObsConstants[psr] = float(oc)
    psrStartingObsTimes[psr] = float(st)


# work out angles and hd values ahead of time
angles = {}
hdValues = {}
for i, ipsr in enumerate(psrNames):
    onePSRAng = {}
    onePSRHDs = {}
    for j, jpsr in enumerate(psrNames):
        if jpsr==ipsr:
            angle = 0
            hd = None
        else: 
            elati, elongi = psrData['ELAT'][i], psrData['ELONG'][i]
            elatj, elongj = psrData['ELAT'][j], psrData['ELONG'][j]
            angle = snrFunctions.h2(elati,elatj,elongi,elongj)
            hd = snrFunctions.hellings_downs(angle)
        onePSRAng[jpsr] = angle
        onePSRHDs[jpsr] = hd
    angles[ipsr]   = onePSRAng
    hdValues[ipsr] = onePSRHDs


oneYearInSeconds = (365.25*24.*60.*60.)

Ts = np.linspace(0,15,50)
TInSeconds = Ts * oneYearInSeconds

A = 2.E-15
beta = 13./3
alpha = (3.-beta)/2.
fref = 1./oneYearInSeconds

c = 26./oneYearInSeconds



snr = np.zeros(len(Ts))
for i, Ti in enumerate(TInSeconds):
    snr[i] = snrFunctions.avePTASNR(psrNames,\
                                    psrObsConstants,\
                                    hdValues,\
                                    psrStartingObsTimes,\
                                    A,alpha,beta,fref,Ti,c)



import matplotlib.pyplot as plt
plt.plot(Ts,snr)
plt.grid()
plt.show()

