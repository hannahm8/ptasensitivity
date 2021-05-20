import numpy as np
import matplotlib.pyplot as plt

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

# workgin out total time used 
totalTime=0
for ipsr in psrNames:
    totalTime+=psrStartingObsTimes[ipsr]


oneYearInSeconds = (365.25*24.*60.*60.)

T = np.linspace(1., 11., 50)
TInSeconds = T * oneYearInSeconds

A = 2.E-15
beta = 13./3
alpha = (3.-beta)/2.

fref = 1./oneYearInSeconds

c = 26./oneYearInSeconds

snr = np.zeros(len(T))
for i, Ti in enumerate(TInSeconds):

    snr[i] = snrFunctions.avePTASNR(psrNames,\
                                  psrObsConstants,\
                                  hdValues,\
                                  psrStartingObsTimes,\
                                  A,alpha,beta,fref,Ti,c)
plt.plot(T,snr,label='Tobs in spreadsheet')

##########################################################################
# read in data and compute angles etc
psrNames, \
psrObsConstants, \
psrStartingObsTimes, \
angles, \
hdValues = readInData.readDataIntoDicts('timeShuffleResults1.dat')

# workgin out total time used 
totalTime=0
for ipsr in psrNames:
    totalTime+=psrStartingObsTimes[ipsr]

snr = np.zeros(len(T))
for i, Ti in enumerate(TInSeconds):

    snr[i] = snrFunctions.avePTASNR(psrNames,\
                                  psrObsConstants,\
                                  hdValues,\
                                  psrStartingObsTimes,\
                                  A,alpha,beta,fref,Ti,c)
plt.plot(T,snr,label='Time shuffle results 1')






"""
snr = np.zeros(len(T))
psrTimesDouble = psrStartingObsTimes.copy()
psrTimesDouble.update((x,y*2) for x,y in psrTimesDouble.items())
for i, Ti in enumerate(TInSeconds):

    snr[i] = snrFunctions.avePTASNR(psrNames,\
                                  psrObsConstants,\
                                  hdValues,\
                                  psrTimesDouble,\
                                  A,alpha,beta,fref,Ti,c)
plt.plot(T,snr,label='2*Tobs')

snr = np.zeros(len(T))
psrTimesHalf = psrStartingObsTimes.copy()
psrTimesHalf.update((x,y*0.5) for x,y in psrTimesHalf.items())
for i, Ti in enumerate(TInSeconds):

    snr[i] = snrFunctions.avePTASNR(psrNames,\
                                  psrObsConstants,\
                                  hdValues,\
                                  psrTimesHalf,\
                                  A,alpha,beta,fref,Ti,c)
plt.plot(T,snr,label='0.5*Tobs')


snr = np.zeros(len(T))

psrTimesShuffle = psrStartingObsTimes.copy()
print(psrTimesShuffle['J1757-5322'])
psrTimesShuffle['J1757-5322']+=psrTimesShuffle['J1730-2304']
psrTimesShuffle['J1730-2304']=0
for i, Ti in enumerate(TInSeconds):

    snr[i] = snrFunctions.avePTASNR(psrNames,\
                                  psrObsConstants,\
                                  hdValues,\
                                  psrTimesShuffle,\
                                  A,alpha,beta,fref,Ti,c)
plt.plot(T,snr,label='shuffle result')

"""
plt.legend()
plt.xlabel('Time (years)')
plt.ylabel('SNR')
plt.savefig('SNRTestDoubleTime.png')
plt.show()
