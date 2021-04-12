import numpy as np
import sys
import matplotlib.pyplot as plt

import snrFunctions
import readInData

# command line arguments
originalFile  = sys.argv[1]
originalLabel = sys.argv[2]

shuffleFile   = sys.argv[3]
shuffleLabel  = sys.argv[4]

# data file
#psrDataFile = '../data/psrDetails.dat'

# original data 
# read in data and compute angles etc
psrNames, \
psrObsConstants, \
psrStartingObsTimes, \
angles, \
hdValues = readInData.readDataIntoDicts(originalFile)

# shuffled times 
#psrTimeShuffleDataNames = np.genfromtxt('oneToOneShuffle/shuffle_14.dat',usecols=0,dtype=str)
#psrTimeShuffleDataTimes = np.genfromtxt('oneToOneShuffle/shuffle_14.dat',usecols=1)
psrTimeShuffleDataNames = np.genfromtxt(shuffleFile,usecols=0,dtype=str)
psrTimeShuffleDataTimes = np.genfromtxt(shuffleFile,usecols=1)


# create dictionary
psrShuffleTimes = {}
for psrName, psrTime in zip(psrTimeShuffleDataNames, psrTimeShuffleDataTimes):
    psrShuffleTimes[psrName] = psrTime



# plot

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
plt.plot(T,snr,label=originalLabel)


snr = np.zeros(len(T))
for i, Ti in enumerate(TInSeconds):

    snr[i] = snrFunctions.avePTASNR(psrNames,\
                                  psrObsConstants,\
                                  hdValues,\
                                  psrShuffleTimes,\
                                  A,alpha,beta,fref,Ti,c)

plt.plot(T,snr,label=shuffleLabel)



plt.xlabel('Time (years)')
plt.ylabel('SNR')
plt.legend()

plt.show()
