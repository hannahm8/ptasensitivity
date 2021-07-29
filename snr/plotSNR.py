import numpy as np
import matplotlib.pyplot as plt
import snrFunctions
import readInData


# red in th deat
psrDataFile = '../data/psrDetails.dat'

psrNames, \
psrObsConstants, \
psrStartingObsTimes, \
angles, \
hdValues = readInData.readDataIntoDicts(psrDataFile)

totalTime=0
for ipsr in psrNames:
    totalTime+=psrStartingObsTimes[ipsr]


oneYearInSeconds = (365.25*24.*60.*60.)

Ts = np.linspace(1,15,10)
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



# the worst constant 
psrW, psrWName = max(zip(psrObsConstants.values(), psrObsConstants.keys()))
# the best constnt
psrB, psrBName = min(zip(psrObsConstants.values(), psrObsConstants.keys()))




# update
psrStartingObsTimes[psrWName] = 0
snrBestSoFar = 0
betterCount=0
for ipsr in psrNames: 

    print(ipsr)
    if ipsr==psrWName:
        print('pass')
        pass
    else:
        # make a copy 
        editedTimes = psrStartingObsTimes
        editedTimes[ipsr]+=psrW

        snrTrials = np.zeros(len(Ts))
        for i, Ti in enumerate(TInSeconds):
            snrTrials[i] = snrFunctions.avePTASNR(psrNames,\
                                            psrObsConstants,\
                                            hdValues,\
                                            editedTimes,\
                                            A,alpha,beta,fref,Ti,c)

        plt.plot(Ts,snrTrials,color='k',alpha=0.5)
        if snrTrials[-1]>snrBestSoFar:
            snrBestSoFar = snrTrials[-1]
            print('oooo',snrTrials[-1])
            betterCount+=1



plt.plot(Ts,snr,label='v0',color='y')


plt.legend()


plt.grid()
plt.show()

