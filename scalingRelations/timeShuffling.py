import numpy as np
import matplotlib.pyplot as plt
import snrFunctions
import readInData


def plotResult(newTimes,hdValues,psrObsConstants,shuffleNumber,resultsDir):

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
                                  newTimes,\
                                  A,alpha,beta,fref,Ti,c)
    plt.plot(T,snr,label='shuffle {}'.format(shuffleNumber))

    plt.savefig('{}/shuffleSNR_{}.png'.format(resultsDir,shuffleNumber))
    plt.close()
    plt.clf()
    return None
    



psrDataFile = '../data/psrDetails.dat'
#psrDataFile = '../data/trialDataFix.dat'
dataOriginalFormat = np.genfromtxt(psrDataFile, names=True)

psrNames, \
psrObsConstants, \
psrStartingObsTimes, \
angles, \
hdValues = readInData.readDataIntoDicts(psrDataFile)



totalTime=0
for ipsr in psrNames:
    totalTime+=psrStartingObsTimes[ipsr]


oneYearInSeconds = (365.25*24.*60.*60.)

T = 10.
TInSeconds = T * oneYearInSeconds

A = 2.E-15
beta = 13./3
alpha = (3.-beta)/2.
fref = 1./oneYearInSeconds

c = 26./oneYearInSeconds

startSNR = snrFunctions.avePTASNR(psrNames,\
                                  psrObsConstants,\
                                  hdValues,\
                                  psrStartingObsTimes,\
                                  A,alpha,beta,fref,TInSeconds,c)

#print(psrStartingObsTimes.values)
#print(sum(psrStartingObsTimes.values()))
#exit()



#orderedPSRList = sorted(zip(psrObsConstants.values(),psrObsConstants.keys()),reverse=True)

# go through list and allocate time to other pulsars
bestSNRRatioSoFar = 1.
psrDataFileUpdate = psrDataFile

resultsDir = 'halfTimeShuffle'

improvement = True
count = 0

psrTimeShuffle = psrStartingObsTimes.copy()
psrNames = list(psrNames)

logFile=open('{}/shuffleLog.dat'.format(resultsDir),'w')

while improvement==True:

    print(count)    
    check=0
    
    for ipsr in psrNames:
        for jpsr in psrNames: 

            # give ipsr time to jpsr
            editedTimes = psrTimeShuffle.copy()
            timeToShift = editedTimes[ipsr]/2.
            editedTimes[ipsr] -= timeToShift
            editedTimes[jpsr] += timeToShift

            snr = snrFunctions.avePTASNR(psrNames,\
                                         psrObsConstants,\
                                         hdValues,\
                                         editedTimes,\
                                         A,alpha,beta,fref,TInSeconds,c)
            currentRatio = snr/startSNR
            #print(currentRatio)
            if currentRatio > 1 and currentRatio>bestSNRRatioSoFar: 
                bestSNRRatioSoFar = currentRatio
                bestSNR = snr
                psrToGive = ipsr
                psrToTake = jpsr
                #print(improvement)
                check=1
            else:
                pass

    if check==0: 
        print('no improvement found')
        break


    # psrTimeShuffle
    timeToShift = psrTimeShuffle[psrToGive]/2.
    psrTimeShuffle[psrToGive] -= timeToShift
    psrTimeShuffle[psrToTake] += timeToShift

    # remove pulsar
    #print(psrToGive)
    #print(psrNames,len(psrNames))
    #psrNames.remove(psrToGive)

    print(psrNames,len(psrNames))
    logFile.write("""
    NPSRs {}
    Best outcome:
    give {} time to {}
    SNR: {}
    SNR ratio so far: {}
    Total time: {}
    \n\n
    """.format(len(psrNames), psrToGive,psrToTake,bestSNR,bestSNRRatioSoFar,sum(psrTimeShuffle.values())))


    
    # save outcome  
    shuffleResult = open('{}/shuffle_{}.dat'.format(resultsDir,count),'w')
    for i,psr in enumerate(psrNames): 
        shuffleResult.write('{}\t{}\n'.format(psr, psrTimeShuffle[psr]))
    shuffleResult.close()
    plotResult(psrTimeShuffle,hdValues,psrObsConstants,count,resultsDir)
    # plot shuffle result
    
    count+=1
        
#PSR	RA	DEC	IntTime	ExpPrecision
logFile.close()  

"""
Best outcome - only takes one out 
    NPSRs 81
    Best outcome:
    give J2322+2057 time to J1757-5322
    SNR: 15.111226590047005
    SNR radio so far: 1.009914243988591

"""


