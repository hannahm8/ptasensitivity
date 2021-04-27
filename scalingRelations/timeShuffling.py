import numpy as np
import matplotlib.pyplot as plt
import snrFunctions
import readInData


def plotResult(newTimes,angCorrValues,psrObsConstants,\
               shuffleNumber,resultsDir,redAs,redGs):


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

        snr[i] = snrFunctions.avePTASNR_incRedNoise(psrNames,\
                                  psrObsConstants,\
                                  angCorrValues,\
                                  newTimes,\
                                  redAs, redGs, \
                                  A,alpha,beta,fref,Ti,c)
    plt.plot(T,snr,label='shuffle {}'.format(shuffleNumber))

    plt.savefig('{}/shuffleSNR_{}.png'.format(resultsDir,shuffleNumber))
    plt.close()
    plt.clf()
    return None
    



psrDataFile = '/fred/oz005/users/hmiddlet/ptasensitivity/data/psrDetails.dat'
#psrDataFile = '/home/hannahm/repositories/ptasensitivity/data/trialPSRDataShort.dat'
dataOriginalFormat = np.genfromtxt(psrDataFile, names=True)



#redNoisePath='/home/hannahm/repositories/ptasensitivity/data/redNoise.dat'
"""
We are using the dp-hd value here!
For shuffling to avoid losing pulsars at larger angular separation
"""
psrNames, \
psrObsConstants, \
psrStartingObsTimes, \
angles, \
angCorrValues, \
redAmps, \
redGammas = readInData.readDataIntoDicts_dphdDiff(psrDataFile,\
                                                  redNoiseFile=redNoisePath)


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

startSNR = snrFunctions.avePTASNR_incRedNoise(psrNames,\
                                              psrObsConstants,\
                                              angCorrValues,\
                                              psrStartingObsTimes,\
                                              redAmps, redGammas, \
                                              A,alpha,beta,fref,TInSeconds,c)




#orderedPSRList = sorted(zip(psrObsConstants.values(),psrObsConstants.keys()),reverse=True)

# go through list and allocate time to other pulsars
bestSNRRatioSoFar = 1.
psrDataFileUpdate = psrDataFile

resultsDir = './times'

improvement = True
count = 0

psrTimeShuffle = psrStartingObsTimes.copy()
psrNames = list(psrNames)

# clear log file
logFile=open('shuffleLog.dat'.format(resultsDir),'w')
logFile.close()

while improvement==True:

    #print(count)    
    check=0
    
    for ipsr in psrNames:
        for jpsr in psrNames: 

            # give ipsr time to jpsr
            editedTimes = psrTimeShuffle.copy()
            timeToShift = editedTimes[ipsr]/4.
            editedTimes[ipsr] -= timeToShift
            editedTimes[jpsr] += timeToShift

            snr = snrFunctions.avePTASNR_incRedNoise(psrNames,\
                                                     psrObsConstants,\
                                                     angCorrValues,\
                                                     editedTimes,\
                                                     redAmps,redGammas, \
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
        logFile=open('shuffleLog.dat','a')
        logFile.write('\n\nno improvement found\n\n')
        logFile.close()	
        break


    # psrTimeShuffle
    timeToShift = psrTimeShuffle[psrToGive]/4.
    psrTimeShuffle[psrToGive] -= timeToShift
    psrTimeShuffle[psrToTake] += timeToShift

    # remove pulsar
    #print(psrToGive)
    #print(psrNames,len(psrNames))
    #psrNames.remove(psrToGive)

    print(psrNames,len(psrNames))
    logFile=open('shuffleLog.dat'.format(resultsDir),'a')
    logFile.write("""
    Shuffle {}
    NPSRs {}
    Best outcome:
    give {} time to {}
    SNR: {}
    SNR ratio so far: {}
    Total time: {}
    \n\n
    """.format(count,len(psrNames), psrToGive,psrToTake,bestSNR,bestSNRRatioSoFar,sum(psrTimeShuffle.values())))
    logFile.close()

    
    # save outcome  
    shuffleResult = open('{}/shuffle_{}.dat'.format(resultsDir,count),'w')
    for i,psr in enumerate(psrNames): 
        shuffleResult.write('{}\t{}\n'.format(psr, psrTimeShuffle[psr]))
    shuffleResult.close()
    
    count+=1


plotResult(psrTimeShuffle,angCorrValues,psrObsConstants,count,resultsDir,redAmps,redGammas)
       

"""
Best outcome - only takes one out 
    NPSRs 81
    Best outcome:
    give J2322+2057 time to J1757-5322
    SNR: 15.111226590047005
    SNR radio so far: 1.009914243988591

"""


