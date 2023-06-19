import numpy as np
import matplotlib.pyplot as plt
import argparse
import snrFunctions
import readInData

import datetime
import copy




def getSNR(givePSR,takePSR,modifiedIntegrals,obsConsts,times,\
           redAmps,redGammas,dmAmps,dmGammas,jitters,constants):

    A, alpha, beta, TInSeconds, c, fref = constants


    total=0
    countCalls = 0
    for i,ipsr in enumerate(psrNames):
        for j,jpsr in enumerate(psrNames):
            if (i>j): # no double counting
                corr = angCorrValues[ipsr][jpsr]    


                if ipsr==givePSR or ipsr==takePSR or jpsr==givePSR or jpsr==takePSR:
                    #print('\t',ipsr,jpsr)

                    sigI = obsConsts[ipsr] / np.sqrt(times[ipsr])
                    sigJ = obsConsts[jpsr] / np.sqrt(times[jpsr])

                    # integral including red, dm, jitter nosie

                    intValue = snrFunctions.get_integral_rnoise_dmnoise_jitter(c,fref,\
                                                                       sigI,\
                                                                       redAmps[ipsr],\
                                                                       redGammas[ipsr],\
                                                                       dmAmps[ipsr],\
                                                                       dmGammas[ipsr],\
                                                                       jitters[ipsr],\
                                                                       sigJ,\
                                                                       redAmps[jpsr],\
                                                                       redGammas[jpsr],\
                                                                       dmAmps[jpsr],\
                                                                       dmGammas[jpsr],\
                                                                       jitters[jpsr],\
                                                                       A,alpha,beta,TInSeconds)
                    '''                                                   
                    intValue = snrFunctions.get_integral_rnoise_jitter(c,fref,\
                                                                       sigI,\
                                                                       redAmps[ipsr],\
                                                                       redGammas[ipsr],\
                                                                       jitters[ipsr],\
                                                                       sigJ,\
                                                                       redAmps[jpsr],\
                                                                       redGammas[jpsr],\
                                                                       jitters[jpsr],\
                                                                       A,alpha,beta,TInSeconds)
                    '''
                    # update the precompute values

                    #print(intValue, modifiedIntegrals[ipsr][jpsr])
                    modifiedIntegrals[ipsr][jpsr] = intValue
                    #print(ipsr,jpsr,modifiedIntegrals[ipsr][jpsr])
                    countCalls+=1
                else: pass

                integral = modifiedIntegrals[ipsr][jpsr]

                aveSNRSinglePulsarPair = 2*TInSeconds*corr*corr*integral
                total+=aveSNRSinglePulsarPair

    snr = np.sqrt(total)
    #print('count calls: ', countCalls)
    return snr, modifiedIntegrals




def get_arguments():

    
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--psrDataFile', '-d', dest='psrDataFile', 
                        required=True, type=str,
                        help='psr data path (PSR RA DEC IntTime ExpPrecision)')
    parser.add_argument('--redNoiseFile', '-r', dest='psrRedNoiseFile', default=None,
                        required=False, type=str,
                        help='red noise data path (PSR ASN gammaSN Ref?)')
    parser.add_argument('--dmNoiseFile', '-dm', dest='psrDMNoiseFile', default=None,
                        required=False, type=str,
                        help='DM noise data path (PSR ADM gammaDM Ref?)')
    parser.add_argument('--jitterNoiseFile', '-j', dest='psrJitterNoiseFile',
                        required=False, type=str, default=None,
                        help='jitter noise data path (PSR jitter error)')
    parser.add_argument('--whichCorrelation', '-c', dest='whichCorrFunc',
                        required=True, type=str,
                        help='choose correlation function (HD, DPHD, EQUAL)')
    args = parser.parse_args()

    return args
    
start = datetime.datetime.now()



args = get_arguments()
psrDataFile     = args.psrDataFile
chooseCorrFunc  = args.whichCorrFunc
redNoisePath    = args.psrRedNoiseFile
dmNoisePath     = args.psrDMNoiseFile
jitterNoisePath = args.psrJitterNoiseFile

print("""
data {}
corr {}
redn {}
jitn {}
""".format(psrDataFile, chooseCorrFunc, redNoisePath, dmNoisePath, jitterNoisePath))




# read in the data
psrNames, \
psrObsConstants, \
psrStartingObsTimes, \
angles, \
angCorrValues, \
redAmps, \
redGammas, \
dmAmps, \
dmGammas, \
jitters = readInData.readDataIntoDicts(psrDataFile,\
                                       chooseCorrFunc,\
                                       redNoiseFile=redNoisePath,\
                                       dmNoiseFile=dmNoisePath,\
                                       jitterNoiseFile=jitterNoisePath)

oneYearInSeconds = (365.25*24.*60.*60.)

T = 10.
TInSeconds = T * oneYearInSeconds

A = 2.E-15
beta = 13./3
alpha = (3.-beta)/2.
fref = 1./oneYearInSeconds

c = 26./oneYearInSeconds

constants = A, alpha, beta, TInSeconds, c, fref

# precompute the integrals for each pulsar pair. This is done once
intContributions = {}
for ipsr in psrNames:
    thisInt = {}
    for jpsr in psrNames: 

        sigI = psrObsConstants[ipsr] / np.sqrt(psrStartingObsTimes[ipsr])
        sigJ = psrObsConstants[jpsr] / np.sqrt(psrStartingObsTimes[jpsr])

        thisInt[jpsr] = snrFunctions.get_integral_rnoise_dmnoise_jitter(c,fref,\
                                                                sigI,\
                                                                redAmps[ipsr],
                                                                redGammas[ipsr],\
                                                                dmAmps[ipsr],\
                                                                dmGammas[ipsr],\
                                                                jitters[ipsr],\
                                                                sigJ,\
                                                                redAmps[jpsr],\
                                                                redGammas[jpsr],\
                                                                dmAmps[jpsr],\
                                                                dmGammas[jpsr],\
                                                                jitters[jpsr],\
                                                                A,alpha,beta,TInSeconds)

    intContributions[ipsr] = thisInt



startSNR = snrFunctions.avePTASNR(psrNames,\
                                  psrObsConstants,\
                                  angCorrValues,\
                                  psrStartingObsTimes,\
                                  redAmps, redGammas, \
                                  dmAmps, dmGammas, \
                                  jitters, \
                                  A,alpha,beta,fref,TInSeconds,c)



print('startSNR', startSNR)


total=0
for i,ipsr in enumerate(psrNames):
    for j,jpsr in enumerate(psrNames):
        if (i>j): # no double counting
            corr = angCorrValues[ipsr][jpsr]    
            inte = intContributions[ipsr][jpsr]
            #print(TInSeconds,corr,inte,type(TInSeconds),type(corr),type(inte))
            aveSNRSinglePulsarPair = 2*TInSeconds*corr*corr*inte
            total+=aveSNRSinglePulsarPair

snr = np.sqrt(total)
print('snr', snr)

# what's the nones about? XX does this do anything at the moment?
snr,mod = getSNR('None','None',intContributions,psrObsConstants,psrStartingObsTimes,\
                  redAmps,redGammas,dmAmps,dmGammas,jitters,constants)
print(snr)

############################################################

bestSNRRatioSoFar = 1.
psrDataFileUpdate = psrDataFile

resultsDir = './times'

improvement = True
count = 0

psrTimeShuffle = psrStartingObsTimes.copy()
psrNames = list(psrNames)

intContribShuffle = intContributions.copy()

# clear log file
logFile=open('shuffleLog.dat'.format(resultsDir),'w')
logFile.close()

shuffleHowMuch=8. # eighth

minimumTime = 256.
minTimeCondition = minimumTime + (minimumTime/(shuffleHowMuch-1.))

while improvement==True: 

    print('shuffle: ', count)
    check=0

    for ipsr in psrNames: 
        for jpsr in psrNames: 
            # this will make the minimum time 256 seconds (342 - 341/4)
            if (ipsr!=jpsr) and psrTimeShuffle[ipsr]>minTimeCondition: 
                #print(ipsr,jpsr)
                # updating the times for the trial run 

                editedTimes = psrTimeShuffle.copy()
                timeToShift = editedTimes[ipsr]/shuffleHowMuch
                editedTimes[ipsr] -= timeToShift
                editedTimes[jpsr] += timeToShift
                

                editedIntContrib = copy.deepcopy(intContribShuffle) 
                
                # get the snr using the partial precompute method. 
                snr,modInt = getSNR(ipsr,jpsr,editedIntContrib,\
                                    psrObsConstants,editedTimes,\
                                    redAmps,redGammas,\
                                    dmAmps,dmGammas,\
                                    jitters,constants)  

                
                # check if better
                currentRatio=snr/startSNR
                if currentRatio>1 and currentRatio>bestSNRRatioSoFar: 
                    bestSNRRatioSoFar = currentRatio
                    bestSNR = snr
                    psrToGive = ipsr    
                    psrToTake = jpsr
                    intValuesSave = copy.deepcopy(modInt)
                    check=1 # an improvement was made with this suffle round
                else:   
                    pass

    if check==0:
        logFile=open('shuffleLog.dat','a') 
        logFile.write('\n\nno improvement found\n\n')
        logFile.close()
        break # there is not more to be gained with this method

    # if improvement has been made, update the times and the intergrals
    timeToShift = psrTimeShuffle[psrToGive]/shuffleHowMuch
    psrTimeShuffle[psrToGive] -= timeToShift    
    psrTimeShuffle[psrToTake] += timeToShift
    
    # update the integrals too 
    for ipsr in psrNames:   
        for jpsr in psrNames: 
            intContribShuffle[ipsr][jpsr] = intValuesSave[ipsr][jpsr]
            
    print(bestSNR,psrToGive,psrToTake)
    #print(psrNames,len(psrNames))
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




end = datetime.datetime.now()


print('time: ', end-start)




exit()


###################################################################

start = datetime.datetime.now()
# check the total average snr
total=0
for i,ipsr in enumerate(psrNames):
    for j,jpsr in enumerate(psrNames):
        if (i>j): # no double counting
            corr = angCorrValues[ipsr][jpsr]    
            inte = intContributions[ipsr][jpsr]
            aveSNRSinglePulsarPair = 2*TInSeconds*corr*corr*inte
            total+=aveSNRSinglePulsarPair

snr = np.sqrt(total)

end = datetime.datetime.now()
preComputeTime=end-start
print(snr, preComputeTime)


start = datetime.datetime.now()
startSNR = snrFunctions.avePTASNR(psrNames,\
                                  psrObsConstants,\
                                  angCorrValues,\
                                  psrStartingObsTimes,\
                                  redAmps, redGammas, \
                                  jitters, \
                                  A,alpha,beta,fref,TInSeconds,c)
end = datetime.datetime.now()
oldWay = end-start
print(startSNR, oldWay)

print('\n percentage speed up ', 100*(oldWay - preComputeTime)/oldWay)




