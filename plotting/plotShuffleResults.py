import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.append('/home/hannahm/repositories/ptasensitivity/snr/')
sys.path.append('/fred/oz005/users/hmiddlet/ptasensitivity/snr/')
sys.path.append('/home/ADF/middlehr/repositories/ptasensitivity/snr/')
import snrFunctions
import readInData

def plotSNRVTimeCompare(psrNames,psrObsConstants,hdValues,psrTimes,\
                        redAs,redGammas,dmAs,dmGammas, \
                        jitters,label,linestyle='solid'): 

    # general stuff 
    oneYearInSeconds = (365.25*24.*60.*60.)
    A = 2.E-15
    beta = 13./3
    alpha = (3.-beta)/2.
    fref = 1./oneYearInSeconds
    c = 26./oneYearInSeconds

    T = np.linspace(1., 10., 50)
    TInSeconds = T * oneYearInSeconds


    snr = np.zeros(len(T))
    for i, Ti in enumerate(TInSeconds):

        snr[i] = snrFunctions.avePTASNR(psrNames,\
                                        psrObsConstants,\
                                        hdValues,\
                                        psrTimes,\
                                        redAs, redGammas,\
                                        dmAs, dmGammas, \
                                        jitters, \
                                        A,alpha,beta,fref,Ti,c)

    T10Yr = 10.*oneYearInSeconds
    snr10Year = snrFunctions.avePTASNR(psrNames,\
                                       psrObsConstants,\
                                       hdValues,\
                                       psrTimes,\
                                       redAs,redGammas,\
                                       dmAs, dmGammas,\
                                       jitters,\
                                       A,alpha,beta,fref,T10Yr,c)
    print("\nSNR at 10 years: {} \n".format(snr10Year))
    plt.plot(T,snr,label=label,ls=linestyle)
    return None


def obsTime(psrNames,psrObsTimes,label,marker='.'):
    
    # move from dict to list
    obsTimes = np.zeros(len(psrNames))
    for i,psr in enumerate(psrNames):
        obsTimes[i] = psrObsTimes[psr]
    print(len(obsTimes),len(psrNames))
    plt.scatter(obsTimes,psrNames,label=label,marker=marker)
    return None
  

def readShuffle(shuffleFile):

    """
    reads the shuffle times into a dictionary
    """
    psrTimeShuffleDataNames = np.genfromtxt(shuffleFile,usecols=0,dtype=str)
    psrTimeShuffleDataTimes = np.genfromtxt(shuffleFile,usecols=1)

    psrShuffleTimes = {}
    for psrName, psrTime in zip(psrTimeShuffleDataNames, psrTimeShuffleDataTimes):
        psrShuffleTimes[psrName] = psrTime
    return psrShuffleTimes


    

def timeComparison(psrNames,psrStartingObsTimes,psrShuffleTimes):

    # scatter plot  
    startT, shuffleT = np.zeros(len(psrNames)), np.zeros(len(psrNames))
    for i,psr in enumerate(psrNames):
        startT[i]   = psrStartingObsTimes[psr] 
        shuffleT[i] = psrShuffleTimes[psr]

    maximumTime = max(max(startT),max(shuffleT)) 
    timeLim = maximumTime *1.1
    plt.scatter(startT,shuffleT)    
    plt.xlim(0,timeLim)
    plt.ylim(0,timeLim)
    plt.xlabel('original obs. time (s)')
    plt.ylabel('new obs. time (s)') 
    plt.savefig('{}/newTVoldT.png'.format(outputDir))
    plt.show()

    plt.clf()

    fracChange = (shuffleT - startT)/startT
    plt.figure(figsize = (4,18))
    plt.scatter(fracChange, psrNames)   
    plt.axvline(0)
    plt.tight_layout()
    plt.legend()
    plt.xlabel('(new tobs - old tobs) / old tobs')
    plt.savefig('{}/timeDiffPerPSR.png'.format(outputDir))
    plt.show()


    plt.clf()
    plt.figure(figsize = (4,18))
    plt.scatter(startT, psrNames, label='original')
    plt.scatter(shuffleT,psrNames, label='new')
    plt.xlabel('tobs (s)')
    plt.tight_layout()
    plt.legend()
    plt.savefig('{}/obsTimes.png'.format(outputDir))
    plt.show()
    
    return startT,shuffleT


# run this like this:
#python ../plotShuffleResults.py ../../data/psrDetails.dat original ./shuffle_3.dat shuffle /path/to/rednoise/file.dat tmp  

# command line arguments
originalFile  = sys.argv[1]
originalLabel = sys.argv[2]

shuffleFile1   = sys.argv[3]
shuffleLabel1  = sys.argv[4]

shuffleFile2   = sys.argv[5]
shuffleLabel2  = sys.argv[6]

shuffleFile3   = sys.argv[7]
shuffleLabel3  = sys.argv[8]

redNoiseFile  = sys.argv[9]

dmNoiseFile = sys.argv[10]

jitterNoiseFile = sys.argv[11]

outputDir = sys.argv[12]


print(shuffleFile2,shuffleFile3)


# data file
#psrDataFile = '../data/psrDetails.dat'

# original data 
# read in data and compute angles etc

whichCorrelationFunction='HD'

psrNames, \
psrObsConstants, \
psrStartingObsTimes, \
angles, \
hdValues, \
ampRed, \
gammaRed, \
ampDM, \
gammaDM, \
jitterNoise = readInData.readDataIntoDicts(originalFile, \
                                        whichCorrelationFunction, \
                                        redNoiseFile=redNoiseFile, \
                                        dmNoiseFile=dmNoiseFile, \
                                        jitterNoiseFile=jitterNoiseFile)


 
# shuffled times 
#psrTimeShuffleDataNames = np.genfromtxt('oneToOneShuffle/shuffle_14.dat',usecols=0,dtype=str)
#psrTimeShuffleDataTimes = np.genfromtxt('oneToOneShuffle/shuffle_14.dat',usecols=1)

# version with multiple files 

if shuffleFile2!='None' and shuffleFile3!='None':
    """
    comparison plot version
    """
    
    plt.rcParams.update({'font.size': 12})
    # read in the data
    shuffleTimes1 = readShuffle(shuffleFile1)
    shuffleTimes2 = readShuffle(shuffleFile2)
    shuffleTimes3 = readShuffle(shuffleFile3)

    # plot the original data
    # set up for plotting 
    times = [psrStartingObsTimes,shuffleTimes1,shuffleTimes2,shuffleTimes3]
    label = [originalLabel,shuffleLabel1,shuffleLabel2,shuffleLabel3]
    linestyles = ['solid','dotted','dashed','dashdot']
    markers = ['o','x','+','*']
     
    for psrTime,runLabel,ls in zip(times,label,linestyles):

        plotSNRVTimeCompare(psrNames, \
                            psrObsConstants, \
                            hdValues, \
                            psrTime, \
                            ampRed, \
                            gammaRed, \
                            ampDM, \
                            gammaDM, \
                            jitterNoise, \
                            runLabel, \
                            linestyle=ls)

    plt.legend()
    plt.xlabel('Time (years)')
    plt.ylabel('Average signal-to-noise ratio')
    plt.savefig('{}/compareSNRVTime.png'.format(outputDir))
    plt.savefig('{}/compareSNRVTime.pdf'.format(outputDir))
    plt.show()
    

    psrNamesFirstHalf  = psrNames[:45]
    psrNamesSecondHalf = psrNames[45:]
    for i,psrNamesHalf in enumerate([psrNamesFirstHalf,psrNamesSecondHalf]):
      plt.clf()
      plt.figure(figsize=(5,10))
      for psrTime,runLabel,mark in zip(times,label,markers):
        obsTime(psrNamesHalf,psrTime,runLabel,marker=mark)
      plt.legend()
      plt.xlabel('Integration time (s)')
      plt.tight_layout()
      plt.savefig('{}/timeComparison-{}.pdf'.format(outputDir,i))
      plt.savefig('{}/timeComparison-{}.png'.format(outputDir,i))
      plt.show()
    exit()

else: 
    shuffleFile  = shuffleFile1
    shuffleLabel = shuffleLabel1



print('plotting snr')
psrTimeShuffleDataNames = np.genfromtxt(shuffleFile,usecols=0,dtype=str)
psrTimeShuffleDataTimes = np.genfromtxt(shuffleFile,usecols=1)

# create dictionary
psrShuffleTimes = {}
for psrName, psrTime in zip(psrTimeShuffleDataNames, psrTimeShuffleDataTimes):
    psrShuffleTimes[psrName] = psrTime



# plotting the SNR vs Time
plotSNRVTimeCompare(psrNames, \
                    psrObsConstants, \
                    hdValues, \
                    psrStartingObsTimes, \
                    ampRed, \
                    gammaRed, \
                    jitterNoise, \
                    originalLabel)

plotSNRVTimeCompare(psrNames, \
                    psrObsConstants, \
                    hdValues, \
                    psrShuffleTimes, \
                    ampRed, \
                    gammaRed, \
                    jitterNoise, \
                    shuffleLabel)


plt.xlabel('Time (years)')
plt.ylabel('SNR')
plt.legend()
plt.savefig('{}/SNRVTime.png'.format(outputDir))
plt.show()



# compare the old and new times for each pulsar & how each changed
startTime, shuffleTime = timeComparison(psrNames,psrStartingObsTimes,psrShuffleTimes)
obsConstants = [ psrObsConstants[ipsr] for ipsr in psrNames]

plt.clf()
plt.scatter(obsConstants,(shuffleTime-startTime)/startTime)
plt.xlabel('tobs-sigma constant')
plt.ylabel('(new time - old time) / old time')
plt.savefig('{}/fractionalTimeDifferenceVObsConstant.png'.format(outputDir))
plt.show()

