"""
plotting the PTA sensitivity to the GWB if there was 
a 10% improvement for J1909 only and for all PSRs
"""


import numpy as np
import matplotlib.pyplot as plt
import snrFunctions
import readInData


def plotSNRVTime(psrNames,sigmas,angCorrValues,redAs,redGammas,jitters,pltLabel,ls='-'):

    oneYearInSeconds = (365.25*24.*60.*60.)

    A = 2.E-15                   # Amplitude of background
    beta = 13./3
    alpha = (3.-beta)/2.
    fref = 1./oneYearInSeconds   # reference frequency
    c = 26./oneYearInSeconds     # cadence


    # PTA duration for the plot
    T = np.linspace(1., 11., 20)
    TInSeconds = T * oneYearInSeconds


    snr10y = snrFunctions.avePTASNR_from_sigmas(psrNames,sigmas, \
                                                angCorrValues, \
                                                redAs,redGammas,jitters,\
                                                A,alpha,beta,fref,10*oneYearInSeconds,c)
    print('value after 10 years', pltLabel, snr10y)
    

    snr = np.zeros(len(T))
    for i, Ti in enumerate(TInSeconds):

        snr[i] = snrFunctions.avePTASNR_from_sigmas(psrNames,sigmas, \
                                                    angCorrValues, \
                                                    redAs,redGammas,jitters,\
                                                    A,alpha,beta,fref,Ti,c)

        
    plt.plot(T,snr,label=pltLabel,ls=ls)  
    print(snr[-1])
    return None



# data files to use 
orginalFile     = '../data/psrDetails.dat'
redNoiseFile    = '../data/redNoise.dat'
jitterNoiseFile = '../data/jitterNoise.dat'

whichCorrelationFunction='HD'

# read in data using 
psrNames, \
_, \
_, \
angles, \
hdValues, \
ampRed, \
gammaRed, \
jitterNoise = readInData.readDataIntoDicts(orginalFile, \
                                           whichCorrelationFunction, \
                                           redNoiseFile=redNoiseFile, \
                                           jitterNoiseFile=jitterNoiseFile)


# get sigmas from file directly 
tmpPSRNames = np.genfromtxt('../data/psrDetails.dat',usecols=0,dtype=str)
data = np.genfromtxt('../data/psrDetails.dat',names=True)
sigmasOriginal = {}
tmpSigmas = data['ExpPrecision']
for psr,sigInMS in zip(tmpPSRNames,data['ExpPrecision']):
    sigmasOriginal[psr] = sigInMS * 1.E-6



plt.clf()

#### original sigmas from spreadsheet
plotSNRVTime(psrNames,sigmasOriginal,hdValues, \
             ampRed,gammaRed,jitterNoise,'original')



#### update J1909 with 10% improvement
sigmaJ1909 = sigmasOriginal['J1909-3744']
sigmasJ1909Improvement = sigmasOriginal.copy()
sigmasJ1909Improvement['J1909-3744'] = sigmaJ1909*0.9
plotSNRVTime(psrNames,sigmasJ1909Improvement,hdValues, \
             ampRed,gammaRed,jitterNoise,'J1909 10% improvement',ls='--')



#### update improver list with 10% improvement
improvers = np.genfromtxt('../data/improvedPSRsFromMatt.dat',dtype=str)
sigmasSelectedImprovement = sigmasOriginal.copy()
for psrI in improvers:
    print(psrI)
    if psrI=='J1804-2717': pass
    else: sigmasSelectedImprovement[psrI] = sigmasOriginal[psrI]*.9
plotSNRVTime(psrNames,sigmasSelectedImprovement,hdValues, \
             ampRed,gammaRed,jitterNoise,'Selected 10% improvement')

#### update everyone with 10% improvement
sigmasAllImprovement = sigmasOriginal.copy()
for psr in psrNames:
    sigmasAllImprovement[psr] = sigmasOriginal[psr]*0.9
plotSNRVTime(psrNames,sigmasAllImprovement,hdValues, \
             ampRed,gammaRed,jitterNoise,'All 10% improvement')

plt.legend()
plt.ylabel('PTA SNR')
plt.xlabel('Time (years)')
plt.savefig('snr-improvement.png')
plt.show()


