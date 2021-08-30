"""
plotting the PTA sensitivity to the GWB if there was 
a 10% improvement for J1909 only and for all PSRs
"""


import numpy as np
import matplotlib.pyplot as plt
import copy
import sys

sys.path.append("../snr/")
import snrFunctions
import readInData



def plotSNRVTime(psrNames,sigmas,angCorrValues,redAs,redGammas,jitters,pltLabel,ls='-'):
    """ 
    plotting the average PTA SNR 
    """

    # constants 
    oneYearInSeconds = (365.25*24.*60.*60.)

    A = 2.E-15                   # Amplitude of background
    beta = 13./3
    alpha = (3.-beta)/2.
    fref = 1./oneYearInSeconds   # reference frequency
    c = 26./oneYearInSeconds     # cadence


    # PTA duration for the plot
    T = np.linspace(1., 10., 20)
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

    return snr[-1]



# data files to use 
originalFile     = '../data/psrDetails.dat'
redNoiseFile    = '../data/redNoise.dat'
jitterNoiseFile = '../data/jitterNoise.dat'

whichCorrelationFunction='HD'

# read in data using 
psrNames, \
psrObsConstants, \
psrStartObsTimes, \
angles, \
hdValues, \
ampRed, \
gammaRed, \
jitterNoise = readInData.readDataIntoDicts(originalFile, \
                                           whichCorrelationFunction, \
                                           redNoiseFile=redNoiseFile, \
                                           jitterNoiseFile=jitterNoiseFile)


# get the precisions (sigmas) needed for plotting
tmpPSRNames = np.genfromtxt(originalFile,usecols=0,dtype=str)
data = np.genfromtxt(originalFile,names=True)
sigmasOriginal = {}
for psr,sigInMilliSec in zip(tmpPSRNames,data['ExpPrecision']):
    sigmasOriginal[psr] = sigInMilliSec * 1E-6

""""
we want to plot four scenarios
  - PTA's orignal timing precision from Ryan's list
  - with J1909-3744 when have 10% timing precision improvement
  - with 30 selected pulsars when have 10% improvement 
  - all 89 pulsars when have 10% improvement
"""

plt.clf()
plt.rcParams.update({'font.size': 12})

"""
original
"""
plotSNRVTime(psrNames,
             sigmasOriginal,
             hdValues,
             ampRed,gammaRed,
             jitterNoise,
             'Original timing precision',
             ls='solid')


"""
J1909 improved precision by 10%
"""
sigJ1909Original = sigmasOriginal['J1909-3744']
sigmasJ1909Improvement = copy.deepcopy(sigmasOriginal)
sigmasJ1909Improvement['J1909-3744'] = sigJ1909Original*0.9


print(sigmasOriginal['J1909-3744'])
print(sigmasJ1909Improvement['J1909-3744'])
plotSNRVTime(psrNames,
             sigmasJ1909Improvement,
             hdValues,
             ampRed,gammaRed,
             jitterNoise,
            'J1909+3744 10% precision improvement',
             ls='dashed')


"""
selected pulsars improvement 
"""
improvers = np.genfromtxt('./improvedPSRsFromMatt.dat',dtype=str)
sigmasSelectedImprovement = copy.deepcopy(sigmasOriginal)
for psrImprover in improvers:
    sigmasSelectedImprovement[psrImprover] = sigmasOriginal[psrImprover]*.9

plotSNRVTime(psrNames,
             sigmasSelectedImprovement,
             hdValues,
             ampRed,gammaRed,
             jitterNoise,
             'Selected 10% precision improvement',
             ls='dotted')


"""
all pulsars 10% improvement
"""
sigmasAllImprovement = copy.deepcopy(sigmasOriginal)
for psr in psrNames: 
    sigmasAllImprovement[psr] = sigmasOriginal[psr]*0.9
plotSNRVTime(psrNames,
             sigmasAllImprovement,
             hdValues,
             ampRed,gammaRed,
             jitterNoise,
             'All 10% precision improvement',
             ls='dashdot')



plt.ylim(0,14)
plt.xlim(1,10)
plt.legend(fontsize=10)
plt.ylabel('Average pulsar timing array S/N')
plt.xlabel('Time (years)')
plt.tight_layout()
plt.savefig('PTASNR.pdf')
plt.show()

