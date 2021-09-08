import numpy as np
import sys

sys.path.append('/home/hannahm/repositories/ptasensitivity/snr/')
sys.path.append('/fred/oz005/users/hmiddlet/ptasensitivity/snr/')
import snrFunctions
import readInData


def lnlike(theta, psrNames, psrConstants, angCorrelationValues, rAs, gammas, jitters): 



    oneYearInSeconds = (365.25*24.*60.*60.)

    T = 10.
    TInSeconds = T * oneYearInSeconds

    A = 2.E-15
    beta = 13./3
    alpha = (3.-beta)/2.
    fref = 1./oneYearInSeconds

    c = 26./oneYearInSeconds

    obsTimes = {}
    for psr,time in zip(psrNames,theta):
        obsTimes[psr] = time 
    

    snr = snrFunctions.avePTASNR(psrNames,psrConstants,angCorrelationValues, \
                                 obsTimes,rAs,gammas,jitters,\
                                 A,alpha,beta,fref,TInSeconds,c)
    print(snr)
    return snr



def lnprior(theta):
    minTime = 256
    maxTime = 2560
    theta = np.atleast_1d(theta)
    maxTotal = 44096.

    for t in theta:

        if minTime < float(t) < maxTime and sum(theta) <= maxTotal:

            return 0.0

    return -np.inf




def lnprob(theta, psrNames, psrObsConstants, angCorrValues, rAs, gammas, jitters):
    ll = lnlike(theta,psrNames, psrObsConstants, angCorrValues, rAs, gammas, jitters)


    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + ll


psrDataFile = '../data/psrDetails.dat'
redNoisePath = '../data/redNoise.dat'
jitterNoisePath = '../data/jitterNoise.dat'
chooseCorrFunc = 'HD'

# read in the data
psrNames, \
psrObsConstants, \
psrStartingObsTimes, \
angles, \
angCorrValues, \
redAmps, \
redGammas, \
jitters = readInData.readDataIntoDicts(psrDataFile,\
                                       chooseCorrFunc,\
                                       redNoiseFile=redNoisePath,\
                                       jitterNoiseFile=jitterNoisePath)






startingPositions = [ psrStartingObsTimes[psr] for psr in psrNames]

print(lnprob(startingPositions, psrNames, psrObsConstants, angCorrValues, redAmps, redGammas, jitters))



nPSRs = len(psrNames)
ndim, nwalkers = nPSRs, 100

pos = [ startingPositions + np.random.rand(ndim) for i in range(nwalkers) ]

import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(psrNames, psrObsConstants, angCorrValues, redAmps, redGammas, jitters))

sampler.run_mcmc(pos,50)

samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

#saving the samples 
savehere = open('savesamples.dat', 'w')
for i in (samples):
    savehere.write(i)


