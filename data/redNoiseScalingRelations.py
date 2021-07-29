import numpy as np


def gammaToGammaSN():

    """
    convert gamma to gammaR - see p17 M. Lam+2017 and also 
    Arzoumanian et al. 2015b    
    This value is a fit fom the population of pulsars in the    
    paper.
    """ 
    gamma = 2.3
    gammaSN = 2.*gamma + 1.

    return gammaSN



def sigmaSNToAr(sigmaSN): 
    """ 
    from equation 15 in M. Lan+2017 
    """
    gammaSN = gammaToGammaSN()
    Tyr = 10.

    ANeedToSortScale = (sigmaSN / 3.E-9) \
                       * (gammaSN - 1.)**0.5 \
                       * (Tyr**((gammaSN-1.)/2.))

    A = ANeedToSortScale*3E-3 # in mu s yr^{1/2}
    A = A * 1E-6 # in s yr^(1/2)
    A = A /(365.25*60*60*24) #in yr^{3/2}
    
    return np.log10(A), gammaSN
    #return A, gammaSN


def convertF0F1ToSigmaSN(f,fDotm15,Tyr):

    """
    functoin takes
    f frequency of pulsar
    fdotm15 frequency derivative of pulsar in units of 10^{-15} s^{-2}
    Tyr Timespan of the dataset in years
    """ 
    
    modFdot = abs(fDotm15)

    # vales from NANO+MSP10PPTA fit in Lam_17 ua(table 4)
    #C2    = np.exp(-1.2)
    #alpha = -0.9
    #beta  = 0.8
    #gamma = 2.3
    C2 = np.exp(-16.)
    alpha = 0.8
    beta = 0.9
    gamma = 4.9

    sigmaSN = C2 * (f**alpha) * (modFdot**beta) * (Tyr**gamma)

    return sigmaSN *1E-6




"""
mu = 339.3156872882446
mudot = -1.614817E-15/ 1E-15
Tyr = 10.

sigmaSN = redNoiseScalingRelation(mu,mudot,Tyr)
print(sigmaRNToAr(sigmaSN))
"""
