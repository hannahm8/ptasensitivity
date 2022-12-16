import numpy as np
import snrFunctions


def readRedNoise(fileName,allPSRNames):
    
    redPSRNames = np.genfromtxt(fileName,usecols=0,dtype=str)
    redData = np.genfromtxt(fileName,names=True)
    
    As     = {}
    gammas = {}
    for rpsr, a, g in zip(redPSRNames,redData['ASN'],redData['gammaSN']):
        As[rpsr]     = a
        gammas[rpsr] = g
    
    ARed = {}
    gRed = {}
    for psr in allPSRNames: 
        try: 
            ARed[psr] = 10.**As[psr] 
            gRed[psr] = gammas[psr]
        except:
            ARed[psr] = 0.
            gRed[psr] = 1.

    return ARed, gRed 


def readDMNoise(fileName,allPSRNames):

    dmPSRNames = np.genfromtxt(fileName,usecols=0,dtype=str)
    dmData = np.genfromtxt(fileName,names=True)
    
    As     = {}
    gammas = {}
    for dmpsr, a, g in zip(dmPSRNames,dmData['ADM'],dmData['gammaDM']):
        As[dmpsr]     = a
        gammas[dmpsr] = g
        
    ADM = {}
    gDM = {}
    for psr in allPSRNames:
        try: 
            ADM[psr] = 10.**As[psr]
            gDM[psr] = gammas[psr]
        except: 
            ADM[psr] = 0.
            gDM[psr] = 1.
            
    return ADM, gDM




def readJitterNoise(fileName,allPSRNames):

    jitterPSRNames = np.genfromtxt(fileName,usecols=0,dtype=str)
    jitterData = np.genfromtxt(fileName,names=True)

    jitters = {}
    for jpsr, j in zip(jitterPSRNames, jitterData['jitter']):
        jitters[jpsr] = j

    jit = {}
    for psr in allPSRNames:
        try: 
            jit[psr] = jitters[psr]*1.E-9
        except:     
            jit[psr] = 0.
    
    return jit






def HDCorrs(psrNames,psrData):

    """
    Computes the HD correlation for each pulsar pair. 
    """
    angles = {}
    hdValues = {}
    for i, ipsr in enumerate(psrNames):
        onePSRAng = {}
        onePSRHDs = {}
        for j, jpsr in enumerate(psrNames):
            if jpsr==ipsr:
                angle = 0
                hd = None
            else: 
                rai, deci = psrData['RA'][i], psrData['DEC'][i]
                raj, decj = psrData['RA'][j], psrData['DEC'][j]
                angle = snrFunctions.h2(rai,raj,deci,decj)
                hd = snrFunctions.hellings_downs(angle)
            onePSRAng[jpsr] = angle
            onePSRHDs[jpsr] = hd
        angles[ipsr]   = onePSRAng
        hdValues[ipsr] = onePSRHDs

    return angles, hdValues




def DPHDDiffCorrs(psrNames,psrData):

    """
    Computes the difference between the HD and dipole correlation for each
    pulsar pair. 
    Reduces preference for the galactic centre pulsars.
    """
    angles = {}
    dphdDiffValues = {}
    for i, ipsr in enumerate(psrNames):
        onePSRAng = {}
        onePSRDiffs = {}
        for j, jpsr in enumerate(psrNames):
            if jpsr==ipsr:
                angle = 0
                diff = None
            else: 
                rai, deci = psrData['RA'][i], psrData['DEC'][i]
                raj, decj = psrData['RA'][j], psrData['DEC'][j]
                angle = snrFunctions.h2(rai,raj,deci,decj)
                hd = snrFunctions.hellings_downs(angle)
                dp = snrFunctions.dipole(angle)
                diff = dp-(2.*hd) # normalise both to 1 
            onePSRAng[jpsr] = angle
            onePSRDiffs[jpsr] = diff
        angles[ipsr] = onePSRAng
        dphdDiffValues[ipsr] = onePSRDiffs

    return angles, dphdDiffValues




def equalCorrs(psrNames,psrData):
    """ 
    This returns equal correlation values for all pulsar pairs.
    The idea is to remove the sky posistion bias from the shuffle.
    """
    
    value = 1./float(len(psrNames))

    angles = {}
    equalCorrValues = {}
    for i, ipsr in enumerate(psrNames):
        onePSRAng = {}
        onePSREqual = {}
        for j, jpsr in enumerate(psrNames):
            if jpsr==ipsr:
                angle = 0
                equal = None
            else: 
                rai, deci = psrData['RA'][i], psrData['DEC'][i]
                raj, decj = psrData['RA'][j], psrData['DEC'][j]
                angle = snrFunctions.h2(rai,raj,deci,decj)
                equal = value
            onePSRAng[jpsr] = angle
            onePSREqual[jpsr] = equal
        angles[ipsr] = onePSRAng
        equalCorrValues[ipsr] = onePSREqual

    return angles, equalCorrValues




def readDataIntoDicts(psrDataFile,\
                      whichCorrelationFunction,\
                      redNoiseFile=None,\
                      dmNoiseFile=None,\
                      jitterNoiseFile=None):


    # get the pulsar names and the data 
    psrNames = np.genfromtxt(psrDataFile,usecols=0,dtype=str)
    psrData = np.genfromtxt(psrDataFile,names=True)


    """
    compute constants for each pulsar from the observation times and precision.
    """
    obsConstants = [ sig * np.sqrt(intT) for sig, intT in \
                                         zip (psrData['ExpPrecision']*1.E-6, \
                                              psrData['IntTime'])]
    psrObsConstants = {}
    psrStartingObsTimes = {}
    for psr, oc, st in zip(psrNames, obsConstants, psrData['IntTime']):
        psrObsConstants[psr] = float(oc)
        psrStartingObsTimes[psr] = float(st)



    """
    work out angles and correlation values -> only needs to be computed once 
    before the shuffle
    """
    if whichCorrelationFunction=='HD':
        angles, correlationValues = HDCorrs(psrNames,psrData)
    elif whichCorrelationFunction=='DPHDDiff':
        angles, correlationValues = DPHDDiffCorrs(psrNames,psrData)
    elif whichCorrelationFunction=='EQUAL':
        angles, correlationValues = equalCorrs(psrNames,psrData) # does not exist yet
    else: 
        print('Error: you have not chosen an available correlation option')
        exit()


    """
    read in the red noise parameters if available. Returns: 
        ampRed = 0 
        gammaRed = 1 
    if not available
    """
    if redNoiseFile!=None:
        ampRed, gammaRed = readRedNoise(redNoiseFile,psrNames)
    else: 
        ampRed,gammaRed = {}, {}
        for psr in psrNames:
            ampRed[psr] = 0
            gammaRed[psr] = 1


    """
    read in the dm noies parameter if available. Returns:
        ampDM = 0
        gammaDM = 1
    if not available.
    """    
    if dmNoiseFile!=None:
        ampDM, gammaDM = readDMNoise(dmNoiseFile,psrNames)
    else:
        ampDM, gammaDM = {}, {}
        for psr in psrNames: 
            ampDM[psr] = 0
            gammaDM[psr] = 1


    """
    read in the jitter noise values if available. Returns
        jitter = 0
    for all if not available
    """
    if jitterNoiseFile!=None: 
        jitterNoise = readJitterNoise(jitterNoiseFile,psrNames)
    else: 
        jitterNoise = {}
        for psr in psrNames:    
            jitterNoise[psr] = 0



    return psrNames, \
           psrObsConstants, \
           psrStartingObsTimes, \
           angles, \
           correlationValues, \
           ampRed, gammaRed, \
           ampDM, gammaDM, \
           jitterNoise 


