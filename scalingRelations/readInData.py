import numpy as np
import snrFunctions


def readDataIntoDicts(psrDataFile):

    psrNames = np.genfromtxt(psrDataFile,usecols=0,dtype=str)
    psrData = np.genfromtxt(psrDataFile,names=True)

    # compute obs constants 
    obsConstants = [ sig * np.sqrt(intT) for sig, intT in \
                                         zip (psrData['ExpPrecision']*1.E-6, \
                                              psrData['IntTime'])]
    psrObsConstants = {}
    psrStartingObsTimes = {}
    for psr, oc, st in zip(psrNames, obsConstants, psrData['IntTime']):
        psrObsConstants[psr] = float(oc)
        psrStartingObsTimes[psr] = float(st)


    # work out angles and hd values ahead of time
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
                elati, elongi = psrData['ELAT'][i], psrData['ELONG'][i]
                elatj, elongj = psrData['ELAT'][j], psrData['ELONG'][j]
                angle = snrFunctions.h2(elati,elatj,elongi,elongj)
                hd = snrFunctions.hellings_downs(angle)
            onePSRAng[jpsr] = angle
            onePSRHDs[jpsr] = hd
        angles[ipsr]   = onePSRAng
        hdValues[ipsr] = onePSRHDs

    return psrNames, psrObsConstants, psrStartingObsTimes, angles, hdValues
