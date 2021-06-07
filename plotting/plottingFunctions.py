"""
Functions for plotting results
"""




def plotSNRVTimeCompare(psrNames,psrObsConstants,hdValues,psrTimes,\
                        redAs,redGammas,jitters,label,linestyle='solid'): 

    # general stuff 
    oneYearInSeconds = (365.25*24.*60.*60.)
    A = 2.E-15
    beta = 13./3
    alpha = (3.-beta)/2.
    fref = 1./oneYearInSeconds
    c = 26./oneYearInSeconds

    T = np.linspace(1., 11., 50)
    TInSeconds = T * oneYearInSeconds


    snr = np.zeros(len(T))
    for i, Ti in enumerate(TInSeconds):

        snr[i] = snrFunctions.avePTASNR(psrNames,\
                                        psrObsConstants,\
                                        hdValues,\
                                        psrTimes,\
                                        redAs, redGammas,jitters, \
                                        A,alpha,beta,fref,Ti,c)
    plt.plot(T,snr,label=label,ls=linestyle)
    return None







def newOldTimesSteps(psrNames,psrStartingObsTimes,psrShuffleTimes,step):

    # scatter plot  
    startT, shuffleT = np.zeros(len(psrNames)), np.zeros(len(psrNames))
    for i,psr in enumerate(psrNames):
        startT[i]   = psrStartingObsTimes[psr] 
        shuffleT[i] = psrShuffleTimes[psr]

    
    plt.clf()
    plt.figure(figsize = (4,18))
    plt.scatter(startT, psrNames, label='original')
    plt.scatter(shuffleT,psrNames, label='new {}'.format(int(step)))
    plt.xlabel('tobs (s)')
    plt.tight_layout()
    plt.legend()

    if   int(step<10):  strStep = '00'+str(int(step))
    elif int(step<100): strStep = '0'+str(int(step))
    else: strStep = str(int(step))

    plt.savefig('times-{}.png'.format(strStep))

    return 
    

