import numpy as np
import scipy.special as sc
from scipy.integrate import quad

"""
 Class. Quantum Grav. 30 224015
"""


def h2(lat1d,lat2d,lon1d,lon2d):

    """
    To get the angle between two pulsars from their ELAT and ELONG in deg

    haversine function is defined as 

                  1 - cos(theta)
    hav(theta) = ---------------- = sin (theta/2) * sin(theta/2)
                        2

    the haversine formula is 

    hav(theta) = hav(lat2 - lat1) + cos(lat1)*cos(lat2)*hav(long2-long1)

    then convert with 
             d
    theta = --- = 2 arcsin(sqrt(hav(theta)))
             r
    """
        
    # convert to radians
    lat1r = np.deg2rad(lat1d)
    lat2r = np.deg2rad(lat2d)

    lon1r = np.deg2rad(lon1d)
    lon2r = np.deg2rad(lon2d)


    # differences
    deltaLatr = lat2r-lat1r
    deltaLonr = lon2r-lon1r

    # compute the angle between two positions on a sphere
    havTheta = (1.-np.cos(deltaLatr)) / 2. \
               + np.cos(lat1r) * np.cos(lat2r) * ((1.-np.cos(deltaLonr))/2.)
        
    thetar = 2. * np.arcsin(np.sqrt(havTheta))

    # returns theta in radians
    return thetar




def get_integral(sigI,sigJ,deltat,b,fL,fH,beta):
    """
    Calculate equation 18 from Siemens+2013
    """

    sigIsigI = sigI*sigI
    sigJsigJ = sigJ*sigJ

    #integrand if sigI==sigJ
    def integrandEqual(f,sig,deltat,b,beta):  
        I = 1./(1.+(sig*sig*deltat)/(b*f**-beta))**2.
        return I 

    if abs(sigI-sigJ)<0.000000001:  
        integral = quad(integrandEqual, fL, fH, args=(sigI,deltat,b,beta))
        #if sigI==sigJ: 
        #    print('same', integral)
        return integral[0] 

    x1 = - (deltat * sigIsigI) / ( b * fH**-beta )
    x2 = - (deltat * sigIsigI) / ( b * fL**-beta )
    x3 = - (deltat * sigJsigJ) / ( b * fH**-beta )
    x4 = - (deltat * sigJsigJ) / ( b * fL**-beta )

    value = (  fH * sigIsigI * sc.hyp2f1(1, beta**-1, 1 + beta**-1, x1) \
             - fL * sigIsigI * sc.hyp2f1(1, beta**-1, 1 + beta**-1, x2) \
             - fH * sigJsigJ * sc.hyp2f1(1, beta**-1, 1 + beta**-1, x3) \
             + fL * sigJsigJ * sc.hyp2f1(1, beta**-1, 1 + beta**-1, x4) ) \
             / ( sigIsigI - sigJsigJ )
    #print(value)
    return float((value))



def get_b(A,fref,alpha):
    """
    computes constant
    """
    b = (A*A) / (24.*np.pi*np.pi) * (1./fref)**(2.*alpha)
    return b



def hellings_downs(angle):
    """
    Hellings-Downs value for angle
    """
    hdCurve = (1./2.) \
              - (1./4.)*((1.-np.cos(angle))/2.) \
              + (3./2.)*((1.-np.cos(angle))/2.)*np.log((1.-np.cos(angle))/2.)

    return hdCurve



def dipole(angle):
    """
    Dipole values
    """ 
    dipValue = (-3./2.) * (np.cos(0) + np.cos(0)) \
               * ( np.cos(angle) - (4./3.) \
                  - 4.*((np.tan(angle/2.))**2.)*np.log(np.sin(angle/2.)) )
    return dipValue




#def avePTASNR(sigmaIs,sigmaJs,angles,A,alpha,beta,fref,T,c):
def avePTASNR(psrNames,psrConstants,hdValues,obsTimes,A,alpha,beta,fref,T,c):

    """ 
    This is from equation 17 of: Siemens+2013
    """
    fL=1./T
    fH=0.5*c #1./c #check
    deltat = 1./c 

    b=get_b(A,fref,alpha)
    total=0
   
    for i,ipsr in enumerate(psrNames):
      for j,jpsr in enumerate(psrNames):
        if (i>j):  # no double counting

          #hd = hellings_downs(angle[ipsr][jpsr])
          hd = hdValues[ipsr][jpsr]

          if obsTimes[ipsr]==0 or obsTimes[jpsr]==0:
            sigI, sigJ=0.0, 0.0
          else:
            sigI = psrConstants[ipsr] / np.sqrt(obsTimes[ipsr])
            sigJ = psrConstants[jpsr] / np.sqrt(obsTimes[jpsr])

            integral = get_integral(sigI,sigJ,deltat, b, fL, fH, beta)
  
            aveSNRSinglePulsarPair = 2.*T*hd*hd*integral 
            total+=aveSNRSinglePulsarPair
            if aveSNRSinglePulsarPair<0: 
              print(sigI,sigJ, integral)
    return np.sqrt(total)

"""

    for sigI,sigJ,ang in zip(sigmaIs,sigmaJs,angles): 

        hd = hellings_downs(ang)

        integral = get_integral(sigI,sigJ,deltat,b,fL,fH,beta)
        #print(count)


        aveSNRSinglePair = 2.*T*hd*hd*integral

        total+=aveSNRSinglePair

    return np.sqrt(total)


"""




# think about whether to use these: 
#def scalingIntermediate(M,c,A,sigma,T,beta):
#    snr = (M*c*A*A*T**beta)/(sigma*sigma)
#    return snr


def scalingIntermediate(angles,sig,T,alpha,beta,fref,c,A):
    T = 1.55*T
    hdTotal = 0
    for ang in angles:
        hdTotal += hellings_downs(ang)
    hdContribution = np.sqrt(hdTotal)

    b=get_b(A,fref,alpha)
    snr = (hdContribution * b * c * T**beta) / (sig*sig * np.sqrt(4.*beta - 2) )
    
    return snr



def scalingLoud(angles,alpha,beta,sig,T,c,fref,A):
    T = 1.55*T
    hdTotal = 0
    for ang in angles:
        hdTotal += hellings_downs(ang)
    hdContribution = np.sqrt(hdTotal)
    b=get_b(A,fref,alpha)
    #print(alpha, ((b*c)/(2.*sig*sig)), 2.*alpha*T * ((b*c)/(2.*sig*sig))**(1./beta))
    #snr = hdContribution * np.sqrt( 2.*alpha*T * ((b*c)/(2.*sig*sig))**(1./beta) )
    bracket1 = (b*c)/(2.*sig*sig)
    bracket2 = alpha * bracket1**(1./beta) * T - 1.
    snr = hdContribution * np.sqrt( 2. * bracket2 )

    return snr
#def scalingLoud(M,c,A,sigma,T,beta):
#    
#    fraction = (c*A*A) / (sigma*sigma)
#    snr = M*fraction**(1./(2.*beta))*T**(1./2.)
#    
#    return snr


def transition(c,A,T,alpha,beta,fref,sigma):
    value = A/(np.pi * fref**alpha) * np.sqrt((c*T**beta)/24.)
    if sigma>value: 
        return 'w'
    elif sigma<value: 
        return 'i'
    return None
