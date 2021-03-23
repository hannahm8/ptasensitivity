import numpy as np
import scipy.special as sc

"""
 Class. Quantum Grav. 30 224015
"""


def get_integral(sigI,sigJ,deltat,b,fL,fH,beta):
    """
    Calculate equation 18 from Siemens+2013
    """

    if sigI==sigJ: return 0 # check this is okay 

    x1 = - (deltat * sigI * sigI) / ( b * fH**-beta )
    x2 = - (deltat * sigI * sigI) / ( b * fL**-beta )
    x3 = - (deltat * sigJ * sigJ) / ( b * fH**-beta )
    x4 = - (deltat * sigJ * sigJ) / ( b * fL**-beta )

    value = (  fH * sigI*sigI * sc.hyp2f1(1, beta**-1, 1 + beta**-1, x1) \
             - fL * sigI*sigI * sc.hyp2f1(1, beta**-1, 1 + beta**-1, x2) \
             - fH * sigJ*sigJ * sc.hyp2f1(1, beta**-1, 1 + beta**-1, x3) \
             + fL * sigJ*sigJ * sc.hyp2f1(1, beta**-1, 1 + beta**-1, x4) ) \
             / ( sigI*sigI - sigJ*sigJ )
    
    return value



def get_b(A,fref,alpha):
    b = (A*A) / (24.*np.pi*np.pi) * (1./fref)**(2.*alpha)
    return b



def hellings_downs(angle):
    hdCurve = (1./2.) \
              - (1./4.)*((1.-np.cos(angle))/2.) \
              + (3./2.)*((1.-np.cos(angle))/2.)*np.log((1.-np.cos(angle))/2.)
    return hdCurve




def avePTASNR(sigmaIs,sigmaJs,angles,A,alpha,beta,fref,T,c):
    """ 
    This is from equation 17 of: Siemens+2013
    """
    fL=1./T
    fH=1./c #check
    deltat = 1./c #check

    b=get_b(A,fref,alpha)

    total=0
   
    count = 0
    for sigI,sigJ,ang in zip(sigmaIs,sigmaJs,angles): 

        hd = hellings_downs(ang)

        integral = get_integral(sigI,sigJ,deltat,b,fL,fH,beta)
        #print(count)
        count+=1

        aveSNRSinglePair = 2.*T*hd*hd*integral

        total+=aveSNRSinglePair

    return np.sqrt(total)







# think about whether to use these: 
def scalingIntermediate(M,c,A,sigma,T,beta):

    snr = (M*c*A*A*T**beta)/(sigma*sigma)

    return snr


def scalingLoud(M,c,A,sigma,T,beta):
    
    fraction = (c*A*A) / (sigma*sigma)
    snr = M*fraction**(1./(2.*beta))*T**(1./2.)
    
    return snr


def transition(c,A,T,alpha,beta,fref,sigma):
    value = A/(np.pi * fref**alpha) * np.sqrt((c*T**beta)/24.)
    if sigma>value: 
        return 'w'
    elif sigma<value: 
        return 'i'
    return None
