import numpy as np
"""
based off of original psrDetails.dat, this will update the 
integration time and expected precision if we allocate more
time to the array

Currently set to increase integration time by 20%
"""

psrNames = np.genfromtxt('psrDetails.dat',usecols=0,dtype=str)
data = np.genfromtxt('psrDetails.dat',names=True)

outFile = 'psrDetailsExtraTime1p2.dat'
extraTimeFile = open(outFile,'w')
extraTimeFile.write("#PSR\tRA\tDEC\tIntTime\tExpPrecision\n")

for i,psr in enumerate(psrNames):

    sigmaOriginal = data['ExpPrecision'][i]*1.E-6
 
    constant =  sigmaOriginal * np.sqrt(data['IntTime'][i]) 

    intTimeNew = data['IntTime'][i]*1.2
    sigmaNew = constant / np.sqrt(intTimeNew)


    extraTimeFile.write("{}\t{}\t{}\t{}\t{}\n".format(psr,\
                                                      data['RA'][i],\
                                                      data['DEC'][i],\
                                                      intTimeNew,\
                                                      sigmaNew/1.E-6))

extraTimeFile.close()

