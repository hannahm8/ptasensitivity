import numpy as np
import scalingRelFunctions

import matplotlib.pyplot as plt



# loop over the T's to plot snr against time


# observation time in years
Ts = np.logspace(0,1.5,100)
TInSeconds = Ts * (365.25*24.*60.*60)

# cadence
c = 20. / (365.25*24.*60.*60.)



# signal details
A=1.E-15
beta=13./3
alpha=(3.-beta)/2.
fref = 1./(365.25*24.*60.*60.)



# number of pulsar pairs - will improve later
nPairs = 190


# random set of properties for the pulsar pairs 
angles = np.random.rand(nPairs)*(np.pi/2.)
sigmaIs = 1E-7 + np.random.rand(nPairs)*1.1E-7
sigmaJs = 1E-7 + np.random.rand(nPairs)*1.1E-7


snrs = np.zeros(len(Ts))

for i, Ti in enumerate(TInSeconds):


    snrs[i] = scalingRelFunctions.avePTASNR(sigmaIs,sigmaJs,angles,
                                            A, alpha, beta, fref,Ti,c)

    
M=20
snrI = [ scalingRelFunctions.scalingIntermediate(M,c,A,1E-7,Ti,beta) for Ti in TInSeconds] 

snrL = [ scalingRelFunctions.scalingLoud(M,c,A,1E-7,Ti,beta) for Ti in TInSeconds]


check='i'
for Ti in TInSeconds:
    transition=scalingRelFunctions.transition(c,A,Ti,alpha,beta,fref,1E-7)
    if transition==check:   
        transitionTime=Ti



# next sort out the angles situation
scaleI = snrI[0]/snrs[0]
scaleL = snrL[-1]/snrs[-1]

print(transitionTime/(365.25*24.*60.*60), scaleI, scaleL)

plt.loglog(Ts,snrs)
#plt.loglog(Ts,snrL)#/scaleL)
#plt.loglog(Ts,snrI)#/scaleI)
plt.axvline(transitionTime)
plt.xlim(1,30)
plt.ylim(1E-2,1E2)
plt.show()
