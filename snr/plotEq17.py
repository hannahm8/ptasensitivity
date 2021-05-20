import numpy as np
import scalingRelFunctions

import matplotlib.pyplot as plt





# loop over the T's to plot snr against time


# observation time in years
Ts = np.logspace(0,2.5,100)
TInSeconds = Ts * (365.25*24.*60.*60)

# cadence
c = 20. / (365.25*24.*60.*60.)



# signal details
A=1.E-15
beta=13./3
alpha=(3.-beta)/2.
fref = 1./(365.25*24.*60.*60.)



# number of pulsar pairs - will improve later
#nPairs = 190


# get pulsar pair data ->
import angleBetweenPSRs
angles, sigmaIs, sigmaJs = angleBetweenPSRs.get_combos_as_is()
sigmaIs = sigmaIs
sigmaJs = sigmaJs

snrs = np.zeros(len(Ts))

for i, Ti in enumerate(TInSeconds):


    snrs[i] = scalingRelFunctions.avePTASNR(sigmaIs,sigmaJs,angles,
                                            A, alpha, beta, fref,Ti,c)


# senario 2 -> using 20 'best' observed for four times as long
howMany2=35
tf2=2
angles, sigmaIs, sigmaJs = angleBetweenPSRs.get_combos_x_best(nToUse=howMany2,timeFactor=tf2)

snrs2 = np.zeros(len(Ts))

for i, Ti in enumerate(TInSeconds):

    snrs2[i] = scalingRelFunctions.avePTASNR(sigmaIs,sigmaJs,angles,
                                            A, alpha, beta, fref,Ti,c)


# senario 2 -> using 20 'best' observed for four times as long
howMany3=15
tf3=4
angles, sigmaIs, sigmaJs = angleBetweenPSRs.get_combos_x_best(nToUse=howMany3,timeFactor=tf3)

snrs3 = np.zeros(len(Ts))

for i, Ti in enumerate(TInSeconds):

    snrs3[i] = scalingRelFunctions.avePTASNR(sigmaIs,sigmaJs,angles,
                                            A, alpha, beta, fref,Ti,c)


"""    
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
"""


#print(transitionTime/(365.25*24.*60.*60), scaleI, scaleL)
print(snrs)

plt.loglog(Ts,snrs, label='current')
plt.loglog(Ts,snrs2,label='{} best {}T'.format(howMany2,tf2))
plt.loglog(Ts,snrs3,label='{} best {}T'.format(howMany3,tf3))
#plt.loglog(Ts,snrL)#/scaleL)
#plt.loglog(Ts,snrI)#/scaleI)
#plt.axvline(transitionTime)
#plt.xlim(1,30)
#plt.ylim(1E-2,1E2)
plt.ylabel('SNR')
plt.xlabel('Time (years)')
plt.grid()
plt.legend()
plt.show()
