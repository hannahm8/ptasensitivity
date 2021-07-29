import numpy as np
import snrFunctions

import matplotlib.pyplot as plt



no longer works!!!! Do not use!! 

# loop over the T's to plot snr against time


# observation time in years
Ts = np.logspace(0,2.5,20)
Ts = np.linspace(0,15,50)
TInSeconds = Ts * (365.25*24.*60.*60)

# cadence
c = 26. / (365.25*24.*60.*60.)



# signal details
A=2.E-15
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

scaleInt = np.zeros(len(Ts))
scaleLoud = np.zeros(len(Ts))

for i, Ti in enumerate(TInSeconds):


    snrs[i] = snrFunctions.avePTASNR(sigmaIs,sigmaJs,angles,
                                            A, alpha, beta,fref,Ti,c)

    scaleInt[i] = snrFunctions.scalingIntermediate(angles,1.E-6,Ti,alpha,beta,fref,c,A)


    scaleLoud[i] = snrFunctions.scalingLoud(angles,alpha,beta,1.E-6,Ti,c,fref,A)

# senario 2 -> using 20 'best' observed for four times as long
howMany2=35
tf2=2
angles, sigmaIs, sigmaJs = angleBetweenPSRs.get_combos_x_best(nToUse=howMany2,timeFactor=tf2)

snrs2 = np.zeros(len(Ts))

for i, Ti in enumerate(TInSeconds):

    snrs2[i] = snrFunctions.avePTASNR(sigmaIs,sigmaJs,angles,
                                            A, alpha, beta, fref,Ti,c)



# senario 2 -> using 20 'best' observed for four times as long
howMany3=15
tf3=4
angles, sigmaIs, sigmaJs = angleBetweenPSRs.get_combos_x_best(nToUse=howMany3,timeFactor=tf3)

snrs3 = np.zeros(len(Ts))

for i, Ti in enumerate(TInSeconds):

    snrs3[i] = snrFunctions.avePTASNR(sigmaIs,sigmaJs,angles,
                                            A, alpha, beta, fref,Ti,c)


"""    
M=20
snrI = [ scalingRelFunctions.scalingIntermediate(M,c,A,1E-7,Ti,beta) for Ti in TInSeconds] 

snrL = [ scalingRelFunctions.scalingLoud(M,c,A,1E-7,Ti,beta ) for Ti in TInSeconds]


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

plt.title('Using cadence = {} per year'.format(int(c*(365.25*24.*60.*60.))))
#plt.loglog(Ts,snrs, label='Full spreadsheet list (minus ~5 )')
#plt.loglog(Ts,snrs2,label='{} PSRs observed for {}x as long'.format(howMany2,tf2))
#plt.loglog(Ts,snrs3,label='{} PSRs observed for {}x as long'.format(howMany3,tf3))

plt.plot(Ts,snrs, label='Full spreadsheet list (minus ~5 )')
plt.plot(Ts,scaleInt,label='scaling rel intermediate')
plt.plot(Ts,scaleLoud,label='scaling rel Loud')
#plt.plot(Ts,snrs2,label='{} PSRs observed for {}x as long'.format(howMany2,tf2))
#plt.plot(Ts,snrs3,label='{} PSRs observed for {}x as long'.format(howMany3,tf3))
print(scaleLoud)
#plt.loglog(Ts,snrL)#/scaleL)
#plt.loglog(Ts,snrI)#/scaleI)
#plt.axvline(transitionTime)
#plt.xlim(1,30)
#plt.ylim(1E-2,1E2)
plt.yscale('log')
plt.xscale('log')
plt.ylim(0,100)
plt.ylabel('SNR')
plt.xlabel('Time (years)')
plt.grid()
plt.legend()
plt.savefig('SNRVTime_2E-15_c{}_15yrs_all.png'.format(int(c*(365.25*24.*60.*60.))))
plt.show()
