import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.append('/home/hannahm/repositories/ptasensitivity/snr/')
sys.path.append('/fred/oz005/users/hmiddlet/ptasensitivity/snr/')
sys.path.append('/home/ADF/middlehr/repositories/ptasensitivity/snr/')
import snrFunctions
import readInData

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Helvetica",
    "font.size":18,
})

# run like this:
# python plotCorrelationFunctions.py ../data/psrDetails.dat ../data/realistic/redNoise.dat ../data/realistic/dmNoise.dat ../data/realistic/jitterNoise.dat 

originalFile = sys.argv[1]
redNoiseFile = sys.argv[2]
dmNoiseFile  = sys.argv[3]
jitNoiseFile = sys.argv[4]
    
whichCorrelationFunction='HD'

psrNames, \
psrObsConstants, \
psrStartingObsTimes, \
angles, \
hdValues, \
ampRed, \
gammaRed, \
ampDM, \
gammaDM, \
jitterNoise = readInData.readDataIntoDicts(originalFile, \
                                        whichCorrelationFunction, \
                                        redNoiseFile=redNoiseFile, \
                                        dmNoiseFile=dmNoiseFile, \
                                        jitterNoiseFile=jitNoiseFile)

time=0
for psr in psrNames:
    time+=psrStartingObsTimes[psr]
print(time/(60*60))


exit()                                        

# get the correlation functions for a range of angles
step=0.01
plottingAngles = np.arange(step,np.pi-step,step)
anglesInDegrees = [ a*(180./np.pi) for a in plottingAngles]

hd = snrFunctions.hellings_downs(plottingAngles)
dp = snrFunctions.dipole(plottingAngles)
dpMinus2hd = dp - (2.*hd)


orColour, hdColour, dpColour, eqColour = '#0072B2', '#E69F00', '#009E73', '#D55E00'
orLine, hdLine, dpLine, eqLine = 'solid','dotted','dashed','dashdot'


fig = plt.figure(figsize=(10,9))
ax = fig.add_subplot()
ax.plot(anglesInDegrees,hd,
         color=hdColour,lw=3,ls=hdLine,zorder=300,label='Hellings-Downs (HD)')
ax.plot(anglesInDegrees,dpMinus2hd,
         color=dpColour,lw=3,ls=dpLine,zorder=200,label='DP-2HD')
ax.plot(anglesInDegrees,dp,
         color='k',lw=3,alpha=0.25,zorder=100,label='Dipole')

for psrI in psrNames: 
    for psrJ in psrNames:
   
        if psrI==psrJ:
            pass
        else: 
            angle = angles[psrI][psrJ]
            angleDeg = angle*(180./np.pi)
            ax.scatter(angleDeg,-0.4,marker='|',color='k',s=200,alpha=0.01)
plt.legend()    
plt.xlim(0,180)
plt.ylabel('Correlation',fontsize=22)
plt.xlabel('Angle (degrees)',fontsize=22) 
plt.tight_layout()
plt.savefig('correlationFunctions.pdf')
plt.savefig('correlationFunctions.png')
plt.show()


