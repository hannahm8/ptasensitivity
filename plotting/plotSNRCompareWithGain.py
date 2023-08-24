import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Helvetica",
    "font.size":18,
})

def getGain(original,new):
    difference = (new-original)
    percentage = 100.*(difference/original)
    return percentage


original = np.genfromtxt('plots/rerun-v2/original_snrVtime.dat')
DP2HD = np.genfromtxt('plots/rerun-v2/DP-2HD_snrVtime.dat')
HD = np.genfromtxt('plots/rerun-v2/HD_snrVtime.dat')
Equal = np.genfromtxt('plots/rerun-v2/Equal_snrVtime.dat')


originalT, originalSNR = original[:,0], original[:,1]
HDT, HDSNR             = HD[:,0], HD[:,1]
DP2HDT, DP2HDSNR       = DP2HD[:,0], DP2HD[:,1]
EqualT, EqualSNR       = Equal[:,0], Equal[:,1]



# get gains
HDPercentage = getGain(originalSNR,HDSNR)
DP2HDPercentage = getGain(originalSNR,DP2HDSNR)
EqualPercentage = getGain(originalSNR,EqualSNR)

print(DP2HDPercentage)
orColour, hdColour, dpColour, eqColour = '#0072B2', '#E69F00', '#009E73', '#D55E00'
orLine, hdLine, dpLine, eqLine = 'solid','dotted','dashed','dashdot'



fig = plt.figure(figsize=(10,9))
gs = gridspec.GridSpec(2, 1, width_ratios=[1], height_ratios=[1,1/3], wspace=0, hspace=0)

for i,g in enumerate(gs):

    ax = fig.add_subplot(g)    
    if i==0:

        ax.plot(originalT,originalSNR,
                ls=orLine, lw=3, color=orColour, label='Original')
        ax.plot(HDT,HDSNR,
                ls=hdLine, lw=3, color=hdColour, label='HD')
        ax.plot(DP2HDT,DP2HDSNR,
                ls=dpLine, lw=3, color=dpColour, label='DP-2HD')
        ax.plot(EqualT,EqualSNR,
                ls=eqLine, lw=3, color=eqColour, label='Equal')
        ax.legend()
        ax.set_ylabel('Average PTA signal to noise ratio', fontsize=22)
        ax.tick_params(labelbottom = False, bottom = False)

    if i==1:
        ax.plot(HDT, HDPercentage, 
                ls=hdLine, lw=3, color=hdColour, label='HD')
        ax.plot(DP2HDT, DP2HDPercentage,
                ls=dpLine, lw=3, color=dpColour, label='DP-2HD')
        ax.plot(EqualT, EqualPercentage,
                ls=eqLine, lw=3, color=eqColour, label='Equal')
        ax.set_ylim(0,95)
        ax.legend()
        ax.set_ylabel('Percentage gain', fontsize=22)
        ax.set_xlabel('Time (years)', fontsize=22)

    ax.set_xlim(1,10)
plt.tight_layout()
plt.savefig('compareSNRWithGain.pdf')
plt.savefig('compareSNRWithGain.png')
plt.show()
