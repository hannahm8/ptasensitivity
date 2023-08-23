import numpy as np

def getPSRNames():
    data = np.genfromtxt('psrDetails.dat',dtype=str)
    psrNames = data[:,0]
    return psrNames


    
        
def makeDMDict():
    """
    create a dictionary of the puslar DM noise values 
    when the value is unknow, use the minimum of the rest
    for the noise amplitude and the median of the rest for
    the noise slope. 
    """
    # set up dictionaries
    ADMs, gDMs = {}, {}
    
    # set up lists for existing values (used to make the median)
    As, gs = [], []

    # open file and read first line
    f = open('dmNoiseWithErrors.dat','r')
    line = f.readline().split()

    # set up dictionary for DM noise values 
    while line:

        if len(line)>1 and line[0][0]=='J': 
            ADMs[line[0]]=float(line[1])
            gDMs[line[0]]=float(line[4])

            As.append(float(line[1]))
            gs.append(float(line[4]))
            
        elif line[0][0]=='J': 
            ADMs[line[0]]=None
            gDMs[line[0]]=None

        # next line 
        line = f.readline().split()

    # get the minimum value for A
    As = np.atleast_1d(As)    
    minimumA = min(As)
    
    # get the median value for gamma
    gs = np.atleast_1d(gs)
    mediang = np.median(gs)

    # put the minimum/median values in the dictionaries when noise unknown
    # minimum for A
    # median for gamma 
    for psr in ADMs.keys():
        if ADMs[psr]==None:
            ADMs[psr] = minimumA
        else: pass
        if gDMs[psr]==None:
            gDMs[psr] = mediang
            
    return ADMs, gDMs
    
   
def makeSNDict():
    """
    create a dictionary of the puslar spin/red noise values 
    when the value is unknow, use the minium for the amplitude
    and the median for the slope
    """
    ASNs, gSNs = {}, {}
    As, gs = [], []
 
    f = open('redNoiseWithErrors.dat','r')
    line = f.readline().split()
    
    while line:
    
        if len(line)>1 and line[0][0]=='J':
            ASNs[line[0]]=float(line[1])
            gSNs[line[0]]=float(line[4])

            As.append(float(line[1]))
            gs.append(float(line[4]))
            
        elif line[0][0]=='J': 
            ASNs[line[0]]=None
            gSNs[line[0]]=None
        
        # next line 
        line = f.readline().split()
   
    # get minimum value for A
    As = np.atleast_1d(As)
    minimumA = min(As)

    # get median for gamma
    gs = np.atleast_1d(gs)
    mediang = np.median(gs)

    # set the default (minimum or median) values for psrs with no 
    # measured noise
    for psr in ASNs.keys():
        if ASNs[psr]==None:
            ASNs[psr] = minimumA
        else: pass
        
        if gSNs[psr]==None:
            gSNs[psr] = mediang
        else: pass
        
    return ASNs, gSNs


def makeJitterDict(): 
    """ 
    create a dictionary for the jitter noise
    """
    
    jitters = {}
    j = []
    
    f = open('jitterNoiseWithErrors.dat','r')
    line = f.readline().split()
    
    while line:
    
        if len(line)>1 and line[0][0]=='J' and line[2]!='T':
            jitters[line[0]] = float(line[1])
            j.append(float(line[1]))

        elif line[0][0]=='J':
            jitters[line[0]] = None
            
        # next line
        line = f.readline().split()
    
    # get the minimum jitter value
    j = np.atleast_1d(j)
    minimumj = min(j)
    
    # use the minimum value where the jitter is unknown
    for psr in jitters.keys():
        if jitters[psr]==None:
            jitters[psr] = minimumj
        else: pass
        
    return jitters
    
    

psrNames = getPSRNames()



sn_As, sn_gs = makeSNDict()
dm_As, dm_gs = makeDMDict()
js = makeJitterDict()

snFile = open('realistic/redNoise.dat','w')
snFile.write('#PSR\tASN\tgammaSN\n')
dmFile = open('realistic/dmNoise.dat','w')
dmFile.write('#PSR\tADM\tgammaDM\n')
jitFile = open('realistic/jitterNoise.dat','w')
jitFile.write('#PSR\tjitter\n')

for psr in psrNames:
    snFile.write('{}\t{}\t{}\n'.format(psr,sn_As[psr],sn_gs[psr]))
    dmFile.write('{}\t{}\t{}\n'.format(psr,dm_As[psr],dm_gs[psr]))
    jitFile.write('{}\t{}\n'.format(psr,js[psr]))
snFile.close()
dmFile.close()
jitFile.close()









    
"""


f = open('dmNoiseWithErrors.dat','r')
line = f.readline().split()
while line:

    if len(line)>1: pass
    else: print(line)

    line = f.readline().split()
    
f.close()

"""
