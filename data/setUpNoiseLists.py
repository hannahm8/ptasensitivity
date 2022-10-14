import numpy as np

def getPSRNames():
    data = np.genfromtxt('psrDetails_new.dat',dtype=str)
    psrNames = data[:,0]
    return psrNames


#def readDMNoise():
    
        
def makeDMDict():
    """
    create a dictionary of the puslar DM noise values 
    when the value is unknow, use the median of the rest
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

    # get the median values
    As = np.atleast_1d(As)    
    medianA = np.median(As)
    gs = np.atleast_1d(gs)
    mediang = np.median(gs)

    # put the median values in the dictionaries when noise unknown
    for psr in ADMs.keys():
        if ADMs[psr]==None:
            ADMs[psr] = medianA
        else: pass
        if gDMs[psr]==None:
            gDMs[psr] = mediang
            
    return ADMs, gDMs
    
   
def makeSNDict():
    """
    create a dictionary of the puslar DM noise values 
    when the value is unknow, use the median of the rest
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
   
    # get median values
    As,gs = np.atleast_1d(As), np.atleast_1d(gs)
    medianA = np.median(As)
    mediang = np.median(gs)

    for psr in ASNs.keys():
        if ASNs[psr]==None:
            ASNs[psr] = medianA
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
    
    j = np.atleast_1d(j)
    medianj = np.median(j)
    
    # put the median value in where the value is unknown
    for psr in jitters.keys():
        if jitters[psr]==None:
            jitters[psr] = medianj
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
