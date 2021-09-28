import numpy as np
import cpnest.model
import sys

sys.path.append('/home/hannahm/repositories/ptasensitivity/snr/')
sys.path.append('/fred/oz005/users/hmiddlet/ptasensitivity/snr/')
import snrFunctions
import readInData


class snrTime(cpnest.model.Model):


    def __init__(self):

        psrDataFile = '../data/psrDetails.dat'
        redNoisePath = '../data/redNoise.dat'
        jitterNoisePath = '../data/jitterNoise.dat'
        chooseCorrFunc = 'HD'

        # read in the data
        self.psrNames, \
        self.psrConstants, \
        self.psrStartingObsTimes, \
        angles, \
        self.angCorrValues, \
        self.redAmps, \
        self.redGammas, \
        self.jitters = readInData.readDataIntoDicts(psrDataFile,\
                                               chooseCorrFunc,\
                                               redNoiseFile=redNoisePath,\
                                               jitterNoiseFile=jitterNoisePath)
        self.names=[]
        for pname in self.psrNames: 
            self.names.append(str(pname))
        self.dim = len(self.names)

        oneYearInSeconds = (365.25*24.*60.*60.)
        T = 10.
        self.TInSeconds = T * oneYearInSeconds
        self.A = 2.E-15
        self.beta = 13./3
        self.alpha = (3.-self.beta)/2.
        self.fref = 1./oneYearInSeconds
        self.c = 26./oneYearInSeconds

        self.totalTime = 44096.

        self.bounds=[]

        minT,maxT = 256, 2560
        for _ in range(self.dim):
            self.bounds.append([minT,maxT])



    def getValues(self,p):
        theta = ((p[str(self.names[0])]),)
        for i in range(self.dim-1):
            theta = theta + ((p[str(self.names[i+1])]),)
        return theta



    def log_prior(self, p):
        if not self.in_bounds(p): return -np.inf
        theta = self.getValues(p)
        if sum(theta)>self.totalTime: return -np.inf
        return 0
        

    def log_likelihood(self,param):

        #theta = ((param[str(self.names[0])]),)
        #for i in range(self.dim-1):
        #    theta = theta + ((param[str(self.names[i+1])]),)
        theta = self.getValues(param)


        obsTimes = {}
        for t,psr in zip(theta,self.names):
            obsTimes[psr] = i

    
        snr = snrFunctions.avePTASNR(self.psrNames,\
                                     self.psrConstants,\
                                     self.angCorrValues, \
                                     obsTimes, \
                                     self.redAmps,\
                                     self.redGammas,\
                                     self.jitters,\
                                     self.A,\
                                     self.alpha,\
                                     self.beta,\
                                     self.fref,\
                                     self.TInSeconds,\
                                     self.c)
        return snr


mymodel = snrTime()
#nest = cpnest.CPNest(mymodel,maxmcmc=1000,nlive=10000,verbose=3,nthreads=1)
nest = cpnest.CPNest(mymodel,maxmcmc=100,nlive=500,verbose=3,nthreads=8)
nest.run()
cpnest.CPNest.get_posterior_samples(nest)

