import numpy as np



#catalog = np.genfromtxt('psrCatElatElong.dat',skip_header=5,dtype=str)
catalog = np.genfromtxt('psrCatRADDECD.dat',skip_header=74,dtype=str)
catPSRNames = catalog[:,1]

ras   = [ float(catalog[i,3]) for i in range(len(catalog[:,3])) ]
decs  = [ float(catalog[i,5]) for i in range(len(catalog[:,5])) ]



meertimePSRs = np.genfromtxt('meerTimePSRsv2.dat',dtype=str,skip_header=2)



psrListOut = open('psrDetails.dat','w')
psrListOut.write('#PSR\tRA\tDEC\tIntTime\tExpPrecision\n')
trueCount, falseCount = 0, 0

for psrName,intTime,expPrecision in zip(meertimePSRs[:,0], \
                                        meertimePSRs[:,1], \
                                        meertimePSRs[:,2]):

    # check if pulsar in psrcat list
    if psrName in catPSRNames: 
    
        element=np.where(catPSRNames==psrName)[0][0]
        psrListOut.write('{}\t{}\t{}\t{}\t{}\n'.format(psrName, \
                                                           ras[element], \
                                                           decs[element], \
                                                           intTime, \
                                                           expPrecision))

        trueCount+=1
    else: 
        print(psrName, 'False')
        falseCount+=1

psrListOut.close()
print('''
found: {}
couldn't find: {}
'''.format(trueCount, falseCount))
