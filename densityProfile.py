#!/usr/bin/env python3
#coding utf-8

import pprint
pprint = pprint.PrettyPrinter(indent=1).pprint

from src.parsers.parseData import parseData
from options.options import options, universalOptions




def definePhase(atomNumber, systemName):
    if systemName == 'segregated' or systemName == 'mixed':
        systemSize = 31320 / 9
        while atomNumber > systemSize:
            atomNumber -= systemSize
        if atomNumber < 721:
            return 1
        elif atomNumber < 1561:
            return 2
        return 3
    elif systemName == 'long':
        if atomNumber < 720 * 9 + 1:
            return 1
        #elif 720 * 9 < atomNumber < 720 * 9 + 108*70 + 1:
        elif 720 * 9 < atomNumber < 720 * 9 + 1902 * 9:
            return 2
        return 3
    elif systemName == '5x20':
        systemSize = 1560 + 5 * 382
        while atomNumber > systemSize:
            atomNumber -= systemSize
        if atomNumber < 721:
            return 1
        elif atomNumber < 1561:
            return 2
        return 3
    elif systemName == '10x20':
        systemSize = 1560 + 10 * 382
        while atomNumber > systemSize:
            atomNumber -= systemSize
        if atomNumber < 721:
            return 1
        elif atomNumber < 1561:
            return 2
        return 3

def densityProfileMass(folderNum=1, systemName="mixed"):
    """Define the folder number, range for files"""
    fnameOptions = 'options/' + systemName + '.json'
    o = options(fnameOptions)
    uo = universalOptions()
    uo.multip = 10
    folder = ('/media/anton/Seagate Expansion Drive/Article-MMT/Cluster' + 
              ' calculations for article/1st conf/mixed ok/300K/' + 
              'npt after 6th cycle/1440224/')
    densityProfile = [[0, 0, 0, 0] for i in 
                          range(uo.multip * uo.bigNumber)]
    for i in range(1, 51):
        print(i)
        fname = folder + 'co.' + str(i * 50000) + '.data'
        pd = parseData(fname)
        masses = pd.masses
        atomsFull = pd.atomsFull
        lx = pd.xhi - pd.xlo
        ly = pd.yhi - pd.ylo
        for atom in atomsFull:
            atomNumber = atom[0]
            #z = int(9 + (atom[6] - int(pd.zlo)) * uo.multip)
            z = int((atom[6] + 150) * uo.multip)
            atomType = atom[2]
            mass = masses[atomType - 1]
            phase = definePhase(atomNumber, systemName)
            densityProfile[z][0] += mass * 1.66 # coefficient from units conversion
            densityProfile[z][phase] += mass * 1.66
    f = open('log', 'w')
    for i in range(len(densityProfile)):
            if (densityProfile[i][0] > 0):
                f.write(str(i / uo.multip - 150) + ' ' +
                    str(densityProfile[i][0] / 50 / (lx * ly) * uo.multip) + ' ' +
                    str(densityProfile[i][1] / 50 / (lx * ly) * uo.multip) + ' ' +
                    str(densityProfile[i][2] / 50 / (lx * ly) * uo.multip) + ' ' +
                    str(densityProfile[i][3] / 50 / (lx * ly) * uo.multip) + '\n')
        
densityProfileMass()
