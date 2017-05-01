#!/usr/bin/env python3
#coding utf-8


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


def returnCalyRanges(systemName):
    if systemName == 'mixed':
        return (4.8, 12.5)
    elif systemName == 'segregated':
        return (51.7, 59.4)

def modifierDistribution(systemName='mixed'):
    """Change systemName, ParentFolders and Folders and start and end"""
    segregatedParentFolder = ("/media/anton/Seagate Expansion Drive/Article-MMT/" +
                             "Cluster calculations for article/1st conf/" + 
                             "segregated/500K/")
    segregatedFolders = ["1332494 - 2.5ns/",
                         "1338014 - 2.5ns more/", 
                         "1343289 - 2.5ns more/",
                         "1352547 - 2.5 ns more/",
                         "1388142 - 2.5 ns more/", 
                         "1404644/"]
    mixedParentFolder = ("/media/anton/Seagate Expansion Drive/Article-MMT/" +
                         "Cluster calculations for article/1st conf/" + 
                         "mixed ok/500K/")
    mixedFolders = ["all/",
                    "1388141/", 
                    "1393565/",
                    "1397381/",
                    "1400005/",
                    "1404643/"]
    fnameOptions = 'options/' + systemName + '.json'
    o = options(fnameOptions)
    uo = universalOptions()
    folder = o.props['folders'][len(o.props['folders']) - 1][:-6:1]
    [top, bottom] = returnCalyRanges(systemName)
    filenum = 0
    f = open('log' ,'w')
    for j in range(len(mixedFolders)):
        end = 51
        start = 1
        if j == 0:
           end = 24
        if j == 2:
            end = 50
        for i in range(start, end):
            filenum += 1
            fname = (mixedParentFolder +
                     mixedFolders[j] + 'co.' + str(i * 50000) + '.data')
            pd = parseData(fname)
            lz = pd.zhi - pd.zlo
            atomsFull = pd.atomsFull
            distanceSum = 0
            atomsNum = 0
            for atom in atomsFull:
                if definePhase(atom[0], systemName) == 2:
                    atomsNum += 1
                    distance = min((atom[6] - top)**2,
                                   (lz + bottom - atom[6])**2)
                    distanceSum += distance
           
            f.write(str(filenum) + ' ' + str((distanceSum / atomsNum)**0.5) + '\n')

modifierDistribution()
