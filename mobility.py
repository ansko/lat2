#!/usr/bin/env python3
#coding utf-8

import copy
import math

import pprint
pprint = pprint.PrettyPrinter(indent=1).pprint

from src.parsers.parseData import parseData
from src.parsers.parseDump import parseDump
from options.options import options, universalOptions
from src.stressToPressure import stressToPressure
from src.fourier.approFourier import appro

folderNum = 6
systemName = "10x20"

def makeMolecule(groupOfAtoms, lx, ly, lz):
    for i in range(len(groupOfAtoms) - 1):
        parentAtom = groupOfAtoms[i]
        nextAtom = groupOfAtoms[i + 1]
        r_min = 10000
        [x_min, y_min, z_min] = [None, None, None]
        for x in [-1, 0, 1]:
            for y in [-1, 0, 1]:
                for z in [-1, 0, 1]:
                    x_parent = parentAtom[4]
                    y_parent = parentAtom[5]
                    z_parent = parentAtom[6]
                    x_next = nextAtom[4] + x * lx
                    y_next = nextAtom[5] + y * ly
                    z_next = nextAtom[6] + z * lz
                    r = ((x_parent - x_next)**2 +
                         (y_parent - y_next)**2 +
                         (z_parent - z_next)**2)
                    if r < r_min:
                        r_min = r
                        [x_min, y_min, z_min] = [x, y, z]
        nextAtom[4] += x_min * lx
        nextAtom[5] += y_min * ly
        nextAtom[6] += z_min * lz
    return groupOfAtoms


def computeChainRange(moleculeNum, systemName):
    if systemName == 'segregated' or systemName == 'mixed':
        first = (3480 * (moleculeNum // 10) + 
                 1560 + 192 * (moleculeNum % 10))
        last = first + 192
    elif systemName == 'long':
        first = 1560 * 9 + (moleculeNum % 10) * 1902
        last = first + 1902
    elif systemName == '5x20' or systemName == '10x20':
        first = (3480 * (moleculeNum // 10) + 
                 1560 + 382 * (moleculeNum % 10))
        last = first + 382
    return (first, last)

def atomNameByCharge(charge):
    if charge == 0.28 or charge == 0.1:
        return 'H'
    elif (charge == 0.02 or charge == -0.2 or
          charge == 0.38 or charge == -0.3):
        return 'C'
    elif charge == -0.56 or charge == -0.5:
        return 'N'
    elif charge == -0.38:
        return 'O'
    print('\n', charge, '\n')

def outputXYZ(groupOfAtoms, fname):
    f = open(fname, 'w')
    falseEnergyString = 'Energy =     -20.32 kcal/mol\n'
    f.write(str(len(groupOfAtoms)) + '\n')
    f.write(falseEnergyString)
    for atom in groupOfAtoms:
        f.write(str(atomNameByCharge(atom[3])) + ' ' + str(atom[4]) +
                ' ' + str(atom[5]) + ' ' + str(atom[6]) + '\n')

def polymerChainsMobility():
    fnameOptions = 'options/' + systemName + '.json'
    o = options(fnameOptions)
    uo = universalOptions()
    folder = o.props['folders'][len(o.props['folders']) - 1][:-6:1]
    cms = [[] for i in range(int(o.props['polymerChainsNum']))]
    #cms[molecule][filenum] = [x, y, z]
    for i in range(1, 51):
        print(i)
        fname = folder + 'co.' + str(i * 50000) + '.data'
        pd = parseData(fname)
        masses = pd.masses
        atomsFull = pd.atomsFull
        lx = pd.xhi - pd.xlo
        ly = pd.yhi - pd.ylo
        lz = pd.zhi - pd.zlo
        for chainNum in range(int(o.props['polymerChainsNum'])):
            (moleculeStart, moleculeEnd) = computeChainRange(chainNum, 
                                               systemName)
            groupOfAtoms = copy.deepcopy(atomsFull[moleculeStart:
                                                   moleculeEnd:1])
            groupOfAtoms = makeMolecule(groupOfAtoms, lx, ly, lz)
            (x_cm, y_cm, z_cm) = computeCM(groupOfAtoms)
            cms[chainNum].append([x_cm, y_cm, z_cm])
            
    for i in range(49):
        xOld = cms[1][i][0]
        yOld = cms[1][i][1]
        zOld = cms[1][i][2]
        xNew = cms[1][i + 1][0]
        yNew = cms[1][i + 1][1]
        zNew = cms[1][i + 1][2]
        xMin = 0
        yMin = 0
        zMin = 0
        rMin = 10000
        for x in [-1, 0, 1]:
            xNewTmp = xNew + x * lx
            r = xNewTmp - xOld
            if r < rMin:
                rMin = r
                xMin = x
        rMin = 10000
        for y in [-1, 0, 1]:
            yNewTmp = yNew + y * ly
            r = yNewTmp - yOld
            if r < rMin:
                rMin = r
                yMin = y
        rMin = 10000
        for z in [-1, 0, 1]:
            zNewTmp = zNew + z * lz
            r = zNewTmp - zOld
            if r < rMin:
                rMin = r
                zMin = z
        xNew += x * lx
        yNew += y * ly
        zNew += z * lz

    deltaCM = 0
    for i in range(int(o.props['polymerChainsNum'])):
        deltaCM += (cms[i][49][0] - cms[i][0][0]) ** 2
        deltaCM += (cms[i][49][1] - cms[i][0][1]) ** 2
        deltaCM += (cms[i][49][2] - cms[i][0][2]) ** 2
    deltaCM = deltaCM ** 0.5
    deltaCM /= int(o.props['polymerChainsNum'])
    print(deltaCM)


def computeCM(groupOfAtoms):
    x_cm = 0
    y_cm = 0
    z_cm = 0
    for atom in groupOfAtoms:
        x_cm += atom[4]
        y_cm += atom[5]
        z_cm += atom[6]
    x_cm /= len(groupOfAtoms)
    y_cm /= len(groupOfAtoms)
    z_cm /= len(groupOfAtoms)
    return (x_cm, y_cm, z_cm)


def atomicMobility():
    multip = 1
    fnameOptions = 'options/' + systemName + '.json'
    o = options(fnameOptions)
    uo = universalOptions()
    folder = o.props['folders'][0][:-6:1]
    pdStart = parseData(folder + 'co.' + str(50000) + '.data')
    pdFinish = parseData(folder + 'co.' + str(2500000) + '.data')
    masses = pdStart.masses
    atomsFullStart = pdStart.atomsFull
    atomsFullFinish = pdFinish.atomsFull
    lx = pdStart.xhi - pdStart.xlo
    ly = pdStart.yhi - pdStart.ylo
    lz = pdStart.zhi - pdStart.zlo
    mobilityProfile = [[0, 0] for i in range(300 * multip)]
    for i in range(len(atomsFullStart)):
        if masses[atomsFullStart[i][2] - 1] < 5:
            continue
        z = int(100 + multip * atomsFullStart[i][6])
        dx = min(abs(atomsFullFinish[i][4] - atomsFullStart[i][4]),
                 lx - abs(atomsFullFinish[i][4] - atomsFullStart[i][4]))
        dy = min(abs(atomsFullFinish[i][5] - atomsFullStart[i][5]),
                 ly - abs(atomsFullFinish[i][5] - atomsFullStart[i][5]))
        dz = min(abs(atomsFullFinish[i][6] - atomsFullStart[i][6]),
                 lz - abs(atomsFullFinish[i][6] - atomsFullStart[i][6]))
        mobilityProfile[z][0] += 1
        mobilityProfile[z][1] += (dx**2 + dy**2 + dz**2)**0.5
    f = open('log', 'w')
    for i in range(len(mobilityProfile)):
        if mobilityProfile[i][0] != 0:
            f.write(str((i - 100) / multip) + ' ' + str(mobilityProfile[i][0]) +
                    ' ' + str(mobilityProfile[i][1] / mobilityProfile[i][0]) + '\n')


#polymerChainsMobility()
atomicMobility()
