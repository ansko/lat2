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
        first = 720 * 9 + (moleculeNum % 9) * 1902
        last = first + 1902
    elif systemName == '10x20':
        first = (3480 * (moleculeNum // 10) + 
                 1560 + 382 * (moleculeNum % 10))
        last = first + 382
    elif systemName == '5x20':
        first = (3480 * (moleculeNum // 5) + 
                 1560 + 382 * (moleculeNum % 5))
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
    #print('\n', charge, '\n')

def outputXYZ(groupOfAtoms, fname):
    f = open(fname, 'w')
    falseEnergyString = 'Energy =     -20.32 kcal/mol\n'
    f.write(str(len(groupOfAtoms)) + '\n')
    f.write(falseEnergyString)
    for atom in groupOfAtoms:
        f.write(str(atomNameByCharge(atom[3])) + ' ' + str(atom[4]) +
                ' ' + str(atom[5]) + ' ' + str(atom[6]) + '\n')
    f.close()


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


def computeRGyr(groupOfAtoms):
    (x_cm, y_cm, z_cm) = computeCM(groupOfAtoms)
    r_x = r_y = r_z = 0
    N = len(groupOfAtoms)
    for atom in groupOfAtoms:
        r_x += (atom[4] - x_cm)**2
        r_y += (atom[5] - y_cm)**2
        r_z += (atom[6] - z_cm)**2
    r_x /= N
    r_y /= N
    r_z /= N
    return (r_x ** 0.5, r_y ** 0.5, r_z ** 0.5)


def computeRgyr(molecule):
    n = len(molecule)
    x = y = z = 0
    rx = ry = rz = 0
    errx = erry = errz = 0

    for i in range(n):
        x += molecule[i][4]
        y += molecule[i][5]
        z += molecule[i][6]
    x /= n
    y /= n
    z /= n

    for atom in molecule:
        rx += (atom[4] - x)**2
        ry += (atom[5] - y)**2
        rz += (atom[6] - z)**2

    rx /= n
    ry /= n
    rz /= n

    for atom in molecule:
        errx += (abs(atom[4] - x) - math.sqrt(rx))**2
        erry += (abs(atom[5] - y) - math.sqrt(ry))**2
        errz += (abs(atom[6] - z) - math.sqrt(rz))**2

    return [rx**0.5, ry**0.5, rz**0.5, errx / n, erry / n, errz / n]


def rgyr(folderNum=6, systemName = "10x20"):
    fnameOptions = 'options/' + systemName + '.json'
    o = options(fnameOptions)
    uo = universalOptions()

    r_xVSz = [[0, 0] for i in range(1500)]
    r_yVSz = [[0, 0] for i in range(1500)]
    r_zVSz = [[0, 0] for i in range(1500)]
    for i in range(len(o.props['folders'])):
        print(i)
        folder = o.props['folders'][i][:-6:1]
        fname = folder + 'co.' + str(1 * 50000) + '.data'
        pd = parseData(fname)
        masses = pd.masses
        atomsFull = pd.atomsFull
        lx = pd.xhi - pd.xlo
        ly = pd.yhi - pd.ylo
        lz = pd.zhi - pd.zlo
        for chainNum in range(int(o.props['polymerChainsNum'])):
            (moleculeStart, moleculeEnd) = computeChainRange(chainNum, 
                                               systemName)
            groupOfAtoms = copy.deepcopy(atomsFull[moleculeStart:moleculeEnd:1])
            groupOfAtoms = makeMolecule(groupOfAtoms, lx, ly, lz)
            (rx, ry, rz) = computeCM(groupOfAtoms)
            outputXYZ(groupOfAtoms, 'pics/' + systemName + str(chainNum) + '.xyz')
            (rGyr_x, rGyr_y, rGyr_z, ex, ey, ez) = computeRgyr(groupOfAtoms)
            z = int((550 + rz) * uo.multip)
            r_xVSz[z][0] += 1
            r_xVSz[z][1] += rGyr_x
            r_yVSz[z][0] += 1
            r_yVSz[z][1] += rGyr_y
            r_zVSz[z][0] += 1
            r_zVSz[z][1] += rGyr_z

    f = open('log', 'w')
    for i in range(len(r_zVSz)):
        if r_zVSz[i][0] != 0:
            f.write(str(i / uo.multip - 150) + ' ' +
                    str(r_zVSz[i][1] / r_zVSz[i][0]) + '\n')

    r_x_mean = 0
    valuesNum = 0
    for i in range(len(r_xVSz)):
        if r_xVSz[i][0] != 0:
            r_x_mean += r_xVSz[i][1]
            valuesNum += r_xVSz[i][0]
    print(r_x_mean / valuesNum)
    r_y_mean = 0
    valuesNum = 0
    for i in range(len(r_yVSz)):
        if r_yVSz[i][0] != 0:
            r_y_mean += r_yVSz[i][1]
            valuesNum += r_yVSz[i][0]
    print(r_y_mean / valuesNum)
    r_z_mean = 0
    valuesNum = 0
    for i in range(len(r_zVSz)):
        if r_zVSz[i][0] != 0:
            r_z_mean += r_zVSz[i][1]
            valuesNum += r_zVSz[i][0]
    print(r_z_mean / valuesNum)
    print(valuesNum)

rgyr()
