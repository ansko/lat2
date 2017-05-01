#!/usr/bin/env python3
#coding utf-8

import copy
import math

import pprint
pprint = pprint.PrettyPrinter(indent=1).pprint

from src.parsers.parseData import parseData
from options.options import options, universalOptions


def bondsOrientation(systemName="mixed"):
    multip = 1
    fnameOptions = 'options/' + systemName + '.json'
    o = options(fnameOptions)
    uo = universalOptions()
    orientations = [[0, 0] for i in range(1500)]
    end = 51
    if systemName == '5x20':
        end = 50
    if systemName == '10x20':
        end = 51
    if systemName == 'long':
        folder = ('/media/anton/Seagate Expansion Drive/Article-MMT/Cluster' + 
                  ' calculations for article/BiggerSystems/Comp/1chain (long)/' +
                  '1845502 - wiggle1/')
    elif systemName == 'mixed':
        folder = ('/media/anton/Seagate Expansion Drive/Article-MMT/Cluster' + 
                  ' calculations for article/1st conf/mixed ok/300K/' + 
                  'npt after 6th cycle/1440224/')
    for i in range(1, 51):
        print(i)
        pd = parseData(folder + 'co.' + str(i * 50000) + '.data')
        lx = pd.xhi - pd.xlo
        ly = pd.yhi - pd.ylo
        lz = pd.zhi - pd.zlo
        atomsFull = pd.atomsFull
        bonds = pd.bonds
        for bond in bonds:
            if systemName == 'mixed':
                if bond[1] not in [3, 4, 5, 7, 8, 11, 12, 13, 14, 16]:
                #if bond[1] not in [5, 8]:
                    continue
            elif systemName == 'segregated':
                if bond[1] not in [3, 4, 5, 7, 8, 11, 12, 13, 14, 16]:
                    continue
            elif systemName == 'long':
                if bond[1] not in [4, 5, 6, 7, 8, 9, 11, 12, 14, 15]:
                    continue
            elif systemName == '5x20':
                if bond[1] not in [3, 4, 5, 7, 8, 11, 12, 13, 14, 16]:
                    continue
            elif systemName == '10x20':
                if bond[1] not in [3, 4, 5, 7, 8, 11, 12, 13, 14, 16]:
                    continue
            atomOne = atomsFull[bond[2] - 1]
            atomTwo = atomsFull[bond[3] - 1]
            bondx = min(abs(atomOne[4] - atomTwo[4]),
                        abs(lx - abs(atomOne[4] - atomTwo[4])))
            bondy = min(abs(atomOne[5] - atomTwo[5]),
                        abs(ly - abs(atomOne[5] - atomTwo[5])))
            bondz = min(abs(atomOne[6] - atomTwo[6]),
                        abs(lz - abs(atomOne[6] - atomTwo[6])))
            if bondx >= 3 or bondy >= 3 or bondz >= 3:
                print("Long bond ", bondx, bondy, bondz)
            cosTheta = bondz / ((bondx**2 + bondy**2 + bondz**2)**0.5)
            parameter = (3 * cosTheta**2 - 1) / 2
            z = int((250 + atomOne[6]) * multip)
            orientations[z][0] += 1
            orientations[z][1] += parameter
    f = open('log', 'w')
    for i in range(len(orientations)):
        if orientations[i][0] != 0:
            f.write(str(i / multip - 250) + ' ' +
                    str(orientations[i][1] / orientations[i][0]) + ' ' +
                    str(orientations[i][0]) + '\n')


def anglesOrientation(systemName="10x20"):
    multip = 1
    fnameOptions = 'options/' + systemName + '.json'
    o = options(fnameOptions)
    uo = universalOptions()
    folder = o.props['folders'][0][:-6:1]
    orientations = [[0, 0] for i in range(500)]
    end = 51
    if systemName == '5x20':
        end = 50
    if systemName == '10x20':
        end = 51
    for i in range(1, end):
        pd = parseData(folder + 'co.' + str(i * 50000) + '.data')
        lx = pd.xhi - pd.xlo
        ly = pd.yhi - pd.ylo
        lz = pd.zhi - pd.zlo
        atomsFull = pd.atomsFull
        angles = pd.angles
        for angle in angles:
            if systemName == 'mixed':
                if angle[1] not in [1, 2, 3, 7, 9, 10, 16, 19, 20, 21, 22, 25, 26,
                                   28, 29]:
                    continue
            elif systemName == 'segregated':
                if angle[1] not in [1, 2, 3, 7, 9, 10, 16, 19, 20, 21, 22, 25, 26,
                                   28, 29]:
                    continue
            elif systemName == 'long':
                if angle[1] not in [1, 2, 5, 7, 8, 9, 10, 13, 14, 16, 20, 22, 23,
                                    24, 26]:
                    continue
            elif systemName == '5x20':
                if angle[1] not in [1, 2, 3, 7, 9, 10, 16, 19, 20, 21, 22, 25, 26,
                                    28, 29]:
                    continue
            elif systemName == '10x20':
                if angle[1] not in [1, 2, 3, 7, 9, 10, 16, 19, 20, 21, 22, 25, 26,
                                    28, 29]:
                    continue
            atomOne = atomsFull[angle[2] - 1]
            atomMid = atomsFull[angle[3] - 1]
            atomTwo = atomsFull[angle[4] - 1]
            bondx = min(abs(atomOne[4] - atomTwo[4]),
                        lx - abs(atomOne[4] - atomTwo[4]))
            bondy = min(abs(atomOne[5] - atomTwo[5]),
                        lx - abs(atomOne[5] - atomTwo[5]))
            bondz = min(abs(atomOne[6] - atomTwo[6]),
                        lx - abs(atomOne[6] - atomTwo[6]))
            cosTheta = bondz / (bondx**2 + bondy**2 + bondz**2)**0.5
            parameter = (3 * cosTheta**2 - 1) / 2
            z = int((100 + atomMid[6]) * multip)
            orientations[z][0] += 1
            orientations[z][1] += parameter
    f = open('log', 'w')
    for i in range(len(orientations)):
        if orientations[i][0] != 0:
            f.write(str(i / multip - 100) + ' ' +
                    str(orientations[i][1] / orientations[i][0]) + ' ' +
                    str(orientations[i][0]) + '\n')


        
bondsOrientation()
#anglesOrientation()
