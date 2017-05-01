#!/usr/bin/env python3
#coding utf-8

import pprint
pprint = pprint.PrettyPrinter(indent=1).pprint

from src.parsers.parseData import parseData
from src.parsers.parseDump import parseDump
from options.options import options, universalOptions
from src.stressToPressure import stressToPressure
from src.fourier.approFourier import appro


folderNum = 5
systemName = "segregated"


def mainOverall():
    step = 100000
    filesNum = int(2500000 / step) + 1
    fnameOptions = 'options/' + systemName + '.json'
    o = options(fnameOptions)
    uo = universalOptions()
    folder = o.props['folders'][folderNum]
    f = open('log', 'w')
    f.write(folder + '\n')
    f.write('Period Magnitude Error\n')
    lx = o.props['xhi'] - o.props['xlo']
    ly = o.props['yhi'] - o.props['ylo']
    lz = o.props['zhi'] - o.props['zlo']
    pressOneStep = [0 for j in range(uo.multip * uo.bigNumber)]
    press = []
    for i in range(filesNum):
        print(i)
        fname = folder + 'ALLstress.' + str(i * step)
        pd = parseDump(fname)
        stp = stressToPressure(pd.stresses, lx, ly, lz)
        stp.presOverall()
        press.append(stp.overallPressure)
    for k in range(15):
        period = len(press)
        pressTmp = []
        for step in range(len(press)):
            pressTmp.append(press[step])
        a = appro(pressTmp)
        a.error()
        if a.magnitude != 0:
            f.write(str(period) + ' ' +
                    str(a.magnitude) + ' +- '+
                    str(a.error) + '\n')
        for step in range(len(press) - 1):
            press[step] += press[step + 1]
            press[step] /= 2
        firstelement = press[0]
        endelement = press[-1]
        press = press[1:-2:2]
        press.insert(0, firstelement)
        press.append(endelement)
        k += 1

def main1A():
    step = 100000
    filesNum = int(2500000 / step) + 1
    fnameOptions = 'options/' + systemName + '.json'
    o = options(fnameOptions)
    uo = universalOptions()
    folder = o.props['folders'][folderNum]
    f = open('log', 'w')
    f.write(folder + '\n')
    f.write('LayerNum Period Magnitude Error\n')
    lx = o.props['xhi'] - o.props['xlo']
    ly = o.props['yhi'] - o.props['ylo']
    lz = o.props['zhi'] - o.props['zlo']
    pressOneStep = [0 for j in range(uo.multip * uo.bigNumber)]
    pressures = []
    for i in range(filesNum):
        print(i)
        fname = folder + 'ALLstress.' + str(i * step)
        pd = parseDump(fname)
        stp = stressToPressure(pd.stresses, lx, ly, lz)
        stp.pres1A()
        pressures.append(stp.pressOneStep)
    for k in range(15):
        period = len(pressures)
        for layerNum in range(uo.multip * uo.bigNumber):
            press = []
            for step in range(len(pressures)):
                press.append(pressures[step][layerNum])
            pressTmp = []
            for step in range(len(press)):
                pressTmp.append(press[step])
            a = appro(pressTmp)
            a.error()
            if a.magnitude != 0:
                f.write(str(layerNum) + ' ' +
                        str(period) + ' ' +
                        str(a.magnitude) + ' +- '+
                        str(a.error) + '\n')
        for step in range(len(pressures) - 1):
            for layer in range(len(pressures[step])):
                pressures[step][layer] += pressures[step + 1][layer]
                pressures[step][layer] /= 2
        firstelement = pressures[0]
        endelement = pressures[-1]
        pressures = pressures[1:-2:2]
        pressures.insert(0, firstelement)
        pressures.append(endelement)
        k += 1

def main7L():
    step = 100
    filesNum = int(2500000 / step) + 1
    fnameOptions = 'options/' + systemName + '.json'
    o = options(fnameOptions)
    uo = universalOptions()
    folder = o.props['folders'][folderNum]
    f = open('log', 'w')
    f.write(folder + '\n')
    f.write('LayerNum Period Magnitude Error\n')
    lx = o.props['xhi'] - o.props['xlo']
    ly = o.props['yhi'] - o.props['ylo']
    lz = o.props['zhi'] - o.props['zlo']
    pressOneStep = [0 for j in range(uo.multip * uo.bigNumber)]
    pressures = []
    for i in range(1, filesNum):
        print(i)
        fname = folder + 'ALLstress.' + str(i * step)
        pd = parseDump(fname)
        stp = stressToPressure(pd.stresses, lx, ly, lz)
        stp.pres7L(o)
        pressures.append(stp.pressOneStep)
    for k in range(15):
        period = len(pressures)
        for layerNum in range(7):
            press = []
            for step in range(len(pressures)):
                press.append(pressures[step][layerNum])
            pressTmp = []
            for step in range(len(press)):
                pressTmp.append(press[step])
            a = appro(pressTmp)
            a.error()
            if a.magnitude != 0:
                f.write(str(layerNum) + ' ' +
                        str(period) + ' ' +
                        str(a.magnitude) + ' +- '+
                        str(a.error) + '\n')
        for step in range(len(pressures) - 1):
            for layer in range(len(pressures[step])):
                pressures[step][layer] += pressures[step + 1][layer]
                pressures[step][layer] /= 2
        firstelement = pressures[0]
        endelement = pressures[-1]
        pressures = pressures[1:-2:2]
        pressures.insert(0, firstelement)
        pressures.append(endelement)
        k += 1


def main1A7L():
    step = 100
    filesNum = int(2500000 / step) + 1
    fnameOptions = 'options/' + systemName + '.json'
    o = options(fnameOptions)
    uo = universalOptions()
    folder = o.props['folders'][folderNum]
    f = open('log', 'w')
    f.write(folder + '\n')
    f.write('LayerNum Period Magnitude Error\n')
    lx = o.props['xhi'] - o.props['xlo']
    ly = o.props['yhi'] - o.props['ylo']
    lz = o.props['zhi'] - o.props['zlo']
    pressOneStep = [0 for j in range(uo.multip * uo.bigNumber)]
    pressures = []
    for i in range(filesNum):
        print(i)
        fname = folder + 'ALLstress.' + str(i * step)
        pd = parseDump(fname)
        stp = stressToPressure(pd.stresses, lx, ly, lz)
        stp.pres1A()
        pressOneStep7L = [0 for i in range(500)]
        for i in range(59, len(pressOneStep)):
            pressOneStep7L[int((i - 59) / 4.5)] += stp.pressOneStep[i]
        pressures.append(pressOneStep7L)
    for k in range(15):
        period = len(pressures)
        for layerNum in range(len(pressures[0])):
            press = []
            for step in range(len(pressures)):
                press.append(pressures[step][layerNum] / 4.5)
            pressTmp = []
            for step in range(len(press)):
                pressTmp.append(press[step])
            a = appro(pressTmp)
            a.error()
            if a.magnitude != 0:
                f.write(str(layerNum) + ' ' +
                        str(period) + ' ' +
                        str(a.magnitude) + ' +- '+
                        str(a.error) + '\n')
        for step in range(len(pressures) - 1):
            for layer in range(len(pressures[step])):
                pressures[step][layer] += pressures[step + 1][layer]
                pressures[step][layer] /= 2
        firstelement = pressures[0]
        endelement = pressures[-1]
        pressures = pressures[1:-2:2]
        pressures.insert(0, firstelement)
        pressures.append(endelement)
        k += 1
    

#mainOverall()
#main1A()
#main7L()
main1A7L()
