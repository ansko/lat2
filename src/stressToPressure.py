#!/usr/bin/env python3
#coding utf-8

from options.options import options, universalOptions
class stressToPressure():
    def __init__(self, stresses, lx, ly, lz):
        self.uo = universalOptions()
        self.stresses = stresses
        self.lx = lx
        self.ly = ly
        self.lz = lz
        
    def presOverall(self):
        volume = self.lx * self.ly * self.lz
        self.overallPressure = 0
        for j in range(len(self.stresses)):
            self.overallPressure += self.stresses[j][2] / volume

    def pres1A(self):
        self.pressOneStep = [0 for j in range(self.uo.multip *
                                              self.uo.bigNumber)]
        surface = self.lx * self.ly
        for j in range(len(self.stresses)):
            z = int(self.uo.multip * self.stresses[j][1])
            self.pressOneStep[z] += self.stresses[j][2] / surface

    def pres7L(self, o):
        self.pressOneStep = [0 for j in range(7)]
        surface = self.lx * self.ly
        volume1 = (o.props['hi1'] - o.props['lo1']) * surface
        volume2 = (o.props['hi2'] - o.props['lo2']) * surface
        volume3 = (o.props['hi3'] - o.props['lo3']) * surface
        volume4 = (o.props['hi4'] - o.props['lo4']) * surface
        volume5 = (o.props['hi5'] - o.props['lo5']) * surface
        volume6 = (o.props['hi6'] - o.props['lo6']) * surface
        volume7 = (o.props['hi7'] - o.props['lo7']) * surface
        for j in range(len(self.stresses)):
            z = int(self.uo.multip * self.stresses[j][1])
            if o.props['lo1'] < z < o.props['hi1']:
                self.pressOneStep[0] += self.stresses[j][2] / volume1
            elif o.props['lo2'] < z < o.props['hi2']:
                self.pressOneStep[1] += self.stresses[j][2] / volume2
            elif o.props['lo3'] < z < o.props['hi3']:
                self.pressOneStep[2] += self.stresses[j][2] / volume3
            elif o.props['lo4'] < z < o.props['hi4']:
                self.pressOneStep[3] += self.stresses[j][2] / volume4
            elif o.props['lo5'] < z < o.props['hi5']:
                self.pressOneStep[4] += self.stresses[j][2] / volume5
            elif o.props['lo6'] < z < o.props['hi6']:
                self.pressOneStep[5] += self.stresses[j][2] / volume6
            elif o.props['lo7'] < z < o.props['hi7']:
                self.pressOneStep[6] += self.stresses[j][2] / volume7
