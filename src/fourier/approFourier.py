#!/usr/bin/env python3
#coding utf-8

import math

class appro():
    def __init__(self, pressures):
        self.pressures = pressures
        self.approximate()

    def approximate(self):
        pressures = self.pressures
        period = len(pressures)
        i = a0 = a1 = a2 = b1 = b2 = 0
        suma0 = suma1 = suma2 = sumb1 = sumb2 = magnitude = 0
        for i in range(period):
            suma0 += float(pressures[i])
            suma1 += (float(pressures[i]) * 
                            math.cos(2 * math.pi / period * (i + 1)))
            sumb1 += (float(pressures[i]) * 
                            math.sin(2 * math.pi / period * (i + 1)))
            suma2 += (float(pressures[i]) *
                            math.cos(4 * math.pi / period * (i + 1)))
            sumb2 += (float(pressures[i]) *
                            math.sin(4 * math.pi / period * (i + 1)))
        self.a0 = suma0 / period
        self.a1 = 2 * suma1 / period
        self.a2 = 2 * suma2 / period
        self.b1 = 2 * sumb1 / period
        self.b2 = 2 * sumb2 / period
        self.magnitude = 5 * math.sqrt(self.a1**2 + self.b1**2)

    def error(self):
        error = 0
        pressures = self.pressures
        period = len(pressures)
        for i in range(period):
            phase = 2 * math.pi / period * (i + 1)
            error += abs((self.a0 +
                          self.a1 * math.cos(phase) +
                          self.b1 * math.sin(phase)) - pressures[i])**2
        self.error = math.sqrt(error/len(pressures))
