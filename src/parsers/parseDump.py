#!/usr/bin/env python3
#coding utf-8

class parseDump():
    def __init__(self, fname):
        self.parse(fname)

    def parse(self, fname):
        f = open(fname, 'r')
        f.readline()
        self.timeStep = int(f.readline())
        f.readline()
        self.atomsNumber = int(f.readline())
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        self.stresses = []
        for i in range(self.atomsNumber):
            ls = f.readline().split()
            self.stresses.append([int(ls[0]),
                                  float(ls[1]),
                                  float(ls[2])])
