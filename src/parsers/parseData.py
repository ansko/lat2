#!/usr/bin/env python3
#coding utf-8

class parseData():
    def __init__(self, fname):
        self.parse(fname)

    def parse(self, fname):
        f = open(fname, 'r')
        header = f.readline()
        f.readline()
        self.atomsNumber = int(f.readline().split()[0])
        self.atomsTypes = int(f.readline().split()[0])
        self.bondsNumber = int(f.readline().split()[0])
        self.bondsTypes = int(f.readline().split()[0])
        self.anglesNumber = int(f.readline().split()[0])
        self.anglesTypes = int(f.readline().split()[0])
        self.dihedralsNumber = int(f.readline().split()[0])
        self.dihedralsTypes = int(f.readline().split()[0])
        self.impropersNumber = int(f.readline().split()[0])
        self.impropersTypes = int(f.readline().split()[0])
        f.readline()
        ls = f.readline().split()
        (self.xlo, self.xhi) = (float(ls[0]), float(ls[1]))
        ls = f.readline().split()
        (self.ylo, self.yhi) = (float(ls[0]), float(ls[1]))
        ls = f.readline().split()
        (self.zlo, self.zhi) = (float(ls[0]), float(ls[1]))
        a = f.readline()
        if len(a) > 1:
            f.readline()
        f.readline()
        f.readline()
        self.masses = []
        for i in range(self.atomsTypes):
            self.masses.append(float(f.readline().split()[1]))
        f.readline()
        f.readline()
        f.readline()
        self.pairCoeffs = []
        for i in range(self.atomsTypes):
            ls = f.readline().split()
            self.pairCoeffs.append([float(ls[0]), float(ls[1])])
        f.readline()
        f.readline()
        f.readline()
        self.bondCoeffs = []
        for i in range(self.bondsTypes):
            ls = f.readline().split()
            self.bondCoeffs.append([float(ls[0]), float(ls[1])])
        f.readline()
        f.readline()
        f.readline()
        self.angleCoeffs = []
        for i in range(self.anglesTypes):
            ls = f.readline().split()
            self.angleCoeffs.append([float(ls[0]), float(ls[1])])
        f.readline()
        f.readline()
        f.readline()
        self.dihedralCoeffs = []
        for i in range(self.dihedralsTypes):
            ls = f.readline().split()
            self.dihedralCoeffs.append([float(ls[0]),
                                        float(ls[1]),
                                        float(ls[2])])
        f.readline()
        f.readline()
        f.readline()
        self.improperCoeffs = []
        for i in range(self.impropersTypes):
            ls = f.readline().split()
            self.improperCoeffs.append([float(ls[0]),
                                        float(ls[1]),
                                        float(ls[2])])
        f.readline()
        f.readline()
        f.readline()
        self.atomsFull = []
        for i in range(self.atomsNumber):
            ls = f.readline().split()
            self.atomsFull.append([int(ls[0]),
                                   int(ls[1]),
                                   int(ls[2]),
                                   float(ls[3]),
                                   float(ls[4]),
                                   float(ls[5]),
                                   float(ls[6]),
                                   int(ls[7]),
                                   int(ls[8]),
                                   int(ls[9])])
        self.atomsFull.sort()
        f.readline()
        bondsOrVelocities = f.readline()
        if bondsOrVelocities.startswith('Velocities'):
            for i in range(self.atomsNumber + 3):
                f.readline()
        f.readline()
        self.bonds = []
        for i in range(self.bondsNumber):
            ls = f.readline().split()
            self.bonds.append([int(ls[0]), 
                               int(ls[1]),
                               int(ls[2]),
                               int(ls[3])])
        self.bonds.sort()
        f.readline()
        f.readline()
        f.readline()
        self.angles = []
        for i in range(self.anglesNumber):
            ls = f.readline().split()
            self.angles.append([int(ls[0]),
                                int(ls[1]),
                                int(ls[2]),
                                int(ls[3]),
                                int(ls[4])])
        self.angles.sort()
        f.readline()
        f.readline()
        f.readline()
        self.dihedrals = []
        for i in range(self.dihedralsNumber):
            ls = f.readline().split()
            self.dihedrals.append([int(ls[0]),
                                   int(ls[1]),
                                   int(ls[2]),
                                   int(ls[3]),
                                   int(ls[4]),
                                   int(ls[5])])
        self.dihedrals.sort()
        f.readline()
        f.readline()
        f.readline()
        self.impropers = []
        for i in range(self.impropersNumber):
            ls = f.readline().split()
            self.impropers.append([int(ls[0]),
                                   int(ls[1]),
                                   int(ls[2]),
                                   int(ls[3]),
                                   int(ls[4]),
                                   int(ls[5])])
        self.impropers.sort()
