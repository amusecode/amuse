#!/usr/bin/python

from __future__ import division
from matplotlib import pyplot
import os

def logLuminosity(V):
    VminMv = 9.7
    MBolSun = 4.74
    Mv = V - VminMv
    L = (MBolSun - Mv)/2.5
    return L

def logTeff(BminV):
    logT = (14.551 - BminV)/3.684
    if logT > 3.961:
        a,b,c = [.344,-3.402,8.037]
        logT = (-b - (b*b-4*a*c)**.5)/(2.*a)
    return logT 

def Teff(BminV):
    return 10.**logTeff(BminV)

class Cluster():
    def __init__(self) :
        self.n = 0
        self.L = []
        self.Teff = []
        self.BmV = []
        self.V = []

    def __iter__(self):
        return self

    def __repr__(self):
        tt = 'Cluster()' 
        return tt

    def read(self):
        try:
            amusedir = os.environ['AMUSE_DIR']
            dir = amusedir+'/examples/textbook/'
        except:
            print 'Environment variable AMUSE_DIR not set'
            dir = './'
    
        isofile = open(dir+'M67Data.dat')
        lines = isofile.readlines()
        E_BminV = 0.04
        for line in lines:
            self.V.append(float(line.split()[0]))
            self.BmV.append(float(line.split()[1]))
            self.L.append(10**logLuminosity(self.V[-1]-E_BminV))
            self.Teff.append(10**logTeff(self.BmV[-1]-E_BminV))
            self.n = self.n +1

    def plot(self):
        pyplot.xlim(max(self.Teff), min(self.Teff))
        pyplot.scatter(self.Teff, self.L)
        pyplot.xlabel("Teff")
        pyplot.ylabel("L")
        pyplot.show()

if __name__=="__main__":
    cls = Cluster()
    cls.read()
    cls.plot()

