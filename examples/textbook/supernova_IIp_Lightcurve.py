## multi-component fits for various supernovae from
## https://arxiv.org/pdf/1404.2004.pdf
import numpy
from amuse.lab import *
#SN 2012au
#SN 2013ej
class Supernova_Ibc:
    def __init__(self, tstart=0|units.day):
        #SN2005hm
        self.mass = 6 | units.MSun
        self.A = 113.38 #muJy
        self.A = (6.7e+42) | units.erg/units.s
        self.C = 2.461  #muJy 
        self.t0 = tstart
        self.tfall = 2.54 | units.day
        self.trise = 5.61 | units.day
        
    def luminosity_at_time(self, time):
        if time > self.t0:
            F = self.A * numpy.exp(-(time-self.t0)/self.tfall)\
                         /  (1 + numpy.exp(-(time-self.t0)/self.trise))
        else:
            F = 6.4e-5 * self.A
        return F

class Supernova_IIp:
    def __init__(self, name, tstart=0|units.day):

        self.tstart = tstart
        self.t0 = 0
        self.particle = Particle()
        self.set_supernova_parameters(name)
        self.print_parameters()

    def print_parameters(self):
        print "Supernova parameters: "
        print "Name=", self.name
        print "Etot = ", self.Etot.in_(units.erg)
        print "Lpeak= ", self.Lpeak.in_(units.LSun)
        print "Mass parameters:", self.Mp, self.A, self.B1, self.B2, \
              self.BdN, self.BdC
        print "time-scales parameters:", self.t1, self.tp, self.t2, self.td

    def set_supernova_parameters(self, name):
        self.name = name
        if "12bku" in name:
            print "Supernova template: PS1-12bku" 
            self.mass = 20.0 | units.MSun
            self.Etot = 6.84e+51 | units.erg
            self.Mp = 2.46
            self.Lpeak= ((10**42.62)|units.erg/units.s)/self.Mp
            self.A= numpy.exp(-0.9)
            self.B1=numpy.exp(-0.5)
            self.B2=numpy.exp(-3.5)
            self.BdN=numpy.exp(-3.0)
            self.BdC=numpy.exp(-5.3)
            self.t1=1.0
            self.tp=5.0
            self.t2=105.0
            self.td=10.0
        elif "11aof" in name:
            print "Supernova template: PS1-11aof"
            self.mass = 23.5 | units.MSun
            self.Etot = 6.84e+51 | units.erg
            self.Mp = 2.05
            self.Lpeak= 10**43.06 |units.erg/units.s
            self.A= numpy.exp(-1)
            self.B1=numpy.exp(-6)
            self.B2=numpy.exp(-3.8)
            self.BdN=numpy.exp(-3.0)
            self.BdC=numpy.exp(-4.9)
            self.t1=1.
            self.tp=25.
            self.t2=105.
            self.td=10.
        elif "10a" in name:
            #for SN PS1-10a
            print "Supernova template: PS1-10a (default)"
            self.mass = 6 | units.MSun
            self.Etot = 6.84e+51 | units.erg
            self.Mp = 2.1
            self.Lpeak= 10**42.36 |units.erg/units.s
            self.A= numpy.exp(-1)
            self.B1=numpy.exp(-2.4)
            self.B2=numpy.exp(-3.1)
            self.BdN=numpy.exp(-2.9)
            self.BdC=numpy.exp(-5.0)
            self.t1=1.
            self.tp=5.
            self.t2=106.
            self.td=10.
        else:
            print "Supernova template: PS1-10a (default)"
            self.mass = 20 | units.MSun
            self.Etot = 6.84e+51 | units.erg
            self.Mp = 2.1
            self.Lpeak= 10**42.36 |units.erg/units.s
            self.A= numpy.exp(-1)
            self.B1=numpy.exp(-2.4)
            self.B2=numpy.exp(-3.1)
            self.BdN=numpy.exp(-2.9)
            self.BdC=numpy.exp(-5.0)
            self.t1=1.
            self.tp=5.
            self.t2=106.
            self.td=10.

    # def supernova_IIp_lightcurve(self, time):

    def luminosity_at_time(self, time):
        t = time.value_in(units.day) - self.tstart.value_in(units.day)
        L0 = 6.4e-5
        L = 0
        if t <= self.t0:
            L = L0
        elif t < self.t0+self.t1:
            L = L0 + self.M1() * (t/self.t1)**self.A
        elif t < self.t0+self.t1+self.tp:
            L = L0 + self.M1() * numpy.exp(self.B1*(t-self.t1))
        elif t < self.t0+self.t1+self.tp+self.t2:
            L = L0 + self.Mp * numpy.exp(-self.B2*(t-(self.t1+self.tp)))
        elif t < self.t0+self.t1+self.tp+self.t2+self.td:
            L = L0 + self.M2() \
                      * numpy.exp(-self.BdN*(t-(self.t1+self.t2+self.tp)))
        elif t>=self.t0+self.t1+self.tp+self.t2+self.td:
            L = self.Md() \
                 * numpy.exp(-self.BdC*(t-(self.t1+self.t2+self.tp+self.td)))
        return self.Lpeak * L 

    def M1(self):
        return self.Mp / numpy.exp(self.B1*self.tp)
    def M2(self):
        return self.Mp * numpy.exp(-self.B2*self.t2)
    def Md(self):
        return self.M2() * numpy.exp(-self.BdN*self.td)

if __name__ in ('__main__', '__plot__'):
    
    PS1_12bku = Supernova_IIp("12bku", 10|units.day)
    SN2005hm = Supernova_Ibc(10|units.day)
    PS1_11aof = Supernova_IIp("11aof", 10|units.day)
    PS1_10a  = Supernova_IIp("10a", 10|units.day)

    SN2005hm = Supernova_Ibc(10|units.day)
    
    t = 10**numpy.arange(-2, 2.5, 0.001) | units.day
    LIb = [] | units.erg/units.s
    L12bku = [] | units.erg/units.s
    L10a = [] | units.erg/units.s
    Laof = [] | units.erg/units.s
    for ti in t:
        LIb.append(SN2005hm.luminosity_at_time(ti))
        L10a.append(PS1_10a.luminosity_at_time(ti))
        L12bku.append(PS1_12bku.luminosity_at_time(ti))
        Laof.append(PS1_11aof.luminosity_at_time(ti))
    L10a = numpy.log10(L10a.value_in(units.LSun))
    L12bku = numpy.log10(L12bku.value_in(units.LSun))
    LIb = numpy.log10(LIb.value_in(units.LSun))
    Laof = numpy.log10(Laof.value_in(units.LSun))
    from matplotlib import pyplot
    pyplot.scatter(t.value_in(units.day), LIb)
    pyplot.plot(t.value_in(units.day), L10a, c='r')
    pyplot.plot(t.value_in(units.day), L12bku, c='b')
    pyplot.plot(t.value_in(units.day), Laof, c='g')
    pyplot.show()


