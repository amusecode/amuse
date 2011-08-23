import struct
import numpy

from collections import namedtuple

from amuse.support.data import core
from amuse.support.core import late
from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.io import base


ioversion=2

nheader=8
nihead=32
nrhead=32
nphead=64

# name - type - index in phead - target particle type - unit
variables=[
  [ "mass", 'd', 1, "all", nbody_system.mass],
  [ "x", 'd', 2, "all", nbody_system.length],
  [ "y", 'd', 2, "all", nbody_system.length],
  [ "z",'d',2,"all", nbody_system.length],
  [ "vx",'d',3,"all", nbody_system.speed],
  [ "vy",'d',3,"all", nbody_system.speed],
  [ "vz",'d',3,"all", nbody_system.speed],
  [ "radius",'d',4,"all", nbody_system.length],
  [ "tform",'d',5,"all", nbody_system.time],
  [ "accx",'d',6,"all", nbody_system.acceleration],
  [ "accy",'d',6,"all", nbody_system.acceleration],
  [ "accz",'d',6,"all", nbody_system.acceleration],
  [ "rho",'d',10,"gas", nbody_system.mass / nbody_system.length ** 3],
  [ "ethermal",'d',11,"gas", nbody_system.speed**2],
  [ "entropy",'d',12,"gas", units.none],   # not sure about this unit
  [ "h_smooth",'d',13,"gas", nbody_system.length],
  [ "fuvheat",'d',14,"gas", units.none],  # special case -> units of Habing field 
  [ "esnthdt",'d',15,"gas", nbody_system.speed**2/nbody_system.time],
  [ "tcollaps",'d',16,"gas", nbody_system.time],
  [ "temperat",'d',17,"gas", units.K],
  [ "elecfrac",'d',18,"gas", units.none],
  [ "csound",'d',19,"gas", nbody_system.speed],
  [ "pressure",'d',20,"gas", nbody_system.mass/nbody_system.length/nbody_system.time**2],
  [ "hsmdivv",'d',21,"gas", nbody_system.speed],
  [ "mumaxdvh",'d',22,"gas", nbody_system.speed],
  [ "hsmcurlv",'d',23,"gas", nbody_system.speed],
  [ "vdisp",'d',24,"gas", nbody_system.speed],
  [ "h2frac",'d',25,"gas", units.none],
  [ "dethdt",'d',26,"gas", nbody_system.speed**2/nbody_system.time],
  [ "dentdt",'d',27,"gas", units.none],   # not sure about this unit
  [ "starfuv",'d',34,"stars", units.erg/units.s],    # special case
  [ "snentropy",'d',35,"stars", units.none],  # not sure about this unit
  [ "pot", 'd', 40, "all", nbody_system.speed**2],
  [ "extpot", 'd', 41, "all", nbody_system.speed**2],
  [ "id", 'i', 42, "all", units.none],
  [ "timestep", 'i', 43, "all", units.none] ]

class FiFileFormatProcessor(base.FortranFileFormatProcessor):

    def load_header(self,file):
        self.ioversion = self.read_fortran_block_ints(file)[0]
        if self.ioversion != ioversion:
            raise Exception
        self.header=self.read_fortran_block_ints(file) 
        self.ihead=self.read_fortran_block_ints(file)
        self.rhead=self.read_fortran_block_doubles(file)
        self.phead=self.read_fortran_block_ints(file)
        
        self.ngas=self.ihead[1]
        self.ndm=self.ihead[2]
        self.nstar=self.ihead[3]
        
        self.unitm_in_msun= self.rhead[16]
        self.unitl_in_kpc= self.rhead[17]/3.086e21 
        
        if self.unitm_in_msun==0. or self.unitl_in_kpc==0.:
            self.unitm_in_msun= 1.e9
            self.unitl_in_kpc= 1.0 
               
                  
        self.convert=nbody_system.nbody_to_si( self.unitm_in_msun | units.MSun, self.unitl_in_kpc | units.kpc)
        
        
        print self.unitl_in_kpc,numpy.log10(self.unitl_in_kpc)
        self.flxscale=-41.282-2*numpy.log10(self.unitl_in_kpc)
        self.heatconst=2.0693e14*self.unitl_in_kpc**2.5/self.unitm_in_msun**1.5
        
       
    def load_body(self,file):
        self.gas=core.Particles(self.ngas)
        self.dark=core.Particles(self.ndm)
        self.star=core.Particles(self.nstar)
        
        for name, type, index, part, unit in variables:
            if self.phead[index - 1] == 1:
                if type == 'd':
                    tmp = self.read_fortran_block_doubles(file)
                if type == 'i':        
                    tmp = self.read_fortran_block_ints(file)
# special cases:
            if index == 14:
                tmp = tmp / self.heatconst              
            if index == 34:
                tmp = tmp / self.heatconst / 10 ** self.flxscale


            if part == "all":
                if self.ngas > 0:
                    setattr(self.gas, name, self.convert.to_si(unit.new_quantity(tmp[0:self.ngas])))
                if self.ndm > 0:
                    setattr(self.dark, name, self.convert.to_si(unit.new_quantity(tmp[self.ngas:self.ngas + self.ndm])))
                if self.nstar > 0:
                    setattr(self.star, name,
                            self.convert.to_si(unit.new_quantity(tmp[self.ngas + self.ndm:self.ngas + self.ndm + self.nstar])))              
            if part == "gas":
                if self.ngas > 0:
                    setattr(self.gas, name, self.convert.to_si(unit.new_quantity(tmp[0:self.ngas])))
            if part == "stars":    
                if self.nstar > 0:
                    setattr(self.star, name,
                            self.convert.to_si(unit.new_quantity(tmp[self.ngas + self.ndm:self.ngas + self.ndm + self.nstar])))              
    

    def load_file(self, file):
        self.load_header(file)
        self.load_body(file)
        return self.gas,self.dark,self.star
        
if __name__=="__main__":
    p=base.read_set_from_file("test2",format=FiFileFormatProcessor)  
    print p[0]
    print len(p[0]),len(p[1]),len(p[2])

