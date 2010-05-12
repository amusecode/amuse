from amuse.support.data import core
from amuse.support.core import late, OrderedMultiDictionary
from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.io import base

from collections import namedtuple


import struct
import array
import numpy


class NemoItemType(type):
    mapping = {}
    def __new__(metaclass, name, bases, dict):
        if 'datatype' in dict:
            if not dict['datatype'] is None:
                dict['datatype'] = numpy.dtype((dict['datatype'], 1,))
        result =  type.__new__(metaclass, name, bases, dict)
        if 'typecharacter' in dict:
            metaclass.mapping[dict['typecharacter']] = result
        return result
    
    @classmethod
    def new_item(metaclass, typecharacter, tagstring, dimensions, mustswap = False):
        return metaclass.mapping[typecharacter](tagstring, dimensions, mustswap)
        
class NemoItem(object):
    __metaclass__ = NemoItemType
    
    def __init__(self, tagstring, dimensions, mustswap = False):
        self.tagstring = tagstring
        self.dimensions = dimensions 
        self.mustswap = mustswap
        self.data = None
        
    @late
    def number_of_values(self):
        return numpy.prod(self.dimensions)
        
    def read(self, nemofile):
        if self.datatype is None:
            pass
        else:
            self.data = nemofile.read_fixed_array(self.datatype, self.number_of_values)
        self.postprocess()
    
    def postprocess(self):
        if self.mustswap:
            self.data.byteswap()
    
    
    def isEndOfSet(self):
        return False
        
    def isEndOfHistory(self):
        return False
        
    def __str__(self):
        return 'nemoitem({0},{1})'.format(self.tagstring, self.dimensions)
        
    def __repr__(self):
        return '<{0!s} {1},{2}>'.format(type(self), self.tagstring, self.dimensions)
        
class AnyItem(NemoItem):
    """anything at all"""
    typecharacter = "a"
    datatype = numpy.byte
    
class CharItem(NemoItem):
    """printable chars"""
    typecharacter = "c"
    datatype = numpy.character
    
    def postprocess(self):
        self.data = ''.join(self.data[:-1])
        
class ByteItem(NemoItem):
    """unprintable chars"""
    typecharacter = "b"
    datatype = numpy.character
    
    
    def postprocess(self):
        self.data = self.data.reshape(self.dimensions)
    
class ShortItem(NemoItem):
    """  short integers """
    typecharacter = "s"
    datatype = numpy.int16
    
    def postprocess(self):
        self.data = self.data.reshape(self.dimensions)
    
class IntItem(NemoItem):
    """  standard integers """
    typecharacter = "i"
    datatype = numpy.int32
    
    def postprocess(self):
        self.data = self.data.reshape(self.dimensions)
    
class LongItem(NemoItem):
    """  long integers """
    typecharacter = "l"
    datatype = numpy.int64
    
    def postprocess(self):
        self.data = self.data.reshape(self.dimensions)
    
class HalfpItem(NemoItem):
    """  half precision floating """
    typecharacter = "h"
    datatype = numpy.int16
    
class FloatItem(NemoItem):
    """  short floating """
    typecharacter = "f"
    datatype = numpy.float32
    
    def postprocess(self):
        self.data = self.data.reshape(self.dimensions)
    
class DoubleItem(NemoItem):
    """  long floating """
    typecharacter = "d"
    datatype = numpy.float64
    
    def postprocess(self):
        self.data = self.data.reshape(self.dimensions)
    
class SetItem(NemoItem):
    """  begin compound item """
    typecharacter = "("
    datatype = None
    
    
    def read(self, nemofile):
        self.data = OrderedMultiDictionary()
        subitem = nemofile.read_item()
        while not subitem.isEndOfSet():
            self.data[subitem.tagstring] = subitem
            subitem = nemofile.read_item()
    
class TesItem(NemoItem):
    """  end of compound item """
    typecharacter = ")"
    datatype = None
    
    def isEndOfSet(self):
        return True
        
    def read(self, file):
        pass
    
#  Experimental Kludge 
class StoryItem(NemoItem):
    """  begin of a story item (see starlab) """
    typecharacter = "["
    typecharacter = "["
    datatype = None
    
    def read(self, nemofile):
        self.data = OrderedMultiDictionary()
        subitem = nemofile.read_item()
        while not subitem.isEndOfHistory():
            self.data[subitem.tagstring] = subitem
            subitem = nemofile.read_item()
    
class YrotsItem(NemoItem):
    """  end of a story item (see starlab) """
    typecharacter = "]"
    datatype = None
        
    def isEndOfHistory(self):
        return True

 
        
class NemoBinaryFile(object):
    
    def __init__(self, file):
        self.file = file
    
    SingMagic  = ((011<<8) + 0222)
    PlurMagic  = ((013<<8) + 0222)

    def _byteswap(self, value, type = 'H'):
        x = array.array('H', [value])
        x.byteswap()
        return x[0]
        
    @late
    def reversed_SingMagic(self):
        return self._byteswap(self.SingMagic)
        
    @late
    def reversed_PlurMagic(self):
        return self._byteswap(self.PlurMagic)
    
    def read_magic_number(self):
        nbytes = 2
        bytes = self.file.read(nbytes)
        if not bytes or len(bytes) < nbytes:
            return None
        return array.array('h', bytes)[0]
        
    def read_array(self, typetag):
        result = array.array(typetag)
        counter = result.itemsize
        must_loop = True
        while must_loop:
            block = self.file.read(result.itemsize)
            result.fromstring(block)
            must_loop = result[-1] != 0
            
        result.pop()
        return result
        
    def read_string(self):
        return self.read_array('b').tostring()

    def read_fixed_array(self, datatype, count):
        return numpy.fromfile(self.file, dtype=datatype, count=count)


    def get_item_header(self):
        
        magic_number = self.read_magic_number()
        if magic_number is None:
            return (None, None, None, None)
        mustswap = False
        if magic_number == self.reversed_SingMagic:
            magic_number == self.SingMagic
            mustswap = True
        elif magic_number == self.reversed_PlurMagic:
            magic_number == self.PlurMagic
            mustswap = True
            
        if not (magic_number == self.SingMagic or magic_number == self.PlurMagic):
            raise Exception("Item does not have a valid header")
        
        typecharacter = self.read_string()
        if not typecharacter == TesItem.typecharacter:
            tagstring = self.read_string()
        else:
            tagstring = ''
        if magic_number ==  self.PlurMagic:
            dim = self.read_array('i')
            if mustswap:
                dim.byteswap()
            dim = dim.tolist()
        else:
            dim = [1]
            
        return (typecharacter, tagstring, dim, mustswap)
        
    def read_item(self):
        typecharacter, tagstring, dim, mustswap = self.get_item_header()
        if typecharacter is None:
            return None
        result = NemoItemType.new_item(typecharacter, tagstring, dim, mustswap)
        result.read(self)
        return result
    
    def read(self):
        result = OrderedMultiDictionary()
        
        item = self.read_item()
        while not item is None:
            result[item.tagstring] = item
            item = self.read_item()
        
        return result
    
    
class NemoBinaryFileFormatProcessor(base.BinaryFileFormatProcessor):
    provided_formats = ['nemobin']
    
    
        
    def load_file(self, file):
        nemofile = NemoBinaryFile(file)
        data = nemofile.read()
        result = None
        for snapshot in data['SnapShot']:
            if not 'Particles' in snapshot.data:
                continue
            
            parameters = snapshot.data['Parameters'][0].data
            nparticles = parameters['Nobj'][0].data[0]
            time = parameters['Time'][0].data[0]
            
            if result is None:
                result = core.Particles(nparticles)
                
            particlesitem = snapshot.data['Particles'][0]
            print particlesitem.data.keys()
            if 'PhaseSpace' in particlesitem.data:
                positions_and_velocities = particlesitem.data['PhaseSpace'][0].data
                positions = positions_and_velocities[...,0,...]
                velocities = positions_and_velocities[...,1,...]
            else:
                positions =  particlesitem.data['Position'][0].data
                velocities = particlesitem.data['Vecolity'][0].data
            
            result.position = nbody_system.length.new_quantity(positions)
            result.velocity = nbody_system.speed.new_quantity(velocities)
            
            if 'Mass' in particlesitem.data:
                result.mass = nbody_system.mass.new_quantity(particlesitem.data['Mass'][0].data)
                
                
            result.savepoint(time | nbody_system.time)
        return result