import struct
import array
import numpy

from collections import namedtuple

from amuse.support.data import core
from amuse.support.core import late, OrderedMultiDictionary
from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.io import base




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
        return metaclass.mapping[typecharacter](tagstring, dimensions, mustswap = mustswap)
        
class NemoItem(object):
    __metaclass__ = NemoItemType
    
    def __init__(self, tagstring, dimensions = [1], data = None, mustswap = False):
        self.tagstring = tagstring
        self.dimensions = dimensions 
        self.mustswap = mustswap
        self.data = data
    
    def is_plural(self):
        if len(self.dimensions) == 1 and self.dimensions[0] <= 1:
            return False
        else:
            return True
            
    @late
    def number_of_values(self):
        return numpy.prod(self.dimensions)
        
    def read(self, nemofile):
        if self.datatype is None:
            pass
        else:
            self.data = nemofile.read_fixed_array(self.datatype, self.number_of_values)
        self.postprocess()
    
    def write(self, nemofile):
        if self.datatype is None:
            pass
        else:
            nemofile.write_fixed_array(self.preprocess(), self.datatype)
        
    def postprocess(self):
        if self.mustswap:
            self.data.byteswap()
    
    def preprocess(self):
        return self.data
    
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
        
    def preprocess(self):
        result = numpy.array(list(self.data), numpy.character)
        result = numpy.append(result, '\x00')
        return result
        
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
    
    def __init__(self, tagstring, dimensions = [1], data = None, mustswap = False):
        if data is None:
            data = OrderedMultiDictionary()
        NemoItem.__init__(self, tagstring, dimensions, data, mustswap)
    
    def read(self, nemofile):
        self.data = OrderedMultiDictionary()
        subitem = nemofile.read_item()
        while not subitem.isEndOfSet():
            self.data[subitem.tagstring] = subitem
            subitem = nemofile.read_item()
    
    def write(self, nemofile):
        for x in self.data.values():
            nemofile.write_item(x)
        nemofile.write_item(TesItem(self.tagstring, [1]))
        
    def add_item(self, item):
        self.data[item.tagstring] = item
        
    
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
            
    def write(self, nemofile):
        for x in self.data.values():
            nemofile.write_item(x)
        nemofile.write_item(YrotsItem(self.tagstring, [1]))
    
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
        bytes = self.file.read(datatype.itemsize * count)
        return numpy.fromstring(bytes, dtype=datatype,)


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
            raise base.IoException("Item does not have a valid header")
        
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
    
    def write_magic_number(self, is_plural):
        if is_plural:
            magic_number = self.PlurMagic
        else:
            magic_number = self.SingMagic
            
        x = array.array('h', [magic_number])
        self.file.write(x.tostring())
        
    
    def write_array(self, typetag, data):
        x = array.array(typetag, data)
        x.append(0)
        self.file.write(x.tostring())
        
    def write_string(self, string):
        return self.write_array('b', string)
        
    def write_item_header(self, item):
        self.write_magic_number(item.is_plural())
        self.write_string(item.typecharacter)
        if not item.typecharacter == TesItem.typecharacter:
            self.write_string(item.tagstring)
        if item.is_plural():
            self.write_array('i', item.dimensions)
        
    def write_item(self, item):
        self.write_item_header(item)
        item.write(self)
        
    def write_fixed_array(self, data, datatype):
        temp  = numpy.array(data, dtype=datatype)
        self.file.write(temp.tostring())
    
    
    def write(self, data):
        for x in data.values():
            self.write_item(x)
    
    
    
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
            nparticles = int(parameters['Nobj'][0].data[0])
            time = parameters['Time'][0].data[0]
            
            if result is None:
                result = core.Particles(nparticles)
                
            particlesitem = snapshot.data['Particles'][0]
            if 'PhaseSpace' in particlesitem.data:
                positions_and_velocities = particlesitem.data['PhaseSpace'][0].data
                positions = positions_and_velocities[...,0,...]
                velocities = positions_and_velocities[...,1,...]
            else:
                positions =  particlesitem.data['Position'][0].data
                positions = positions.reshape(positions.shape[:-1])
                velocities = particlesitem.data['Vecolity'][0].data
                velocities = velocities.reshape(velocities.shape[:-1])
            
            result.position = nbody_system.length.new_quantity(positions)
            result.velocity = nbody_system.speed.new_quantity(velocities)
            
            if 'Mass' in particlesitem.data:
                mass = particlesitem.data['Mass'][0].data
                result.mass = nbody_system.mass.new_quantity(mass)
                
                
            result.savepoint(time | nbody_system.time)
        return result
    
    def store_file(self, file):
        
        nemofile = NemoBinaryFile(file)
        item = SetItem('SnapShot')
        parameters_item = SetItem('Parameters')
        item.add_item(parameters_item)
        parameters_item.add_item(DoubleItem('Nobj', data = len(self.set)))
        if self.set.get_timestamp() is None:
            parameters_item.add_item(DoubleItem('Time', data = 0.0))
        else:
            parameters_item.add_item(DoubleItem('Time', data = set.get_timestamp().value_in(nbody_system.time)))
        
        particles_item = SetItem('Particles')
        item.add_item(particles_item)
        particles_item.add_item(
            DoubleItem(
                'Position', 
                dimensions=(len(self.set), 3, 1,), 
                data = self.set.position.value_in(nbody_system.length)
            )
        )
        particles_item.add_item(DoubleItem('Vecolity', dimensions=(len(self.set), 3, 1,), data = self.set.velocity.value_in(nbody_system.speed)))
        particles_item.add_item(DoubleItem('Mass', dimensions=(len(self.set), ), data = self.set.mass.value_in(nbody_system.mass)))
        nemofile.write({'a':item})
