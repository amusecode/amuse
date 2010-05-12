from amuse.support.io import gadget
from amuse.support.io import nemobin
from amuse.support.units import nbody_system
from amuse.test import amusetest
from StringIO import StringIO

import os.path

class GadgetFileFormatProcessorTests(amusetest.TestCase):
    header_parts = ( 
        '\x00\x01\x00\x00 N\x00\x00 \xa1\x07\x00 N\x00' ,
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00' ,
        '\x00\x00\x00\x00\x00\x00\xf6`Q\xc6#\xcc\x9c>\x8d' ,
        '\xed\xb5\xa0\xf7\xc6\x90>\x00\x00\x00\x00\x00\x00\x00\x00' ,
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00' ,
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00' ,
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00 N\x00\x00 ' ,
        '\xa1\x07\x00 N\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00' ,
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00' ,
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00' ,
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00' ,
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00' ,
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00' ,
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00' ,
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00' ,
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00' ,
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00' ,
        '\x00\x00\x00\x00\x00\x00\x01\x00\x00' ,
    )

    def test1(self):
        header = ''.join(self.header_parts)
        x = gadget.GadgetFileFormatProcessor()
        file = StringIO(header)
        x.load_header(file)
        print x.header_struct
        self.assertEquals(x.header_struct.Npart[0], 20000)
        self.assertEquals(x.header_struct.Npart[1], 500000)
        self.assertEquals(x.header_struct.Npart[2], 20000)
        
        
        self.assertEquals(x.header_struct.Massarr[0], 0.0)
        self.assertAlmostRelativeEqual(x.header_struct.Massarr[1], 4.2911501e-07, 8)
        print x.header_struct.Massarr[2]
        self.assertAlmostRelativeEqual(x.header_struct.Massarr[2], 2.5000000e-07, 8)
        
        self.assertEquals(x.header_struct.FlagSfr, 0)
        self.assertEquals(x.header_struct.FlagFeedback, 0)
        self.assertEquals(x.header_struct.FlagAge, 0)
        self.assertEquals(x.header_struct.HubbleParam, 0)
        
    def xtest2(self):
        # turned of as it needs a large file, will make a 
        # small file for testing
        header = ''.join(self.header_parts)
        x = gadget.GadgetFileFormatProcessor()
        file = open('/data2/vanelteren/develop/python/wolk/tureluur/gadget_snapshot','rb')
        x.load_header(file)
        positions = x.read_fortran_block_float_vectors(file)
        file.close()
        print x.total_number_of_particles 
        print positions
        self.assertEquals(x.total_number_of_particles, len(positions))
        

class NemoBinaryFileFormatProcessorTests(amusetest.TestCase):
    
    
    def test1(self):
        directory_name = os.path.dirname(__file__)
        filename = os.path.join(directory_name, 'plummer128.nemo')
        file = open(filename, 'rb')
        nemofile = nemobin.NemoBinaryFile(file)
        tagcharacter, tagstring, dim, mustswap = nemofile.get_item_header()
        self.assertEquals(tagcharacter, 'c')
        self.assertEquals(tagstring, 'Headline')
        self.assertEquals(len(dim), 1)
        self.assertEquals(dim[0], 28)
        file.close()
        
    def test2(self):
        directory_name = os.path.dirname(__file__)
        filename = os.path.join(directory_name, 'plummer128.nemo')
        file = open(filename, 'rb')
        nemofile = nemobin.NemoBinaryFile(file)
        item = nemofile.read_item()
        
        self.assertEquals(item.data, "init_xrandom: seed used 123")
        file.close()
        
    
    def test3(self):
        directory_name = os.path.dirname(__file__)
        filename = os.path.join(directory_name, 'plummer128.nemo')
        file = open(filename, 'rb')
        nemofile = nemobin.NemoBinaryFile(file)
        data = nemofile.read()
        file.close()
        print data
        self.assertEquals(len(data), 3)
        tags = list(data.keys())
        self.assertEquals(tags[0] , 'Headline')
        self.assertEquals(tags[1] , 'History')
        self.assertEquals(tags[2] , 'SnapShot')
        self.assertEquals(data['History'][0].data, 'mkplummer out=plummer128.nemo nbody=128 seed=123 VERSION=2.8b')
        
        self.assertEquals(len(data['SnapShot'][0].data), 2)
        tags = list(data['SnapShot'][0].data.keys())
        self.assertEquals(tags[0] , 'Parameters')
        self.assertEquals(tags[1] , 'Particles')
        
    
    
    def test4(self):
        directory_name = os.path.dirname(__file__)
        filename = os.path.join(directory_name, 'plummer128.nemo')
        file = open(filename, 'rb')
        x = nemobin.NemoBinaryFileFormatProcessor()

        set = x.load_file(file)
        file.close()
        self.assertEquals(len(set), 128)
        self.assertAlmostRelativeEquals(set.kinetic_energy(), 0.230214395174 | nbody_system.energy, 8)
        self.assertAlmostRelativeEquals(set.potential_energy(G=nbody_system.G), -0.473503040144  | nbody_system.energy, 8)        

