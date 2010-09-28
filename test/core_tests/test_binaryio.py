from amuse.support.io import gadget
from amuse.support.io import nemobin
from amuse.support import io
from amuse.support.units import nbody_system
from amuse.test import amusetest
from StringIO import StringIO

import os.path

class GadgetFileFormatProcessorTests(amusetest.TestCase):
    header_parts = ( 
        '\x00\x01\x00\x00 N\x00\x00 \xa1\x07\x00 N\x00',
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00',
        '\x00\x00\x00\x00\x00\x00\xf6`Q\xc6#\xcc\x9c>\x8d',
        '\xed\xb5\xa0\xf7\xc6\x90>\x00\x00\x00\x00\x00\x00\x00\x00',
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00',
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00',
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00 N\x00\x00 ',
        '\xa1\x07\x00 N\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00',
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00',
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00',
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00',
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00',
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00',
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00',
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00',
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00',
        '\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00',
        '\x00\x00\x00\x00\x00\x00\x01\x00\x00',
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
        
    def test2(self):
        directory_name = os.path.dirname(__file__)
        filename = os.path.join(directory_name, 'gadget_snapshot')
        x = gadget.GadgetFileFormatProcessor()
        file = open(filename,'rb')
        result = x.load_file(file)
        file.close()
        self.assertEquals(len(result[0]), 1000)
        self.assertEquals(len(result[1]), 10000)
        
    def test3(self):
        directory_name = os.path.dirname(__file__)
        filename = os.path.join(directory_name, 'gadget_snapshot')
        x = gadget.GadgetFileFormatProcessor()
        result = io.read_set_from_file(filename, format='gadget')
        self.assertEquals(len(result[0]), 1000)
        self.assertEquals(len(result[1]), 10000)
        
    def test4(self):
        directory_name = os.path.dirname(__file__)
        filename = os.path.join(directory_name, 'gassphere_littleendian.dat')
        x = gadget.GadgetFileFormatProcessor()
        result = io.read_set_from_file(filename, format='gadget')
        self.assertEquals(len(result[0]), 1472)
        self.assertEquals(len(result[1]), 0)
        

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
        self.assertEquals(len(data), 3)
        tags = list(data.keys())
        self.assertEquals(tags[0], 'Headline')
        self.assertEquals(tags[1], 'History')
        self.assertEquals(tags[2], 'SnapShot')
        self.assertEquals(data['History'][0].data, 'mkplummer out=plummer128.nemo nbody=128 seed=123 VERSION=2.8b')
        
        self.assertEquals(len(data['SnapShot'][0].data), 2)
        tags = list(data['SnapShot'][0].data.keys())
        self.assertEquals(tags[0], 'Parameters')
        self.assertEquals(tags[1], 'Particles')
        
    
    
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

    def test5(self):
        directory_name = os.path.dirname(__file__)
        filename = os.path.join(directory_name, 'plummer128.nemo')
        file = open(filename, 'rb')
        nemofile = nemobin.NemoBinaryFile(file)
        data = nemofile.read()
        file.close()
        
        outputfile =  StringIO()
        nemooutputfile = nemobin.NemoBinaryFile(outputfile)
        nemooutputfile.write(data)
        string = outputfile.getvalue()
        outputfile.close()
        inputfile =  StringIO(string)
        nemoinputfile = nemobin.NemoBinaryFile(inputfile)
        tagcharacter, tagstring, dim, mustswap = nemoinputfile.get_item_header()
        self.assertEquals(tagcharacter, 'c')
        self.assertEquals(tagstring, 'Headline')
        self.assertEquals(len(dim), 1)
        self.assertEquals(dim[0], 28)
        inputfile.close()
        
    
    def test6(self):        
        inputfile =  StringIO()
        nemoinputfile = nemobin.NemoBinaryFile(inputfile)
        data = nemoinputfile.read()
        self.assertEquals(len(data), 0)
        
    
    def test7(self):
        directory_name = os.path.dirname(__file__)
        filename = os.path.join(directory_name, 'plummer128.nemo')
        file = open(filename, 'rb')
        nemofile = nemobin.NemoBinaryFile(file)
        data = nemofile.read()
        file.close()
        
        outputfile =  StringIO()
        nemooutputfile = nemobin.NemoBinaryFile(outputfile)
        nemooutputfile.write(data)
        string = outputfile.getvalue()
        outputfile.close()
        
        file = open(filename, 'rb')
        original = file.read()
        file.close()
        
        self.assertEquals(len(original), len(string))
        self.assertEquals(original, string)
        
    
    def test8(self):
        directory_name = os.path.dirname(__file__)
        filename = os.path.join(directory_name, 'plummer128.nemo')
        file = open(filename, 'rb')
        x = nemobin.NemoBinaryFileFormatProcessor()

        set = x.load_file(file)
        file.close()
        
        outputfile =  StringIO()

        y = nemobin.NemoBinaryFileFormatProcessor()
        y.set = set
        y.store_file(outputfile)
        string = outputfile.getvalue()
        outputfile.close()
        
        inputfile = StringIO(string)
        x = nemobin.NemoBinaryFileFormatProcessor()
        set = x.load_file(inputfile)
        inputfile.close()
        self.assertEquals(len(set), 128)
        self.assertAlmostRelativeEquals(set.kinetic_energy(), 0.230214395174 | nbody_system.energy, 8)
        self.assertAlmostRelativeEquals(set.potential_energy(G=nbody_system.G), -0.473503040144  | nbody_system.energy, 8)        



        
