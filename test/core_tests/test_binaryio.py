from amuse.test import amusetest
from StringIO import StringIO

import os.path
from amuse import io
from amuse.io import gadget
from amuse.io import nemobin
from amuse.units import nbody_system

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
        
    def test5(self):
        options = io.get_options_for_format('gadget')
        found_has_acceleration = False
        for name, description, defaultval in options:
            if name == 'has_acceleration':
                found_has_acceleration = True   
        
        self.assertTrue(found_has_acceleration)
        
    def test6(self):
        directory_name = os.path.dirname(__file__)
        filename = os.path.join(directory_name, 'gassphere_littleendian.dat')
        x = gadget.GadgetFileFormatProcessor()
        gas, halo, disk, bulge, stars, bndry =  io.read_set_from_file(filename, format='gadget')
        self.assertEquals(len(gas), 1472)
        self.assertEquals(len(halo), 0)
        self.assertEquals(gas[0].key,1)
        self.assertEquals(gas[1].key,2)
        self.assertEquals(gas[2].key,3)
        self.assertEquals(gas[1471].key,1472)
        
    def test7(self):
        """test returned ids from gadget file
        for ticket #245.
        All the 'uneven' particles have key "1", and identical velocities/positions. This is incorrect
        upon further inspection, the test file is incorrect
        """
        directory_name = os.path.dirname(__file__)
        filename = os.path.join(directory_name, 'ticket245.dat')
        gas, halo, disk, bulge, stars, bndry = io.read_set_from_file(filename, format='gadget')
        
        self.assertEquals(len(gas), 0)
        self.assertEquals(len(halo),1324)
        self.assertEquals(len(disk), 0)
        self.assertEquals(len(bulge), 0)
        self.assertEquals(len(stars), 0)
        self.assertEquals(len(bndry), 0)
        self.assertEquals(halo[0].key,544418538)
        self.assertEquals(halo[1].key,0)
        self.assertEquals(halo[2].key,544511335)
        self.assertAlmostRelativeEquals(halo[0].velocity[0], -24.785614 |  nbody_system.speed, 7)
        print halo[1].velocity
        self.assertAlmostRelativeEquals(halo[1].velocity[0], -26.2435913086 |  nbody_system.speed, 7)
        self.assertAlmostRelativeEquals(halo[2].velocity[0], -25.394440 |  nbody_system.speed, 7)
        
    def test8(self):
        """test returned ids from gadget file
        for ticket #245.
        added option to not use the ids as a key, should fix the problem
        for incorrect id's
        """
        directory_name = os.path.dirname(__file__)
        filename = os.path.join(directory_name, 'ticket245.dat')
        gas, halo, disk, bulge, stars, bndry = io.read_set_from_file(filename, format='gadget', ids_are_keys = False)
        
        self.assertEquals(len(gas), 0)
        self.assertEquals(len(halo),1324)
        self.assertEquals(len(disk), 0)
        self.assertEquals(len(bulge), 0)
        self.assertEquals(len(stars), 0)
        self.assertEquals(len(bndry), 0)
        self.assertEquals(halo[0].id,544418538)
        self.assertEquals(halo[1].id,0)
        self.assertEquals(halo[2].id,544511335)
        self.assertAlmostRelativeEquals(halo[0].velocity[0], -24.785614 |  nbody_system.speed, 7)
        print halo[1].velocity
        self.assertAlmostRelativeEquals(halo[1].velocity[0], -26.2435913086 |  nbody_system.speed, 7)
        self.assertAlmostRelativeEquals(halo[2].velocity[0], -25.394440 |  nbody_system.speed, 7)
        
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

