from amuse.test import amusetest
from StringIO import StringIO

import os.path
from amuse import io
from amuse.io import gadget
from amuse.io import nemobin
from amuse.units import nbody_system
from amuse.datamodel import Particles

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
        self.assertEquals(x.header_struct.Npart[3], 0)
        self.assertEquals(x.header_struct.Npart[4], 0)
        self.assertEquals(x.header_struct.Npart[5], 0)
        
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
        self.assertEquals(halo[1].key,544511335)
        self.assertEquals(halo[2].key,544511457)
        self.assertAlmostRelativeEquals(halo[0].velocity[0], -24.785614 |  nbody_system.speed, 7)
        print halo[1].velocity
        self.assertAlmostRelativeEquals(halo[1].velocity[0], -25.346375 |  nbody_system.speed, 7)
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
        self.assertEquals(halo[1].id,544511335)
        self.assertEquals(halo[2].id,544511457)
        self.assertAlmostRelativeEquals(halo[0].velocity[0], -24.785614 |  nbody_system.speed, 7)
        self.assertAlmostRelativeEquals(halo[1].velocity[0], -25.346375 |  nbody_system.speed, 7)
        self.assertAlmostRelativeEquals(halo[2].velocity[0], -25.394440 |  nbody_system.speed, 7)
        
    
    def test9(self):
        class FakeList(object):
            def __init__(self, _len):
                self._len = _len
            def __len__(self):
                return self._len
                
        set = (FakeList(20000), FakeList(500000), FakeList(20000), (), (), ())
        x = gadget.GadgetFileFormatProcessor(set = set)
        x.equal_mass_array=(0.0,4.291150104743886e-07, 2.5e-07,0.0,0.0,0.0) |nbody_system.mass
        file = StringIO()
        x.store_header(file)
        print x.header_struct
        self.assertEquals(x.header_struct.Npart[0], 20000)
        self.assertEquals(x.header_struct.Npart[1], 500000)
        self.assertEquals(x.header_struct.Npart[2], 20000)
        self.assertEquals(x.header_struct.Npart[3], 0)
        self.assertEquals(x.header_struct.Npart[4], 0)
        self.assertEquals(x.header_struct.Npart[5], 0)
        
        print repr(file.getvalue())
        print repr(''.join(self.header_parts))
        self.assertEquals(repr(file.getvalue()[0:30]), repr(''.join(self.header_parts)[0:30]))
    
    def test10(self):
        p = Particles(2)
        p[0].position = [1.0, 2.0, 3.0] | nbody_system.length
        p[1].position = [4.0, 5.0, 6.0] | nbody_system.length
        p[0].velocity = [7.0, 8.0, 10.0] | nbody_system.length / nbody_system.time
        p[1].velocity = [11.0, 12.0, 13.0] | nbody_system.length / nbody_system.time
        p.u = [3,4] | nbody_system.potential
        p.rho = [5,6] | nbody_system.density
        p.mass = [5,6] | nbody_system.mass
        x = gadget.GadgetFileFormatProcessor(set = p)
        file = StringIO()
        x.store_body(file)
        input = StringIO(file.getvalue())
        positions = x.read_fortran_block_float_vectors(input)
        self.assertEquals(positions[0] , [1.0, 2.0, 3.0])
        self.assertEquals(positions[1] , [4.0, 5.0, 6.0])
        velocities = x.read_fortran_block_float_vectors(input)
        self.assertEquals(velocities[0] , [7.0, 8.0, 10.0])
        self.assertEquals(velocities[1] , [11.0, 12.0, 13.0])
        ids = x.read_fortran_block_ulongs(input)
        self.assertEquals(ids[0], p[0].key)
        self.assertEquals(ids[1], p[1].key)
        masses = x.read_fortran_block_floats(input)
        self.assertEquals(masses[0], 5)
        self.assertEquals(masses[1], 6)
        u = x.read_fortran_block_floats(input)
        self.assertEquals(u[0], 3)
        self.assertEquals(u[1], 4)
        
    def test11(self):
        directory_name = os.path.dirname(__file__)
        filename = os.path.join(directory_name, 'gassphere_littleendian.dat')
        gas, halo, disk, bulge, stars, bndry = io.read_set_from_file(filename, format='gadget')
        self.assertEquals(len(gas), 1472)
        self.assertEquals(len(halo), 0)
        self.assertEquals(gas[0].key,1)
        self.assertEquals(gas[1].key,2)
        self.assertEquals(gas[2].key,3)
        self.assertEquals(gas[1471].key,1472)
        self.assertAlmostRelativeEquals(gas[0:5].x,[-0.0713372901082, 0.0713372901082, -0.21178227663, -0.0698266476393, 0.0698266476393] | nbody_system.length, 7)
        self.assertAlmostRelativeEquals(gas[0:5].u, [0.0500000007451, 0.0500000007451, 0.0500000007451, 0.0500000007451, 0.0500000007451] | (nbody_system.length / nbody_system.time)**2, 7 )
       
        outputfilename = 'gadgettest.output'
        try:
            io.write_set_to_file((gas, halo, disk, bulge, stars, bndry), outputfilename, format='gadget', ids_are_long = False)
        
            gas, halo, disk, bulge, stars, bndry = io.read_set_from_file(filename, format='gadget')
            self.assertEquals(len(gas), 1472)
            self.assertEquals(len(halo), 0)
            self.assertEquals(gas[0].key,1)
            self.assertEquals(gas[1].key,2)
            self.assertEquals(gas[2].key,3)
            self.assertEquals(gas[1471].key,1472)
        finally:
            if os.path.exists(outputfilename):
                os.remove(outputfilename)
        
        
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
 
