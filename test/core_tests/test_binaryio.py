from amuse.support.io import gadget
from amuse.test import amusetest
from StringIO import StringIO

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
    
