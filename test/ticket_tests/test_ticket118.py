from amuse.support import io
from amuse.test import amusetest

import os.path

class TestsForTicket118(amusetest.TestCase):
    
    def test1(self):
        filename = os.path.join(os.path.dirname(__file__), 'FinalSnapshot.out')
        
        set = io.read_set_from_file(filename, 'dyn')
        self.assertEquals(len(set), 10)
        self.assertAlmostRelativeEquals(0.000227826766314251919, set[9].mass.number)
        self.assertAlmostRelativeEquals(0.953575109205781479,  set[0].parent.mass.number)
        
