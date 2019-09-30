from amuse.test import amusetest

import os.path
from amuse import io
from amuse.units import units

class TestsForTicket118(amusetest.TestCase):
    
    def test1(self):
        filename = os.path.join(os.path.dirname(__file__), 'FinalSnapshot.out')
        set = io.read_set_from_file(filename, 'dyn')
        
        self.assertEqual(len(set), 10)
        
    def test2(self):
        filename = os.path.join(os.path.dirname(__file__), 'FinalSnapshot.out')
        
        root = io.read_set_from_file(filename, 'dyn', return_children=False)
        
        self.assertFalse(root is None)
        
    def test3(self):
        filename = os.path.join(os.path.dirname(__file__), 'FinalSnapshot.out')
        set = io.read_set_from_file(filename, 'dyn')
        print(set)
        print("set[0].parent.mass", set[0].parent.mass)
        self.assertAlmostRelativeEquals(0.000227826766314251919 * 617.75586357299929284496, set[9].mass.value_in(units.MSun), 12)
        self.assertAlmostRelativeEquals(617.75586357299929284496 * 0.953575109205781479,  set[0].parent.mass.value_in(units.MSun), 12)
        #self.assertTrue(False)
