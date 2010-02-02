from amuse.support.io import text
from amuse.support.units import units
import unittest
import StringIO
import textwrap

class CursorTests(unittest.TestCase):
    
    def test1(self):
        contents = "1\n2\n3"
        data_file = StringIO.StringIO(contents)
        instance = text.LineBasedFileCursor(data_file)
        
        self.assertEquals("1", instance.line())
        self.assertEquals("1", instance.line())
        instance.forward()
        self.assertEquals("2", instance.line())
        instance.forward()
        self.assertEquals("3", instance.line())
        self.assertEquals("3", instance.line())
    
    def test2(self):
        contents = "1\n2\n3"
        data_file = StringIO.StringIO(contents)
        instance = text.LineBasedFileCursor(data_file)
        
        self.assertFalse(instance.is_at_end())
        instance.forward()
        self.assertFalse(instance.is_at_end())
        instance.forward()
        self.assertFalse(instance.is_at_end())
        instance.forward()
        self.assertTrue(instance.is_at_end())
        
class TableFormattedTextTests(unittest.TestCase):
    
    def test1(self):
        contents = "#header\n1 2 3\n4 5 6"
        data_file = StringIO.StringIO(contents)
        instance = text.TableFormattedText("test.txt", data_file)
        instance.attribute_names = ['a', 'b', 'c']
        particles = instance.load()
        
        self.assertEquals(len(particles), 2)
        self.assertEquals(particles[0].a, 1|units.none)
        
        
class Athena3DTextTests(unittest.TestCase):
    
    def test1(self):
        contents = """\
        # Nx1 = 128
        # x1-size = 1
        # Time = 0.103125
        #
        # [1]=i-zone [2]=x1 [3]=d [4]=M1 [5]=M2 [6]=M3 [7]=P [8]=E [9]=B1c [10]=B2c [11]=B3c [12]=B1i [13]=B1i [14]=B1i
        #
          4  3.90625e-03  1.00000e+00 -8.66511e-07  4.08477e-07  1.44419e-07  6.00000e-01  2.52500e+00  1.00000e+00  1.41421e+00  5.00000e-01  1.00000e+00  1.41421e+00  5.00000e-01
        """
        data_file = StringIO.StringIO(textwrap.dedent(contents))
        instance = text.Athena3DText("test.tab", data_file)
        particles = instance.load()
        
        self.assertEquals(len(particles), 1)
        self.assertEquals(particles[0].x1, 3.90625e-03 |units.none)
        self.assertEquals(particles[0].i_zone, 4 |units.none)
        self.assertEquals(particles[0].M1, -8.66511e-07 |units.none)
        
        
        
        
        
        
