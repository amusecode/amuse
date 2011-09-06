from amuse.test import amusetest
import StringIO
import textwrap
import os
from amuse import io
from amuse.io import text
from amuse.units import units
from amuse.units import generic_unit_system
from amuse import datamodel
class CursorTests(amusetest.TestCase):
    
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
        
class TableFormattedTextTests(amusetest.TestCase):
    
    def test1(self):
        contents = "#header\n1 2 3\n4 5 6\n       \n7 8 9\n     "
        data_file = StringIO.StringIO(contents)
        instance = text.TableFormattedText("test.txt", data_file)
        instance.attribute_names = ['a', 'b', 'c']
        particles = instance.load()
        
        self.assertEquals(len(particles), 3)
        self.assertEquals(particles[0].a, 1|units.none)
        
    def test2(self):
        p = datamodel.Particles(2)
        p.a = [1, 4] | units.m
        p.b = [2, 5] | units.m
        p.c = [3, 6] | units.m
        
        data_file = StringIO.StringIO()
        instance = text.TableFormattedText("test.txt", data_file, p)
        instance.attribute_types = [units.m, units.m, units.m]
        instance.store()
        
        contents = data_file.getvalue()
        #print contents       
        self.assertEquals("#a b c\n#m m m\n1.0 2.0 3.0\n4.0 5.0 6.0\n", contents)
    
    def test3(self):
        x = datamodel.Particles(2)
        x.mass = [1.0, 2.0] | units.MSun
        x.radius = [3.0, 4.0] | units.RSun
        io.write_set_to_file(x, "test.csv","txt", attribute_types = (units.MSun, units.RSun))
        with open("test.csv", "r") as f:
            contents = f.read()
        self.assertEquals("#mass radius\n#MSun RSun\n1.0 3.0\n2.0 4.0\n", contents)
    
        os.remove("test.csv")
        

    def test4(self):
        mass = [1.0,2.0,3.0] | generic_unit_system.mass
        length = [3.0,4.0,5.0] | generic_unit_system.length
        
        stream = StringIO.StringIO()
        output = text.TableFormattedText(stream = stream)
        output.quantities = (mass, length)
        output.store()
        contents = stream.getvalue()
        self.assertEquals("#col(0) col(1)\n#mass length\n1.0 3.0\n2.0 4.0\n3.0 5.0\n", contents)
    

class CsvFileTextTests(amusetest.TestCase):
    
    def test0(self):
        self.assertEquals(units.convert_string_to_unit('m'),         units.m)
        self.assertEquals(units.convert_string_to_unit('m * s**-2'), units.m * units.s**-2)
        self.assertEquals(units.convert_csv_string_to_unit('m, kg, s, parsec, hour'), 
            [units.m, units.kg, units.s, units.parsec, units.hour])
    
    def test1(self):
        print "Test 1: Read comma separated values (CSV) - specified attributes"
        contents = "#header\n1,2,3\n4,5,6\n7,8,9\n"
        data_stream = StringIO.StringIO(contents)
        instance = text.CsvFileText(None, data_stream)
        instance.attribute_names = ['a', 'b', 'c']
        instance.attribute_types = [units.none, units.m, units.m/units.s]
        particles = instance.load()
        self.assertEquals(len(particles), 3)
        self.assertEquals(particles.a, [1, 4, 7] | units.none)
        self.assertEquals(particles.b, [2, 5, 8] | units.m)
        self.assertEquals(particles.c, [3, 6, 9] | units.m / units.s)
    
    def test2(self):
        print "Test 2: Read comma separated values (CSV) - attributes defined in header"
        contents = ("#a, b, c\n#no_system.get('none'),system.get('S.I.').base('length'),"
            "(system.get('S.I.').base('length') / system.get('S.I.').base('time'))\n"
            "#none, m, m/s\n1,2,3\n4,5,6\n7,8,9\n")
        data_stream = StringIO.StringIO(contents)
        instance = text.CsvFileText(None, data_stream)
        particles = instance.load()
        self.assertEquals(len(particles), 3)
        self.assertEquals(particles.a, [1, 4, 7] | units.none)
        self.assertEquals(particles.b, [2, 5, 8] | units.m)
        self.assertEquals(particles.c, [3, 6, 9] | units.m / units.s)
    
    def test3(self):
        print "Test 3: Read comma separated values (CSV) - generic units"
        contents = ("#a,b,c\n"
            "#system.get('generic').base('mass'),system.get('generic').base('length'),"
            "(((system.get('generic').base('length')**2) * (system.get('generic').base('time')**-2)) * "
            "system.get('generic').base('mass'))\n"
            "#mass,length,length**2 * time**-2 * mass\n1.0,2.0,3.0\n4.0,5.0,6.0\n")
        data_stream = StringIO.StringIO(contents)
        instance = text.CsvFileText(None, data_stream)
        particles = instance.load()
        self.assertEquals(len(particles), 2)
        self.assertEquals(particles.a, [1, 4] | generic_unit_system.mass)
        self.assertEquals(particles.b, [2, 5] | generic_unit_system.length)
        self.assertEquals(particles.c, [3, 6] | generic_unit_system.energy)
        
    def test4(self):
        print "Test 4: Write comma separated values (CSV) - specified attributes"
        particles = datamodel.Particles(2)
        particles.a = [1, 4] | units.none
        particles.b = [2, 5] | units.m
        particles.c = [3, 6] | units.kg / units.m**3
        
        data_stream = StringIO.StringIO()
        instance = text.CsvFileText(None, data_stream, particles)
        instance.attribute_names = ['a', 'b']
        instance.attribute_types = [units.none, 100*units.cm]
        instance.store()
        
        contents = data_stream.getvalue()
        self.assertEquals("#a,b\n"
            "#no_system.get('none'),(100 * (0.01 * system.get('S.I.').base('length')))\n"
            "#none,100 * cm\n1.0,2.0\n4.0,5.0\n", contents)
    
    def test5(self):
        print "Test 5: Write comma separated values (CSV) - attributes defined automatically"
        particles = datamodel.Particles(2)
        particles.a = [1, 4] | units.none
        particles.b = [2, 5] | units.m
        particles.c = [3, 6] | units.kg / units.m**3
        
        data_stream = StringIO.StringIO()
        instance = text.CsvFileText(None, data_stream, particles)
        instance.store()
        
        contents = data_stream.getvalue()
        self.assertEquals("#a,b,c\n"
            "#no_system.get('none'),system.get('S.I.').base('length'),"
            "((system.get('S.I.').base('length')**-3) * system.get('S.I.').base('mass'))\n"
            "#none,m,m**-3 * kg\n1.0,2.0,3.0\n4.0,5.0,6.0\n", contents)
    
    def test6(self):
        print "Test 6: Write comma separated values (CSV) - generic units"
        particles = datamodel.Particles(2)
        particles.a = [1, 4] | generic_unit_system.mass
        particles.b = [2, 5] | generic_unit_system.length
        particles.c = [3, 6] | generic_unit_system.energy
        particles.d = [4, 7] | generic_unit_system.temperature
        
        data_stream = StringIO.StringIO()
        instance = text.CsvFileText(None, data_stream, particles)
        instance.store()
        
        contents = data_stream.getvalue()
        self.assertEquals("#a,b,c,d\n"
            "#system.get('generic').base('mass'),system.get('generic').base('length'),"
            "(((system.get('generic').base('length')**2) * (system.get('generic').base('time')**-2)) * "
            "system.get('generic').base('mass')),system.get('generic').base('thermodynamic temperature')\n"
            "#mass,length,length**2 * time**-2 * mass,thermodynamic temperature\n1.0,2.0,3.0,4.0\n4.0,5.0,6.0,7.0\n", contents)
    
    def test7(self):
        print "Test 7: User interface (write_set_to_file and read_set_from_file)"
        particles = datamodel.Particles(2)
        particles.a = [1, 4] | units.none
        particles.b = [2, 5] | units.m
        particles.c = [3, 6] | units.kg / units.m**3
        io.write_set_to_file(particles, "test_textio.csv","csv")
        
        read_particles = io.read_set_from_file("test_textio.csv", format = "csv")
        self.assertEquals(len(read_particles), 2)
        self.assertEquals(read_particles.a, [1, 4] | units.none)
        self.assertEquals(read_particles.b, [2, 5] | units.m)
        self.assertEquals(read_particles.c, [3, 6] | units.kg / units.m**3)
        os.remove("test_textio.csv")
    

class Athena3DTextTests(amusetest.TestCase):
    
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
        
        
        
        
        
        
