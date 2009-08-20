import unittest
import sys

#from amuse.legacy.sse import muse_dynamics_mpi as mpi_interface
from amuse.legacy.sse import SSE_muse_interface as f2py_interface
from amuse.legacy.sse import muse_stellar_mpi as mpi_interface

from amuse.support.data.core import Particle
from amuse.support.units import units


class TestF2PYInterface(unittest.TestCase):
    class state(object):
        def __init__(self):
            self.type = 0.0
            self.zams_mass = 0.0
            self.mass = 0.0
            self.radius = 0.0
            self.luminosity  = 0.0
            self.core_mass = 0.0
            self.core_radius = 0.0
            self.envelope_mass = 0.0
            self.envelope_radius = 0.0
            self.spin = 0.0
            self.epoch = 0.0
            self.t_ms = 0.0
            self.sse_age = 0.0
        
    def test1(self):
        sse = f2py_interface
        
        metallicity = 0.02
        
        neta = 0.5
        bwind =  0.0
        hewind =  0.5
        sigma =  190.0
        
        ifflag = 0
        wdflag =  1
        bhflag =  0 
        nsflag =  1
        mxns =  3.0
        
        pts1 = 0.05
        pts2 = 0.01
        pts3 = 0.02


        status = sse.initialize(metallicity,
            neta, bwind, hewind, sigma,
            ifflag, wdflag, bhflag, nsflag, mxns,
            pts1, pts2, pts3)
        self.assertEqual(status,0)
        
    def test2(self):
        sse = f2py_interface
        
        metallicity = 0.02
        
        neta = 0.5
        bwind =  0.0
        hewind =  0.5
        sigma =  190.0
        
        ifflag = 0
        wdflag =  1
        bhflag =  0 
        nsflag =  1
        mxns =  3.0
        
        pts1 = 0.05
        pts2 = 0.01
        pts3 = 0.02
    
    
        status = sse.initialize(metallicity,
            neta, bwind, hewind, sigma,
            ifflag, wdflag, bhflag, nsflag, mxns,
            pts1, pts2, pts3)
        self.assertEqual(status,0)
        new_state = self.state()
        new_state.mass = 1.0
        new_state.zams_mass = 1.0
        new_state.type = 1.0
        new_state.age = 1e-06
        result = sse.evolve(
            new_state.type, 
            new_state.zams_mass, new_state.mass, new_state.radius, 
            new_state.luminosity, new_state.core_mass, new_state.core_radius,
            new_state.envelope_mass, new_state.envelope_radius, new_state.spin,
            new_state.epoch, new_state.t_ms, new_state.sse_age, new_state.age
        )
        updated_state = self.state()
        (updated_state.type, updated_state.zams_mass, updated_state.mass, updated_state.radius, 
            updated_state.luminosity, updated_state.core_mass, updated_state.core_radius,
            updated_state.envelope_mass, updated_state.envelope_radius, updated_state.spin,
            updated_state.epoch, updated_state.t_ms, updated_state.sse_age, updated_state.age) = result
        attributes = ('type', 'zams_mass', 'mass', 'radius', 'luminosity', 'core_mass', 'core_radius',
                'envelope_mass', 'envelope_radius', 'spin', 'epoch', 't_ms', 'sse_age', 'age')
        expected = {
            'zams_mass': '0x1.0000000000000p+0' ,
            'mass': '0x1.0000000000000p+0' ,
            'radius': '0x1.c6c8a1c793bcep-1' ,
            'luminosity': '0x1.653b1b2d0333bp-1' ,
            'core_mass': '0x0.0p+0' ,
            'core_radius': '0x0.0p+0' ,
            'envelope_mass': '0x1.0d6fc100ab510p-5' ,
            'envelope_radius': '0x1.db27631ba0e5ap-3' ,
            'spin': '0x1.07413b0522d63p+10' ,
            'epoch': '0x0.0p+0' ,
            't_ms': '0x1.57d90abe54643p+13' ,
            'sse_age': '0x1.0c6f7a0b5ed8dp-20' ,
            'age': '0x1.0c6f7a0b5ed8dp-20' ,
        };    
        for x in expected:
            self.assertEqual(float.fromhex(expected[x]),getattr(updated_state, x))        
        for x in attributes:
            value = getattr(updated_state, x)
            #if isinstance(value, float):
            #    print "'"+x+"':", value.hex() ,","
        self.assertEquals(updated_state.age, 1e-06)
        dt = sse.get_time_step(updated_state.type,
            updated_state.zams_mass, 
            updated_state.age, 
            updated_state.mass, 
            updated_state.t_ms, 
            updated_state.epoch)
        self.assertAlmostEqual(dt, 550.1565, 2)
        
class TestMPIInterface(unittest.TestCase):
    class state(object):
        def __init__(self):
            self.type = 0.0
            self.zams_mass = 0.0
            self.mass = 0.0
            self.radius = 0.0
            self.luminosity  = 0.0
            self.core_mass = 0.0
            self.core_radius = 0.0
            self.envelope_mass = 0.0
            self.envelope_radius = 0.0
            self.spin = 0.0
            self.epoch = 0.0
            self.t_ms = 0.0
            self.sse_age = 0.0
        
    def test1(self):
        sse = mpi_interface.SSE()
        
        metallicity = 0.02
        
        neta = 0.5
        bwind =  0.0
        hewind =  0.5
        sigma =  190.0
        
        ifflag = 0
        wdflag =  1
        bhflag =  0 
        nsflag =  1
        mxns =  3.0
        
        pts1 = 0.05
        pts2 = 0.01
        pts3 = 0.02
    
    
        status = sse.initialize(metallicity,
            neta, bwind, hewind, sigma,
            ifflag, wdflag, bhflag, nsflag, mxns,
            pts1, pts2, pts3)
        self.assertEqual(status,0)
        del sse
    def test2(self):
        sse = mpi_interface.SSE()
        
        metallicity = 0.02
        
        neta = 0.5
        bwind =  0.0
        hewind =  0.5
        sigma =  190.0
        
        ifflag = 0
        wdflag =  1
        bhflag =  0 
        nsflag =  1
        mxns =  3.0
        
        pts1 = 0.05
        pts2 = 0.01
        pts3 = 0.02
    
    
        status = sse.initialize(metallicity,
            neta, bwind, hewind, sigma,
            ifflag, wdflag, bhflag, nsflag, mxns,
            pts1, pts2, pts3)
        self.assertEqual(status,0)
        new_state = self.state()
        new_state.mass = 1.0
        new_state.zams_mass = 1.0
        new_state.type = 1.0
        new_state.age = 1e-06
        result = sse.evolve(
            new_state.type, 
            new_state.zams_mass, new_state.mass, new_state.radius, 
            new_state.luminosity, new_state.core_mass, new_state.core_radius,
            new_state.envelope_mass, new_state.envelope_radius, new_state.spin,
            new_state.epoch, new_state.t_ms, new_state.sse_age, new_state.age
        )
        updated_state = self.state()
        (updated_state.type,updated_state.zams_mass, updated_state.mass, updated_state.radius, 
            updated_state.luminosity, updated_state.core_mass, updated_state.core_radius,
            updated_state.envelope_mass, updated_state.envelope_radius, updated_state.spin,
            updated_state.epoch, updated_state.t_ms, updated_state.sse_age, updated_state.age) = result
        attributes = ('type', 'zams_mass', 'mass', 'radius', 'luminosity', 'core_mass', 'core_radius',
            'envelope_mass', 'envelope_radius', 'spin', 'epoch', 't_ms', 'sse_age', 'age')
        
         
        expected = {
            'zams_mass': '0x1.0000000000000p+0' ,
            'mass': '0x1.0000000000000p+0' ,
            'radius': '0x1.c6c8a1c793bcep-1' ,
            'luminosity': '0x1.653b1b2d0333bp-1' ,
            'core_mass': '0x0.0p+0' ,
            'core_radius': '0x0.0p+0' ,
            'envelope_mass': '0x1.0d6fc100ab510p-5' ,
            'envelope_radius': '0x1.db27631ba0e5ap-3' ,
            'spin': '0x1.07413b0522d63p+10' ,
            'epoch': '0x0.0p+0' ,
            't_ms': '0x1.57d90abe54643p+13' ,
            'sse_age': '0x1.0c6f7a0b5ed8dp-20' ,
            'age': '0x1.0c6f7a0b5ed8dp-20' ,
        };    
        for x in expected:
            #print x, getattr(updated_state, x).hex()
            self.assertEqual(float.fromhex(expected[x]),getattr(updated_state, x))
            
        self.assertEquals(updated_state.age, 1e-06)
        dt = sse.get_time_step(updated_state.type,
            updated_state.zams_mass, 
            updated_state.age, 
            updated_state.mass, 
            updated_state.t_ms, 
            updated_state.epoch)
        self.assertAlmostEqual(dt, 550.1565, 2)
        del sse
        
        
class TestSSE(unittest.TestCase):
    
    def test1(self):
        sse = mpi_interface.SSE()
        sse.initialize_module_with_default_parameters() 
        star = Particle(0)
        star.mass = 5 | units.MSun
        star.radius = 0.0 | units.RSun
        sse.initialize_star(star)
        previous_type = star.type.value()
        results = []
        t0 = 0 | units.Myr
        while t0 < (125 | units.Myr):
            t0 += sse.get_time_step_for_star(star)
            sse.evolve_star(star, t0)
            
            if not star.type.value() == previous_type:
                results.append((star.current_time.value(), star.mass.value(), star.type.value()))
                previous_type = star.type.value()
                
        self.assertEqual(len(results), 6)
        
        times = ( 
            104.0 | units.Myr, 
            104.4 | units.Myr, 
            104.7 | units.Myr, 
            120.1 | units.Myr,
            120.9 | units.Myr,
            121.5 | units.Myr
        )
        for result, expected in map(None, results, times):
            self.assertAlmostEqual(result[0].number, expected.number, 1)
            
        masses = ( 
            5.000 | units.MSun, 
            5.000 | units.MSun, 
            4.998 | units.MSun, 
            4.932 | units.MSun,
            4.895 | units.MSun,
            0.997 | units.MSun
        )
        for result, expected in map(None, results, masses):
            self.assertAlmostEqual(result[1].number, expected.number, 3)
         
        del sse
            
    def xtest2(self):
        sse = mpi_interface.SSE()
        sse.initialize_module_with_default_parameters() 
        star = Particle(0)
        star.initial_mass = 5 | units.MSun
        star.mass = star.initial_mass
        star.radius = 0.0 | units.RSun
        sse.initialize_star(star)
        previous_type = star.type
        results = []
        t0 = 0 | units.Myr
        while t0 < (100 | units.Myr):
            t0 += sse.get_time_step_for_star(star)
            sse.evolve_star(star, t0)
        
        total_radius = star.core_radius + star.envelope_radius
        print star
        self.assertEqual(total_radius, star.radius)
