import os
import sys
import numpy
from amuse.test.amusetest import TestWithMPI

from amuse.community.mercury.interface import MercuryInterface, Mercury
from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.ext import plummer

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

class TestMPIInterface(TestWithMPI):

    def test1(self):
        instance=MercuryInterface()  
        instance.stop()

    def test2(self):
        instance=MercuryInterface()
        instance.initialize_code()  
        instance.stop()

    def test3(self):
        instance=MercuryInterface()
        instance.initialize_code()  
        instance.cleanup_code()  
        instance.stop()

    def test4(self):
        instance=MercuryInterface()
        instance.initialize_code()  
        mass,err=instance.get_central_mass()
        self.assertEqual(mass,1.0)
        radius,err=instance.get_central_radius()
        self.assertEqual(radius,.005)
        j2,j4,j6,err=instance.get_central_oblateness()
        self.assertEqual(j2,0.0)
        self.assertEqual(j4,0.0)
        self.assertEqual(j6,0.0)
        lx,ly,lz,err=instance.get_central_spin()
        self.assertEqual(lx,0.0)
        self.assertEqual(ly,0.0)
        self.assertEqual(lz,0.0)
        instance.stop()

    def test5(self):
        instance=MercuryInterface()
        instance.initialize_code()  
        instance.set_central_mass(2.0)
        instance.set_central_radius(0.0625)
        instance.set_central_oblateness(0.001,0.002,0.003)
        instance.set_central_spin(-0.1,-0.2,-0.3)


        mass,err=instance.get_central_mass()
        self.assertEqual(mass,2.0)
        radius,err=instance.get_central_radius()
        self.assertEqual(radius,.0625)
        j2,j4,j6,err=instance.get_central_oblateness()
        self.assertEqual(j2,0.001)
        self.assertEqual(j4,0.002)
        self.assertEqual(j6,0.003)
        lx,ly,lz,err=instance.get_central_spin()
        self.assertEqual(lx,-0.1)
        self.assertEqual(ly,-0.2)
        self.assertEqual(lz,-0.3)
        instance.stop()


    def test6(self):
        instance=MercuryInterface()
        instance.initialize_code()
        instance.commit_particles()
        instance.cleanup_code()  
        instance.stop()

    def test7(self):
        instance=MercuryInterface()
        instance.initialize_code()  
        time,err=instance.get_time()  
        self.assertEqual(time,0.0)
        instance.stop()

    def test8(self):
        instance=MercuryInterface()
        instance.initialize_code()  
        n,err=instance.get_number_of_orbiters()  
        self.assertEqual(n,0)
        instance.stop()

    def test9(self):
        instance=MercuryInterface(redirection='none')
        instance.initialize_code()
        pid,err=instance.new_orbiter(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)  
        n,err=instance.get_number_of_orbiters()  
        self.assertEqual(n,1)
        instance.stop()

    def test10(self):
        instance=MercuryInterface()
        instance.initialize_code()
        pid,err=instance.new_orbiter(0.5,1.,2.,3.,4.,5.,6.,7.,1.,2.,3.,8.)  
        mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit,err=instance.get_orbiter_state(pid)
        self.assertEqual(mass,0.5)       
        self.assertEqual(dens,1.)       
        self.assertEqual(x,2.)       
        self.assertEqual(y,3.)       
        self.assertEqual(z,4.)
        self.assertEqual(vx,5.)
        self.assertEqual(vy,6.)
        self.assertEqual(vz,7.)
        self.assertEqual(sx,1.)
        self.assertEqual(sy,2.)
        self.assertEqual(sz,3.)
        self.assertEqual(celimit,8.)       
        instance.stop()

    def test11(self):
        instance=MercuryInterface()
        instance.initialize_code()
        mass=3.04043264264672381E-06
        dens=5.52
        x=2.42093942183383037E-01
        y=-9.87467766698604366E-01
        z=-4.54276292555233496E-06
        vx=1.64294055023289365E-02
        vy=4.03200725816140870E-03
        vz=1.13609607260006795E-08
        sx=sy=sz=0.
        celimit=20.
        pid,err=instance.new_orbiter(mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit)  
        instance.commit_particles()
        de,err=instance.get_energy_deviation()
        self.assertEqual(de,0.0)
        etot,err=instance.get_total_energy()
        ep,err=instance.get_potential_energy()
        ek,err=instance.get_kinetic_energy()
        self.assertAlmostEqual(etot,ep+ek,15)
        instance.stop()

    def test12(self):
        instance=MercuryInterface()
        instance.initialize_code()
        mass=3.04043264264672381E-06
        dens=5.52
        x=2.42093942183383037E-01
        y=-9.87467766698604366E-01
        z=-4.54276292555233496E-06
        vx=1.64294055023289365E-02
        vy=4.03200725816140870E-03
        vz=1.13609607260006795E-08
        sx=sy=sz=0.
        celimit=20.
        pid,err=instance.new_orbiter(mass,dens,x,y,z,vx,vy,vz,sx,sy,sz,celimit)  
        instance.commit_particles()
        etot1,err=instance.get_total_energy()
        err=instance.evolve(8*100000.)
        etot2,err=instance.get_total_energy()
        self.assertAlmostEqual(etot1,etot2,7)
        instance.stop()

class TestMercury(TestWithMPI):
    def new_system_of_sun_and_earth(self):
        stars = core.Stars(2)
        sun = stars[0]
        sun.mass = units.MSun(1.0)
        sun.position = units.m(numpy.array((0.0,0.0,0.0)))
        sun.velocity = units.ms(numpy.array((0.0,0.0,0.0)))
        sun.radius = units.RSun(1.0)
        
        earth = stars[1]
        earth.mass = units.kg(5.9736e24)
        earth.radius = units.km(6371) 
        earth.position = units.km(numpy.array((149.5e6,0.0,0.0)))
        earth.velocity = units.ms(numpy.array((0.0,29800,0.0)))
        
        return stars
    
    def test0(self):
        orbiter = core.Particles(1)
        orbiter.mass = 1.0 | units.MSun
        orbiter.density = 1.0|units.g/units.cm**3
        orbiter.position = [0.0,0.0,0.0] | units.AU
        orbiter.velocity = [100.0,100.0,100.0] | units.AUd
        orbiter.sx = 0.0|units.day
        orbiter.sy = 0.0|units.day
        orbiter.sz = 0.0|units.day
        
        orbiter.celimit = 0.0|units.none

        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        mercury = Mercury(convert_nbody)
        mercury.initialize_code()
        mercury.commit_parameters()
        mercury.particles.add_particles(orbiter)
        self.assertEquals(mercury.get_number_of_orbiters()['norbiters'],1)
        mercury.commit_particles()
        mercury.stop()

    def xtest1(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
    
        hermite = Hermite(convert_nbody)
        hermite.initialize_code()
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
                
        hermite.particles.add_particles(stars)
        
        hermite.evolve_model(365.0 | units.day)
        hermite.particles.copy_values_of_state_attributes_to(stars)
        
        position_at_start = earth.position.value_in(units.AU)[0]
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(position_at_start, position_after_full_rotation, 6)
        
        hermite.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        hermite.particles.copy_values_of_state_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 2)
                
        hermite.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        hermite.particles.copy_values_of_state_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
        
        hermite.cleanup_code()
        
        hermite.stop()
    
