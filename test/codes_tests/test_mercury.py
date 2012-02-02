import os
import sys
import numpy as np
from amuse.test.amusetest import TestWithMPI

from amuse.community.mercury.interface import MercuryInterface, MercuryWayWard,Mercury
from amuse.ext.solarsystem import new_solar_system_for_mercury,new_solar_system
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.ic import plummer
DUMMYID=0

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

class TestMercuryInterface(TestWithMPI):

    
    def is_fortan_version_up_to_date(self):
        try:
            from amuse import config
            is_configured = hasattr(config, 'compilers')
            if is_configured:
                is_configured = hasattr(config.compilers, 'gfortran_version')
        except ImportError:
            is_configured = False
    
        if is_configured:
            if not config.compilers.gfortran_version:
                return True
            
            try:
                parts = [int(x) for x in config.compilers.gfortran_version.split('.')]
            except:
                parts = []
                
            if len(parts) < 2:
                return True
                
            return parts[0] >= 4 and parts[1] > 1
        else:
            return True
            
    def setUp(self):
        super(TestWithMPI, self).setUp()
        self.check_fortran_version()
        
    def check_fortran_version(self):
        if not self.is_fortan_version_up_to_date():
            self.skip('cannot compile, fortran module names cannot be resolved correctly in this gfortran version')
            
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
        mass,err=instance.get_central_mass(DUMMYID)
        self.assertEqual(mass,1.0)
        radius,err=instance.get_central_radius(DUMMYID)
        self.assertEqual(radius,.005)
        j2,j4,j6,err=instance.get_central_oblateness(DUMMYID)
        self.assertEqual(j2,0.0)
        self.assertEqual(j4,0.0)
        self.assertEqual(j6,0.0)
        lx,ly,lz,err=instance.get_central_spin(DUMMYID)
        self.assertEqual(lx,0.0)
        self.assertEqual(ly,0.0)
        self.assertEqual(lz,0.0)
        instance.stop()

    def test5(self):
        instance=MercuryInterface()
        instance.initialize_code()  
        instance.set_central_mass(DUMMYID, 2.0)
        instance.set_central_radius(DUMMYID,0.0625)
        instance.set_central_oblateness(DUMMYID,0.001,0.002,0.003)
        instance.set_central_spin(DUMMYID,-0.1,-0.2,-0.3)


        mass,err=instance.get_central_mass(DUMMYID)
        self.assertEqual(mass,2.0)
        radius,err=instance.get_central_radius(DUMMYID)
        self.assertEqual(radius,.0625)
        j2,j4,j6,err=instance.get_central_oblateness(DUMMYID)
        self.assertEqual(j2,0.001)
        self.assertEqual(j4,0.002)
        self.assertEqual(j6,0.003)
        lx,ly,lz,err=instance.get_central_spin(DUMMYID)
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
        instance=MercuryInterface()
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
        err=instance.evolve_model(8*100000.)
        etot2,err=instance.get_total_energy()
        self.assertAlmostEqual(etot1,etot2,7)
        instance.stop()

    def test13(self):
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
        instance.enable_stopping_condition(3)
        instance.set_stopping_condition_timeout_parameter(0.01)
        etot1,err=instance.get_total_energy()
        err=instance.evolve_model(8*1e100)
        self.assertTrue(instance.is_stopping_condition_set(3)['result']==1)
        etot2,err=instance.get_total_energy()
        self.assertAlmostEqual(etot1,etot2,7)
        instance.stop()

    def test14(self):
        instance=MercuryInterface(redirection='none')
        instance.initialize_code()  
        instance.set_central_particle_state(DUMMYID, 1.,2.,3.,4.,5.,6.,7.,8.)
        mass,radius,j2,j4,j6,lx,ly,lz,err=instance.get_central_particle_state(DUMMYID)
        self.assertEqual(mass, 1.0)
        self.assertEqual(radius, 2.0)
        self.assertEqual(j2, 3.0)
        self.assertEqual(j4, 4.0)
        self.assertEqual(j6, 5.0)
        self.assertEqual(lx, 6.0)
        self.assertEqual(ly, 7.0)
        self.assertEqual(lz, 8.0)
        instance.stop()


class TestMercury(TestWithMPI):
    
    def is_fortan_version_up_to_date(self):
        try:
            from amuse import config
            is_configured = hasattr(config, 'compilers')
            if is_configured:
                is_configured = hasattr(config.compilers, 'gfortran_version')
        except ImportError:
            is_configured = False
    
        if is_configured:
            if not config.compilers.gfortran_version:
                return True
                
            parts = [int(x) for x in config.compilers.gfortran_version.split('.')]
            
            if len(parts) < 2:
                return True
                
            return parts[0] >= 4 and parts[1] > 1
        else:
            return True
            
    def setUp(self):
        super(TestWithMPI, self).setUp()
        self.check_fortran_version()
        
    def check_fortran_version(self):
        if not self.is_fortan_version_up_to_date():
            self.skip('cannot compile, fortran module names cannot be resolved correctly in this gfortran version')
    
    
    def sun_and_earth(self):
        orbiter = datamodel.Particles(1)
        orbiter.mass = 5.97e24 | units.kg
        orbiter.density = 1.0|units.g/units.cm**3
        orbiter.position = [1.0,0.0,0.0] | units.AU
        orbiter.velocity = [0.0, 2.0*3.1415926535*1.0/365, 0.0] | units.AUd
        orbiter.angularmomentum = [1.0,0,0] | units.MSun * units.AU**2/units.day
        orbiter.celimit = 0.0 | units.none
        
        centre = datamodel.Particles(1)
        centre.mass = 1.0 | units.MSun
        centre.radius = .0000001 | units.AU
        centre.j2 = .0001|units.AU**2
        centre.j4 = .0|units.AU**4
        centre.j6 = .0|units.AU**6
        
        centre.angularmomentum = [0.0, 0.0, 0.0] | units.MSun * units.AU**2/units.day

        return centre, orbiter

    def test0(self):
        centre, orbiter = self.sun_and_earth()
        mercury = MercuryWayWard()#,debugger='xterm')
        mercury.initialize_code()
        mercury.central_particle.add_particles(centre)
        mercury.orbiters.add_particles(orbiter)
        mercury.commit_particles()

        self.assertAlmostEqual(mercury.central_particle.j4, .0|units.AU**4)
        self.assertAlmostEqual(mercury.central_particle.mass, 1.98892e+30 |units.kg, 3)
        self.assertAlmostEqual(mercury.central_particle.mass, 1.0 |units.MSun, 3)
        self.assertEquals(mercury.get_number_of_orbiters()['norbiters'],1)
        self.assertEquals(mercury.orbiters.position, [[1,0,0]] | units.AU)
        self.assertEquals(mercury.orbiters.density, 1.0|units.g/units.cm**3 )
        self.assertEquals(mercury.orbiters.angularmomentum, [[1.0, 0.0, 0.0]] | units.MSun*units.AU**2/units.day)

        mercury.evolve_model(365.24 | units.day)

        self.assertAlmostEqual(mercury.orbiters.position, [[1.0, 0.0, 0.0]] | units.AU, 1)
        self.assertAlmostEqual(mercury.kinetic_energy+mercury.potential_energy,mercury.total_energy,3)

        mercury.stop()

    def test1(self):
        centre, orbiters = new_solar_system_for_mercury()

        mercury = MercuryWayWard()
        mercury.initialize_code()

        mercury.central_particle.add_particles(centre)
        mercury.orbiters.add_particles(orbiters)
        mercury.commit_particles()

        start_pos = mercury.orbiters[2].position
        mercury.evolve_model(365.14|units.day)
        self.assertAlmostEqual(mercury.orbiters[2].position, start_pos, 1)
        mercury.stop()

    def test2(self):
        centre, orbiters = new_solar_system_for_mercury()

        mercury = MercuryWayWard()
        mercury.initialize_code()

        mercury.central_particle.add_particles(centre)
        mercury.orbiters.add_particles(orbiters)
        mercury.commit_particles()

        start_pos = mercury.orbiters[2].position
        mercury.evolve_model(16.|units.day)
        self.assertEqual(mercury.get_time(),16.|units.day)
        mercury.stop()

    def test3(self):
        solsys = new_solar_system()

        mercury = Mercury()
        mercury.initialize_code()

        mercury.particles.add_particles(solsys)
        mercury.commit_particles()

        start_pos = mercury.orbiters[2].position
        mercury.evolve_model(365.14|units.day)
        self.assertAlmostEqual(mercury.orbiters[2].position, start_pos, 1)
        mercury.stop()

