import os
import sys
import numpy as np
from amuse.test.amusetest import TestWithMPI

from amuse.community.mercury.interface import MercuryInterface, MercuryWayWard,Mercury
from amuse.ext.solarsystem import new_solar_system_for_mercury,new_solar_system
from amuse.units import nbody_system
from amuse.units import units, constants
from amuse import datamodel
from amuse.ic import plummer
DUMMYID=0

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

class TestMercuryInterface(TestWithMPI):
    
    def setUp(self):
        super(TestWithMPI, self).setUp()
            
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

    def test15(self):
        instance=MercuryInterface(redirection='none')
        instance.initialize_code()  
        i,err=instance.get_integrator()
        self.assertEqual(i,10)
        err=instance.set_integrator(2)
        self.assertEqual(err,0)
        i,err=instance.get_integrator()
        self.assertEqual(i,2)
        err=instance.set_integrator(22)
        self.assertEqual(err,-1)
        instance.stop()

    def test16(self):
        instance=MercuryInterface()
        instance.initialize_code()
        err=instance.set_integrator(2)
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


class TestMercury(TestWithMPI):
                
    def setUp(self):
        super(TestWithMPI, self).setUp()
        
    def sun_and_earth(self):
        orbiter = datamodel.Particles(1)
        orbiter.mass = 5.97e24 | units.kg
        orbiter.density = 1.0|units.g/units.cm**3
        orbiter.position = [1.0,0.0,0.0] | units.AU
        orbiter.velocity = [0.0, 2.0*3.1415926535*1.0/365, 0.0] | units.AUd
        orbiter.angular_momentum = [1.0,0,0] | units.MSun * units.AU**2/units.day
        orbiter.celimit = 0.0 | units.none
        
        centre = datamodel.Particles(1)
        centre.mass = 1.0 | units.MSun
        centre.radius = .0000001 | units.AU
        centre.j2 = .0001|units.AU**2
        centre.j4 = .0|units.AU**4
        centre.j6 = .0|units.AU**6
        
        centre.angular_momentum = [0.0, 0.0, 0.0] | units.MSun * units.AU**2/units.day

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
        self.assertEqual(mercury.get_number_of_orbiters(),1)
        self.assertEqual(mercury.orbiters.position, [[1,0,0]] | units.AU)
        self.assertEqual(mercury.orbiters.density, 1.0|units.g/units.cm**3 )
        self.assertEqual(mercury.orbiters.angular_momentum, [[1.0, 0.0, 0.0]] | units.MSun*units.AU**2/units.day)

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
            
    def test4(self):
        solsys = new_solar_system()

        mercury = Mercury()
        mercury.initialize_code()
        self.assertEqual(mercury.parameters.timestep, 8 | units.day)
        mercury.set_initial_timestep(1 | units.day)
        mercury.parameters.timestep = 1 | units.day
        
        mercury.particles.add_particles(solsys)
        mercury.commit_particles()

        start_pos = mercury.orbiters[2].position
        mercury.evolve_model(365.14|units.day)
        self.assertAlmostEqual(mercury.orbiters[2].position, start_pos, 1)
        mercury.stop()

    def test5(self):
        centre, orbiters = new_solar_system_for_mercury()

        mercury = MercuryWayWard()
        self.assertEqual(mercury.get_name_of_current_state(), 'UNINITIALIZED')
        mercury.initialize_code()
        self.assertEqual(mercury.get_name_of_current_state(), 'INITIALIZED')
        mercury.commit_parameters()
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')
        
        mercury.central_particle.add_particles(centre)
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')
        mercury.orbiters.add_particles(orbiters)
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')

        mercury.commit_particles()
        self.assertEqual(mercury.get_name_of_current_state(), 'RUN')

        start_pos = mercury.orbiters[2].position
        mercury.evolve_model(365.14|units.day)
        self.assertEqual(mercury.get_name_of_current_state(), 'EVOLVED')

        self.assertAlmostEqual(mercury.orbiters[2].position, start_pos, 1)
        mercury.cleanup_code()
        self.assertEqual(mercury.get_name_of_current_state(), 'END')
        mercury.stop()

    def test6(self):
        centre, orbiters = new_solar_system_for_mercury()

        mercury = MercuryWayWard()
        self.assertEqual(mercury.get_name_of_current_state(), 'UNINITIALIZED')
        mercury.initialize_code()
        self.assertEqual(mercury.get_name_of_current_state(), 'INITIALIZED')
        mercury.commit_parameters()
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')
        
        mercury.central_particle.add_particles(centre)
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')
        mercury.orbiters.add_particles(orbiters[4:5])
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')
        mercury.commit_particles()
        self.assertEqual(mercury.get_name_of_current_state(), 'RUN')

        start_pos = mercury.orbiters[0].position
        mercury.evolve_model(11.8618|units.yr)
        self.assertEqual(mercury.get_name_of_current_state(), 'EVOLVED')
        self.assertAlmostEqual(mercury.orbiters[0].position, start_pos, 1)

        mercury.orbiters.add_particles(orbiters[0:1])
        self.assertEqual(mercury.get_name_of_current_state(), 'UPDATE')
        mercury.recommit_particles()

        self.assertAlmostEqual(mercury.orbiters[0].position, start_pos, 1)

        start_pos = mercury.orbiters[0].position
        mercury.evolve_model(2*11.8618|units.yr)
        self.assertEqual(mercury.get_name_of_current_state(), 'EVOLVED')
        self.assertAlmostEqual(mercury.orbiters[0].position, start_pos, 1)

        mercury.cleanup_code()
        self.assertEqual(mercury.get_name_of_current_state(), 'END')
        mercury.stop()

    def test7(self):
        centre, orbiters = new_solar_system_for_mercury()

        mercury = MercuryWayWard()#debugger="gdb")
        self.assertEqual(mercury.get_name_of_current_state(), 'UNINITIALIZED')
        mercury.initialize_code()
        self.assertEqual(mercury.get_name_of_current_state(), 'INITIALIZED')
        mercury.commit_parameters()
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')
        
        mercury.central_particle.add_particles(centre)
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')
        mercury.orbiters.add_particles(orbiters[0:5])
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')
        mercury.commit_particles()
        self.assertEqual(mercury.get_name_of_current_state(), 'RUN')

        start_pos = mercury.orbiters[4].position
        mercury.evolve_model(11.8618|units.yr)
        self.assertEqual(mercury.get_name_of_current_state(), 'EVOLVED')
        self.assertAlmostEqual(mercury.orbiters[4].position, start_pos, 1)

        mercury.orbiters.remove_particles(orbiters[0:1])
        self.assertEqual(mercury.get_name_of_current_state(), 'UPDATE')
        mercury.recommit_particles()

        self.assertAlmostEqual(mercury.orbiters[3].position, start_pos, 1)

        start_pos = mercury.orbiters[3].position
        mercury.evolve_model(2*11.8618|units.yr)
        self.assertEqual(mercury.get_name_of_current_state(), 'EVOLVED')
        self.assertAlmostEqual(mercury.orbiters[3].position, start_pos, 1)

        mercury.cleanup_code()
        self.assertEqual(mercury.get_name_of_current_state(), 'END')
        mercury.stop()


    def test8(self):
        centre, orbiters = new_solar_system_for_mercury()

        mercury = MercuryWayWard()#debugger="gdb")
        self.assertEqual(mercury.get_name_of_current_state(), 'UNINITIALIZED')
        mercury.initialize_code()
        self.assertEqual(mercury.get_name_of_current_state(), 'INITIALIZED')
        mercury.commit_parameters()
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')
        
        mercury.central_particle.add_particles(centre)
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')
        mercury.orbiters.add_particles(orbiters[0:5])
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')
        mercury.commit_particles()
        self.assertEqual(mercury.get_name_of_current_state(), 'RUN')

        start_pos = mercury.orbiters[4].position
        mercury.evolve_model(11.8618|units.yr)
        self.assertEqual(mercury.get_name_of_current_state(), 'EVOLVED')
        self.assertAlmostEqual(mercury.orbiters[4].position, start_pos, 1)

        mercury.orbiters.remove_particles(orbiters[0:4])
        self.assertEqual(mercury.get_name_of_current_state(), 'UPDATE')
        mercury.recommit_particles()

        self.assertAlmostEqual(mercury.orbiters[0].position, start_pos, 1)

        start_pos = mercury.orbiters[0].position
        mercury.evolve_model(2*11.8618|units.yr)
        self.assertEqual(mercury.get_name_of_current_state(), 'EVOLVED')
        self.assertAlmostEqual(mercury.orbiters[0].position, start_pos, 1)

        mercury.cleanup_code()
        self.assertEqual(mercury.get_name_of_current_state(), 'END')
        mercury.stop()

    def test9(self):
        centre, orbiters = new_solar_system_for_mercury()

        mercury = MercuryWayWard()#debugger="gdb")
        self.assertEqual(mercury.get_name_of_current_state(), 'UNINITIALIZED')
        mercury.central_particle.add_particles(centre)
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')
        mercury.orbiters.add_particles(orbiters[0:5])
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')
        start_pos = mercury.orbiters[4].position
        mercury.evolve_model(11.8618|units.yr)
        self.assertEqual(mercury.get_name_of_current_state(), 'EVOLVED')
        self.assertAlmostEqual(mercury.orbiters[4].position, start_pos, 1)

        mercury.orbiters.remove_particles(orbiters[0:4])
        self.assertEqual(mercury.get_name_of_current_state(), 'UPDATE')
        self.assertAlmostEqual(mercury.orbiters[0].position, start_pos, 1)

        start_pos = mercury.orbiters[0].position
        mercury.evolve_model(2*11.8618|units.yr)
        self.assertEqual(mercury.get_name_of_current_state(), 'EVOLVED')
        self.assertAlmostEqual(mercury.orbiters[0].position, start_pos, 1)

    def test10(self):
        centre, orbiters = new_solar_system_for_mercury()

        mercury = MercuryWayWard()
        mercury.central_particle.add_particles(centre)
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')
        mercury.orbiters.add_particles(orbiters[4:5])
        self.assertEqual(mercury.get_name_of_current_state(), 'EDIT')
        start_pos = mercury.orbiters[0].position
        mercury.evolve_model(11.8618|units.yr)
        self.assertEqual(mercury.get_name_of_current_state(), 'EVOLVED')
        self.assertAlmostEqual(mercury.orbiters[0].position, start_pos, 1)

        mercury.orbiters.add_particles(orbiters[0:4])
        self.assertEqual(mercury.get_name_of_current_state(), 'UPDATE')
        self.assertAlmostEqual(mercury.orbiters[0].position, start_pos, 1)
        start_pos = mercury.orbiters[0].position
        mercury.evolve_model(2*11.8618|units.yr)
        self.assertEqual(mercury.get_name_of_current_state(), 'EVOLVED')
        self.assertAlmostEqual(mercury.orbiters[0].position, start_pos, 1)

    def test11(self):
        solsys = new_solar_system()

        mercury = Mercury()
        self.assertEqual(mercury.parameters.timestep, 8 | units.day)
        mercury.parameters.timestep = 1 | units.day
        self.assertEqual(mercury.parameters.timestep, 1 | units.day)
        
        mercury.particles.add_particles(solsys)
        start_pos = mercury.particles[5].position
        mercury.evolve_model(11.8618|units.yr)
        self.assertAlmostEqual(mercury.particles[5].position, start_pos, 1)
        mercury.stop()

    def test12(self):
        solsys = new_solar_system()

        mercury = Mercury()
        mercury.parameters.timestep = 1. | units.day
        
        mercury.particles.add_particles(solsys)
        start_pos = mercury.particles[5].position
        mercury.evolve_model(11.8618|units.yr)
        self.assertAlmostEqual(mercury.particles[5].position, start_pos, 1)
        mercury.particles.remove_particles(mercury.particles[1:5])
        self.assertAlmostEqual(mercury.particles[1].position, start_pos, 1)
        start_pos = mercury.particles[1].position
        mercury.evolve_model(2*11.8618|units.yr)
        self.assertAlmostEqual(mercury.particles[1].position, start_pos, 1)
        mercury.stop()

    def test13(self):
        solsys = new_solar_system()

        mercury = Mercury()        
        mercury.particles.add_particles(solsys)
        idpos1 = [ (p.position - q.position) for p in mercury.particles[1:10] for q in mercury.particles[1:10] ]
        mercury.evolve_model(11.8618|units.yr)
        edpos1 = [ (p.position - q.position) for p in mercury.particles[1:10] for q in mercury.particles[1:10] ]
        mercury.stop()

        centre, orbiters = new_solar_system_for_mercury()

        mercury = MercuryWayWard()
        mercury.central_particle.add_particles(centre)
        mercury.orbiters.add_particles(orbiters)

        idpos2 = [ (p.position - q.position) for p in mercury.orbiters[0:9] for q in mercury.orbiters[0:9] ]
        mercury.evolve_model(11.8618|units.yr)
        edpos2 = [ (p.position - q.position) for p in mercury.orbiters[0:9] for q in mercury.orbiters[0:9] ]
        mercury.stop()

        for d1,d2 in zip(idpos1,idpos2):
          self.assertAlmostEqual(d1,d2, 7)
        for d1,d2 in zip(edpos1,edpos2):
          self.assertAlmostEqual(d1,d2, 7)

    def test14(self):
        centre, orbiters = new_solar_system_for_mercury()
        oneyear = 365.14 | units.day
        halfyear = oneyear / 2.0
        mercury = MercuryWayWard()
        mercury.initialize_code()
        mercury.commit_parameters()

        mercury.central_particle.add_particles(centre)
        mercury.orbiters.add_particles(orbiters)

        start_pos = mercury.orbiters[2].position
        mercury.evolve_model(halfyear)
        central_particles = mercury.central_particle.copy()
        orbiters = mercury.orbiters.copy()
        mercury.stop()
        
        mercury = MercuryWayWard()
        mercury.initialize_code()
        mercury.parameters.begin_time = halfyear
        mercury.commit_parameters()

        mercury.central_particle.add_particles(centre)
        mercury.orbiters.add_particles(orbiters)
        mercury.evolve_model(oneyear)
        self.assertAlmostEqual(mercury.orbiters[2].position, start_pos, 1)
        mercury.stop()

    def test15(self):
        solsys = new_solar_system()

        solsys.x-=1.| units.AU

        p=datamodel.Particles(3)
        p.mass=[1,2,3] | units.MSun
        p.x=[1,10,100] | units.AU
        p.y=[0,0,-10] | units.AU
        p.z=[0,0,10] | units.AU

        pe=p.potential_energy_in_field(solsys)

        mercury = Mercury()
        mercury.particles.add_particles(solsys)
        pot=mercury.get_potential_at_point(p.x*0.,p.x,p.y,p.z)

        self.assertAlmostRelativeEqual((pot*p.mass).sum(),pe,12)

    def test16(self):
        solsys = new_solar_system()

        solsys.x-=1.| units.AU

        p=datamodel.Particles(3)
        p.mass=[1,2,3] | units.MSun
        p.x=[1,10,100] | units.AU
        p.y=[0,0,-10] | units.AU
        p.z=[0,0,10] | units.AU

        from amuse.community.huayno.interface import Huayno
        from amuse.units import nbody_system        
        conv=nbody_system.nbody_to_si(1. | units.AU, 1.| units.MSun)
        h = Huayno(conv)
        h.particles.add_particles(solsys)
        
        ax1,ay1,az1=h.get_gravity_at_point(p.x*0.,p.x,p.y,p.z)
        
        mercury = Mercury()        
        mercury.particles.add_particles(solsys)

        ax2,ay2,az2=mercury.get_gravity_at_point(p.x*0.,p.x,p.y,p.z)


        self.assertAlmostRelativeEqual(ax1,ax2,12)
        self.assertAlmostRelativeEqual(ay1,ay2,12)
        self.assertAlmostRelativeEqual(az1,az2,12)

    def test17(self):
        solsys = new_solar_system()

        mercury = Mercury()
        mercury.particles.add_particles(solsys)
        mercury.evolve_model(5.|units.yr)
        mercury.evolve_model(10.|units.yr)

        dpos1=mercury.particles[0].position-mercury.particles[5].position

        mercury.stop()

        mercury = Mercury()
        mercury.particles.add_particles(solsys)
        mercury.evolve_model(5.|units.yr)
        mercury.particles.x+=1 | units.AU
        mercury.particles.y+=2 | units.AU
        mercury.particles.z+=3 | units.AU
        mercury.particles.vx+=10. | units.kms
        mercury.particles.vy+=20. | units.kms
        mercury.particles.vz+=30. | units.kms
        mercury.evolve_model(10.|units.yr)

        dpos2=mercury.particles[0].position-mercury.particles[5].position
        mercury.stop()
                
        self.assertAlmostEqual(dpos1, dpos2, 12)

    def test18(self):
        solsys = new_solar_system()

        mercury = Mercury()
        mercury.initialize_code()
        self.assertEqual(mercury.parameters.integrator,10)
        mercury.parameters.integrator=2
        self.assertEqual(mercury.parameters.integrator,2)
        mercury.stop()

    def test19(self):
        solsys = new_solar_system()

        mercury = Mercury()
        mercury.initialize_code()
        mercury.parameters.integrator=1

        mercury.particles.add_particles(solsys)
        mercury.commit_particles()

        start_pos = mercury.orbiters[2].position
        mercury.evolve_model(365.14|units.day)
        self.assertAlmostEqual(mercury.orbiters[2].position, start_pos, 1)
        mercury.stop()

    def test20(self):
        solsys = new_solar_system()

        mercury = Mercury()#debugger="gdb")
        mercury.initialize_code()
        mercury.parameters.integrator=2

        mercury.particles.add_particles(solsys)
        mercury.commit_particles()

        start_pos = mercury.orbiters[2].position
        mercury.evolve_model(365.14|units.day)
        self.assertAlmostEqual(mercury.orbiters[2].position, start_pos, 1)
        mercury.stop()

    def test21(self):
        solsys = new_solar_system()

        mercury = Mercury()
        mercury.initialize_code()
        
        names=["elements_file","close_encounters_file","info_file",
         "bigbody_file","smallbody_file","integration_parameters_file","restart_file"]
        
        for name in names:
          self.assertEqual(getattr(mercury.parameters,name),"/dev/null")
        
        for name in names:
          setattr(mercury.parameters,name,os.path.join(mercury.output_directory,name))

        for name in names:
          self.assertEqual(getattr(mercury.parameters,name),os.path.join(mercury.output_directory,name))

        mercury.stop()

    def xtest22(self):
        """ collision test hangs or fails if internal collision detection is enabled """
         
        def collision():
          
            M=1.| units.MSun
            m=1.| units.MJupiter
            r=5. | units.AU
            vcirc=(constants.G*(M+m)/r)**0.5
            
            sys=datamodel.Particles(4)
            sys[0].mass=M
            sys[0].radius=1. | units.RSun
            sys[0].x=0 | units.AU
            sys[0].y=0 | units.AU
            sys[0].z=0 | units.AU
            sys[0].vx=0 | units.kms
            sys[0].vy=0 | units.kms
            sys[0].vz=0 | units.kms
            
            sys[1].mass=m
            sys[1].radius=0.01 | units.RSun
            sys[1].x=r
            sys[1].y=0 | units.AU
            sys[1].z=0 | units.AU
            sys[1].vx=0 | units.kms
            sys[1].vy=vcirc
            sys[1].vz=0 | units.kms
            
            sys[2].mass=m
            sys[2].radius=0.01 | units.RSun
            sys[2].x=-r
            sys[2].y=0 | units.AU
            sys[2].z=0 | units.AU
            sys[2].vx=0 | units.kms
            sys[2].vy=vcirc
            sys[2].vz=0 | units.kms
          
            sys[3].mass=m
            sys[3].radius=0.01 | units.RSun
            sys[3].x=0 | units.AU
            sys[3].y=r
            sys[3].z=0 | units.AU
            sys[3].vx=0 | units.kms
            sys[3].vy=0 | units.kms
            sys[3].vz=0 | units.kms
            
            return sys
  
        code=Mercury()
                
        sys=collision()
        
        code.particles.add_particles(sys)
        
        tend=3.5| units.yr
        dt=100. | units.day
        tnow=code.model_time
                
        while tnow<tend:
          code.evolve_model(tnow+dt)
          tnow=code.model_time
          print(tnow.in_(units.yr))

        code.stop()
        
