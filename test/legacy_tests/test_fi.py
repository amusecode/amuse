import os
import sys
import numpy
from amuse.test.amusetest import TestWithMPI

from amuse.legacy.fi.interface import fi, Fi
from amuse.ext.evrard_test import MakeEvrardTest

from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.data import core
from amuse.legacy.support import channel

class TestMPIInterface(TestWithMPI):

    def test1(self):
        instance=fi()  
        instance.setup_module()
        instance.stop()
    
    def test2(self):
        instance=fi()  
        instance.setup_module()
        
        for x,l in [('usesph',0),('radiate',1),('starform',1),('cosmo',1),
                    ('sqrttstp',1),('acc_tstp',0),('freetstp',1),('usequad',1),
                    ('directsum',1),('selfgrav',0),('fixthalo',1),
                    ('adaptive_eps',1),('gdgop',0),('smoothinput',1),
                    ('consph',0),('sphinit',0),('uentropy',1),('isotherm',1),
                    ('eps_is_h',0)]:
            result,err=eval("instance.get_"+x)()
            self.assertEquals( (x,result),(x,l))
            err=eval("instance.set_"+x)(1)
            result,err=eval("instance.get_"+x)()
            self.assertEquals( (x,result),(x,1))
            err=eval("instance.set_"+x)(0)
            result,err=eval("instance.get_"+x)()
            self.assertEquals((x,result),(x,0))
        
        for x,i in [ ('firstsnap',0),('stepout',5),('steplog',5),('max_tbin',4096),
                     ('minppbin',1),('targetnn',32),('verbosity',0),('nsmooth',64)]:
            result,err=eval("instance.get_"+x)()
            self.assertEquals( (x,result),(x,i))
            err=eval("instance.set_"+x)(12345)
            result,err=eval("instance.get_"+x)()
            self.assertEquals((x,result),(x,12345))
        
        for x,r in [ ('pboxsize',300.),('unitl_in_kpc',1.),('unitm_in_msun',1.e9),('dtime',1.),
                     ('tstepcrit',1.),('tstpcr2',0.25),('freev',0.5),('freea',0.35),('freevexp',0.),
                     ('freeaexp',-1.),('bh_tol',0.5),('gdgtol',0.01),('nn_tol',0.1),
                     ('epsgas',0.005),('gamma',1.6666667),('alpha',0.5),('beta',1.0),('epssph',0.01),
                     ('courant',0.3),('removgas',0.25),('consthsm',0.2),('nsmtol',0.1),
                     ('graineff',0.05),('crionrate',3.6),('heat_par1',0.),('heat_par2',0.),
                     ('cool_par',1.),('optdepth',0.),('tcollfac',1.),('masscrit',1.e5),
                     ('sfeff',0.25),('tbubble',3.e7),('sne_eff',0.),('tsnbeg',3.e6),
                     ('rhomax',100.),('eps',0.)]:
            result,err=eval("instance.get_"+x)()
            self.assertAlmostEquals(result,r,7)
            err=eval("instance.set_"+x)(0.)
            result,err=eval("instance.get_"+x)()
            self.assertEquals(result,0.)
            err=eval("instance.set_"+x)(0.12345)
            result,err=eval("instance.get_"+x)()
            self.assertAlmostEquals(result,0.12345,7)
        
        for x,s in [('halofile','none'),('feedback','fuv'),('sfmode','gerritsen'),
                    ('hupdatemethod','mass'),('sph_visc','sph')]:
            result,err=eval("instance.get_"+x)()
            self.assertEquals((x,result),(x,s))
            err=eval("instance.set_"+x)("123")
            result,err=eval("instance.get_"+x)()
            self.assertEquals((x,result),(x,"123"))
        
        instance.stop()
    
    def test3(self):
        instance=fi()  
        instance.setup_module()
        instance.new_particle(11.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        retrieved_state = instance.get_state(1)
        self.assertEquals(11.0,  retrieved_state['mass'])
        self.assertEquals(2.0, retrieved_state['radius'])
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 1)
        instance.cleanup_module()
        instance.stop()
    
    def test4(self):
        instance=fi()  
        instance.setup_module()
        instance.new_particle(11.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        instance.set_time_step(2.0)
        retrieved_state = instance.get_time_step()
        self.assertEquals(2.0, retrieved_state['time_step'])
        instance.cleanup_module()
        instance.stop()
    
    def test5(self):
        instance=fi()
        instance.setup_module()
        instance.set_eps(0.001)
        instance.set_directsum(0)
        instance.commit_parameters()
        instance.new_particle( 
           [1.0,1.0,1.0],
           [0.0,0.0,0.0],
           [1.0,0.0,-1.0],
           [0.0,0.0,0.0],
           [0.0,0.0,0.0],
           [0.0,1.0,0.0],
           [0.0,0.0,0.0],
           [0.0,0.0,0.0] )
        instance.initialize_particles(0.0)
        self.assertEqual(instance.get_number_of_particles().number_of_particles, 3)
        instance.synchronize_model()        
        Ep=instance.get_potential_energy()['potential_energy']
        Ek=instance.get_kinetic_energy()['kinetic_energy']
        
        self.assertAlmostEqual( Ek, 0.5,10)
        self.assertAlmostEqual( Ep, -2.5,10)
        instance.delete_particle(2)
        instance.reinitialize_particles() 
        instance.synchronize_model()        
        n=instance.get_number_of_particles()['number_of_particles']
        Ep=instance.get_potential_energy()['potential_energy']
        Ek=instance.get_kinetic_energy()['kinetic_energy']
        self.assertEqual( n, 2)
        self.assertAlmostEqual( Ek, 0.,10)
        self.assertAlmostEqual( Ep, -0.5,10)        
        
        instance.cleanup_module()
        instance.stop()

class TestEvrard(TestWithMPI):

    def test0(self):
        evrard=MakeEvrardTest(1000)
        mass,x,y,z,vx,vy,vz,u=evrard.new_model()
        smooth=numpy.zeros_like(mass)
        nb = fi()
        nb.setup_module()

        nb.set_stepout(99999)
        nb.set_steplog(99999)
        nb.set_usesph(0)
        nb.set_radiate(1)
        nb.set_dtime(0.05)
        nb.set_gdgop(0)
        nb.set_uentropy(0)
        nb.set_verbosity(0)
        nb.commit_parameters()
        
        ids,error = nb.new_sph_particle(mass,smooth,x,y,z,vx,vy,vz,u)
        if filter(lambda x: x != 0, error) != []: raise Exception
     
        nb.initialize_particles(0.0)

        nb.synchronize_model()
        Ek,ret=nb.get_kinetic_energy()
        Ep,ret=nb.get_potential_energy()
        Eth,ret=nb.get_thermal_energy()

        self.assertAlmostEqual( Ek, 0.,3)
        self.assertAlmostEqual( Ep, -0.6611,3)        
        self.assertAlmostEqual( Eth, 0.05,3)        

        nb.evolve(0.5)
        nb.synchronize_model()
        Ek,ret=nb.get_kinetic_energy()
        Ep,ret=nb.get_potential_energy()
        Eth,ret=nb.get_thermal_energy()

        self.assertAlmostEqual( Ek, 0.129577,3)
        self.assertAlmostEqual( Ep, -0.831976,3)
        self.assertAlmostEqual( Eth,0.08567999,3)

        del evrard
        nb.stop()
        

class TestFiInterface(TestWithMPI):

    def test0(self):
        instance=Fi()
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.parameters.timestep = 0.5 | units.day
        self.assertEquals(instance.parameters.timestep, 0.5 | units.day)
        instance.commit_parameters()
        instance.commit_particles()
        self.assertEquals(instance.parameters.timestep, 0.5 | units.day)
        instance.stop()
    
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        instance = Fi(convert_nbody)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.00000001 | units.AU**2
        instance.parameters.timestep = 0.5 | units.day
        self.assertEquals(instance.parameters.timestep, 0.5 | units.day)
        instance.commit_parameters()
        
        stars = core.Particles(2)
        
        sun = stars[0]
        sun.mass = units.MSun(1.0)
        sun.position = [0.0,0.0,0.0] | units.m
        sun.velocity = [0.0,0.0,0.0] | units.ms
        sun.radius = units.RSun(1.0)
        
        earth = stars[1]
        earth.mass = units.kg(5.9736e24)
        earth.radius = units.km(6371) 
        earth.position = [149.5e6, 0.0, 0.0] | units.km
        earth.velocity = [0.0, 29800, 0.0] | units.ms
        
        instance.particles.add_particles(stars)
        self.assertEquals(instance.parameters.timestep, 0.5 | units.day)
        
        postion_at_start = earth.position.x
        
        instance.evolve_model(365.0 | units.day)
        
        instance.update_particles(stars)
        
        postion_after_full_rotation = earth.position.x
        
        self.assertAlmostRelativeEqual(postion_at_start, postion_after_full_rotation, 4)
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
       
        instance.update_particles(stars)
        
        postion_after_half_a_rotation = earth.position.x
        
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        instance.update_particles(stars)
        
        postion_after_half_a_rotation = earth.position.x
        self.assertAlmostRelativeEqual(-postion_at_start, postion_after_half_a_rotation, 3)
        
        instance.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
         
        instance.update_particles(stars)
        
        postion_after_half_a_rotation = earth.position.y
        
        self.assertAlmostRelativeEqual(-postion_at_start, postion_after_half_a_rotation, 4)
        
        instance.cleanup_code()
        instance.stop()
    
    def test2(self):
        print "Test 2: testing fi data directory"
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        instance = Fi(convert_nbody, redirection='none')
        error = instance.initialize_code()
        self.assertEquals(0, error)
        self.assertTrue('data/fi/input/' in instance.get_fi_data_directory().fi_data_directory)
        
        self.assertEquals(False, instance.get_radiation_flag())
        instance.set_radiation_flag(True)
        self.assertEquals(True, instance.get_radiation_flag())
        
        self.assertEquals(False, instance.get_star_formation_flag())
        instance.set_star_formation_flag(True)
        self.assertEquals(True, instance.get_star_formation_flag())
        instance.commit_parameters()
        
        stars = core.Particles(2)
        stars.mass = [1.0, 3.0e-6] | units.MSun
        stars.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.AU
        stars.velocity = [[0.0,0.0,0.0], [0.0,29.8,0.0]] | units.km / units.s
        stars.radius = [1.0, 0.01] | units.RSun
        
        instance.particles.add_particles(stars)
        instance.commit_particles()
        instance.evolve_model(365.0 | units.day)
        
        instance.cleanup_code()
        instance.stop()
